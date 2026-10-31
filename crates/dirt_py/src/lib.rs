//! `import dirt`: build DIRT simulations from Python.
//!
//! The App API is the one every GRASS module shares (`grass_pyo3`); this
//! module adds DIRT's plugins, its schedule phases, and a `Particles` view
//! so Python systems can read particle state and apply external forces.
//! A `dirt.App` joins a parent from any other GRASS module (for example a
//! `grass.App` running a Python fluid code) through the stable C ABI.

use std::any::{Any, TypeId};

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use dirt_core::soil_core::{Atom, AtomDataRegistry, ParticleSimScheduleSet as P};
use grass_pyo3::{RustPhase, RustPluginPy, RustView};
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::prelude::*;


/// Local particles of a DIRT App, as NumPy copies. `f_ext` (N x 3, zeros on
/// entry) is added to the particle forces after the system returns, so
/// run the system in `dirt.phases.PostForce`.
#[pyclass]
struct Particles {
    #[pyo3(get)]
    x: Py<PyArray2<f64>>,
    #[pyo3(get)]
    v: Py<PyArray2<f64>>,
    #[pyo3(get)]
    radius: Py<PyArray1<f64>>,
    #[pyo3(get)]
    mass: Py<PyArray1<f64>>,
    #[pyo3(get, set)]
    f_ext: Py<PyArray2<f64>>,
}

fn read(py: Python<'_>, deps: &[&dyn Any]) -> PyResult<PyObject> {
    let atoms = deps[0].downcast_ref::<Atom>().expect("Atom");
    let reg = deps[1].downcast_ref::<AtomDataRegistry>().expect("registry");
    let n = atoms.nlocal() as usize;
    let dem = reg.expect::<DemAtom>("dirt.Particles");
    let rows = |c: &[[f64; 3]]| c[..n].iter().map(|p| p.to_vec()).collect::<Vec<_>>();
    let x = PyArray2::from_vec2(py, &rows(&atoms.pos[..]))?;
    let v = PyArray2::from_vec2(py, &rows(&atoms.vel[..]))?;
    let radius = PyArray1::from_slice(py, &dem.radius[..n]);
    let mass: Vec<f64> = (0..n)
        .map(|i| 4.0 / 3.0 * std::f64::consts::PI * dem.radius[i].powi(3) * dem.density[i])
        .collect();
    let p = Particles {
        x: x.unbind(),
        v: v.unbind(),
        radius: radius.unbind(),
        mass: PyArray1::from_vec(py, mass).unbind(),
        f_ext: PyArray2::zeros(py, [n, 3], false).unbind(),
    };
    Ok(Py::new(py, p)?.into_any())
}

fn write(obj: &Bound<'_, PyAny>, deps: &mut [&mut dyn Any]) -> PyResult<()> {
    let p = obj.downcast::<Particles>()?.borrow();
    let f: PyReadonlyArray2<f64> = p.f_ext.bind(obj.py()).readonly();
    let f = f.as_array();
    let atoms = deps[0].downcast_mut::<Atom>().expect("Atom");
    let n = (atoms.nlocal() as usize).min(f.nrows());
    for i in 0..n {
        for d in 0..3 {
            atoms.force[i][d] += f[[i, d]];
        }
    }
    Ok(())
}

macro_rules! plugin {
    ($fname:ident, $label:literal, $doc:literal, $build:expr) => {
        #[doc = $doc]
        #[pyfunction(name = $label)]
        fn $fname() -> RustPluginPy {
            RustPluginPy::new($label, $build)
        }
    };
}

plugin!(core_plugins, "Core", "Input, comm, domain, neighbor, run ([[run]] stages, dt, steps) and thermo print.", |app| {
    app.add_plugins(CorePlugins);
});
plugin!(granular, "GranularDefault", "DEM atom data, insertion, Verlet, Hertz-Mindlin contact, rotation.", |app| {
    app.add_plugins(GranularDefaultPlugins);
});
plugin!(gravity, "Gravity", "[gravity] body force.", |app| {
    app.add_plugins(GravityPlugin);
});
plugin!(walls, "Wall", "[[wall]] planes.", |app| {
    app.add_plugins(WallPlugin);
});

#[pymodule]
fn dirt(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    grass_pyo3::register(m)?;
    let plugins = m.getattr("plugins")?.downcast_into::<PyModule>()?;
    plugins.add_function(wrap_pyfunction!(core_plugins, &plugins)?)?;
    plugins.add_function(wrap_pyfunction!(granular, &plugins)?)?;
    plugins.add_function(wrap_pyfunction!(gravity, &plugins)?)?;
    plugins.add_function(wrap_pyfunction!(walls, &plugins)?)?;

    let phases = PyModule::new(py, "phases")?;
    for (name, p) in [
        ("InitialIntegration", P::InitialIntegration),
        ("PreForce", P::PreForce),
        ("Force", P::Force),
        ("PostForce", P::PostForce),
        ("FinalIntegration", P::FinalIntegration),
        ("PostFinalIntegration", P::PostFinalIntegration),
    ] {
        phases.add(name, Py::new(py, RustPhase::new(p))?)?;
    }
    m.add_submodule(&phases)?;
    py.import("sys")?.getattr("modules")?.set_item("dirt.phases", &phases)?;

    m.add_class::<Particles>()?;
    grass_pyo3::register_view(RustView {
        py_type: py.get_type::<Particles>().unbind(),
        deps: vec![
            (TypeId::of::<Atom>(), "Atom"),
            (TypeId::of::<AtomDataRegistry>(), "AtomDataRegistry"),
        ],
        read,
        write: Some(write),
    });
    Ok(())
}
