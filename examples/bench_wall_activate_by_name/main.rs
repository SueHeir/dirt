//! bench_wall_activate_by_name — runtime named-wall reactivation check.
//!
//! A single sphere is held at a fixed overlap with a named `dirt_wall` plane.
//! The example samples the wall force while the named wall is active, then calls
//! `Walls::deactivate_by_name`, then `Walls::activate_by_name` on the same name.
//! All geometry and material parameters are fixed; the only changing input is the
//! wall active flag. The recorder writes the per-sample particle force and wall
//! force accumulator to `data/wall_activate_by_name_results.csv`.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::dirt_wall::wall_contact_force;
use dirt_core::prelude::*;
use serde::Deserialize;
use std::fs;
use std::io::Write as IoWrite;

#[derive(Deserialize)]
struct ConfigFile {
    material: MaterialCfg,
    scenario: ScenarioCfg,
}

#[derive(Deserialize)]
struct MaterialCfg {
    name: String,
    youngs_mod: f64,
    poisson_ratio: f64,
    restitution: f64,
    friction: f64,
}

#[derive(Clone, Deserialize)]
struct ScenarioCfg {
    wall_name: String,
    radius: f64,
    overlap: f64,
    dt: f64,
    samples_per_phase: usize,
}

#[derive(Clone)]
struct PhaseSample {
    phase: &'static str,
    active: bool,
    expected: &'static str,
}

struct TogglePlan {
    wall_name: String,
    samples: Vec<PhaseSample>,
    idx: usize,
}

struct Recorder {
    file: fs::File,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let config_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("examples/bench_wall_activate_by_name/config.toml");
    let cfg: ConfigFile = toml::from_str(&fs::read_to_string(config_path).expect("read config"))
        .expect("parse config");

    let mut app = App::new();
    app.add_resource(make_atom(
        cfg.scenario.radius,
        cfg.scenario.overlap,
        cfg.scenario.dt,
    ));
    app.add_resource(make_registry(cfg.scenario.radius));
    app.add_resource(make_material_table(&cfg.material));
    app.add_resource(make_walls(&cfg.scenario));

    let samples = make_samples(cfg.scenario.samples_per_phase);
    let nsamples = samples.len();
    app.add_resource(TogglePlan {
        wall_name: cfg.scenario.wall_name.clone(),
        samples,
        idx: 0,
    });
    app.add_resource(Recorder {
        file: create_results_file(),
    });
    app.add_update_system(control_wall_and_zero, ParticleSimScheduleSet::PreForce);
    app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
    app.add_update_system(record_force, ParticleSimScheduleSet::PostForce);
    app.organize_systems();

    for _ in 0..nsamples {
        app.run();
    }

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let max_abs_fz = atom.force[0][2].abs() as f64;
    println!(
        "bench_wall_activate_by_name wrote examples/bench_wall_activate_by_name/data/wall_activate_by_name_results.csv; last |Fz| = {:.6e} N",
        max_abs_fz
    );
}

fn make_samples(samples_per_phase: usize) -> Vec<PhaseSample> {
    let mut samples = Vec::with_capacity(samples_per_phase * 3);
    for _ in 0..samples_per_phase {
        samples.push(PhaseSample {
            phase: "active_before",
            active: true,
            expected: "nonzero",
        });
    }
    for _ in 0..samples_per_phase {
        samples.push(PhaseSample {
            phase: "deactivated",
            active: false,
            expected: "zero",
        });
    }
    for _ in 0..samples_per_phase {
        samples.push(PhaseSample {
            phase: "reactivated",
            active: true,
            expected: "nonzero",
        });
    }
    samples
}

fn make_atom(radius: f64, overlap: f64, dt: f64) -> Atom {
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);

    let mut atom = Atom::new();
    atom.dt = dt;
    atom.push_test_atom(1, [0.0, 0.0, radius - overlap], radius, mass);
    atom.nlocal = 1;
    atom.natoms = 1;
    atom
}

fn make_registry(radius: f64) -> AtomDataRegistry {
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);

    let mut dem = DemAtom::new();
    dem.radius.push(radius);
    dem.density.push(density);
    dem.inv_inertia.push(1.0 / (0.4 * mass * radius * radius));
    dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
    dem.omega.push([0.0; 3]);
    dem.ang_mom.push([0.0; 3]);
    dem.torque.push([0.0; 3]);
    dem.body_id.push(0.0);

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, 1).unwrap();
    registry
}

fn make_material_table(mat: &MaterialCfg) -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material(
        &mat.name,
        mat.youngs_mod,
        mat.poisson_ratio,
        mat.restitution,
        mat.friction,
        0.0,
        0.0,
    );
    mt.build_pair_tables();
    mt
}

fn make_walls(sc: &ScenarioCfg) -> Walls {
    let plane = WallPlane {
        point_x: 0.0,
        point_y: 0.0,
        point_z: 0.0,
        normal_x: 0.0,
        normal_y: 0.0,
        normal_z: 1.0,
        material_index: 0,
        name: Some(sc.wall_name.clone()),
        bound_x_low: f64::NEG_INFINITY,
        bound_x_high: f64::INFINITY,
        bound_y_low: f64::NEG_INFINITY,
        bound_y_high: f64::INFINITY,
        bound_z_low: f64::NEG_INFINITY,
        bound_z_high: f64::INFINITY,
        velocity: [0.0; 3],
        motion: WallMotion::Static,
        origin: [0.0, 0.0, 0.0],
        force_accumulator: 0.0,
        temperature: None,
    };

    Walls {
        planes: vec![plane],
        active: vec![true],
        cylinders: Vec::new(),
        cylinder_active: Vec::new(),
        spheres: Vec::new(),
        sphere_active: Vec::new(),
        regions: Vec::new(),
        region_active: Vec::new(),
        time: 0.0,
        tangential_springs: std::collections::HashMap::new(),
        rolling_springs: std::collections::HashMap::new(),
    }
}

fn create_results_file() -> fs::File {
    let data_dir = "examples/bench_wall_activate_by_name/data";
    fs::create_dir_all(data_dir).expect("create data dir");
    let path = format!("{data_dir}/wall_activate_by_name_results.csv");
    let mut file = fs::File::create(&path).expect("create results csv");
    writeln!(
        file,
        "sample,phase,wall_active,expected_response,particle_fz,wall_force"
    )
    .unwrap();
    file
}

fn control_wall_and_zero(mut atoms: ResMut<Atom>, mut walls: ResMut<Walls>, plan: Res<TogglePlan>) {
    atoms.force[0] = [0.0; 3];
    for wall in &mut walls.planes {
        wall.force_accumulator = 0.0;
    }

    let sample = &plan.samples[plan.idx];
    if sample.active {
        walls.activate_by_name(&plan.wall_name);
    } else {
        walls.deactivate_by_name(&plan.wall_name);
    }
}

fn record_force(
    atoms: Res<Atom>,
    walls: Res<Walls>,
    mut plan: ResMut<TogglePlan>,
    mut recorder: ResMut<Recorder>,
) {
    let sample = &plan.samples[plan.idx];
    writeln!(
        recorder.file,
        "{},{},{},{},{:.12e},{:.12e}",
        plan.idx,
        sample.phase,
        walls.active[0],
        sample.expected,
        atoms.force[0][2],
        walls.planes[0].force_accumulator
    )
    .unwrap();
    plan.idx += 1;
}
