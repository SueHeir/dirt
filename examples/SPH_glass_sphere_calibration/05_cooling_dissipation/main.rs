//! sphcal_cooling_dissipation — granular-temperature dissipation closure.
//!
//! A freely cooling granular gas of glass beads in a periodic box. With DIRT's
//! velocity-independent restitution (constant `e`), the granular temperature
//! obeys Haff's law  T(t) = T0 / (1 + t/tc)^2 ; the cooling time `tc` encodes
//! the inelastic dissipation rate (the de-fluidization energy-balance closure).
//!
//! Thin recorder only: it measures translational and rotational granular
//! temperature and writes raw rows to `<output_dir>/data/cooling.csv`. All
//! theory (the Haff fit and the dissipation coefficient) lives in `sweep.py`.
//!
//! ```bash
//! cargo run --release --example sphcal_cooling_dissipation --no-default-features \
//!     -- examples/SPH_glass_sphere_calibration/05_cooling_dissipation/config.toml
//! ```

use std::fs;
use std::io::Write;
use std::sync::Once;

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins);

    app.add_update_system(record, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Measure translational + rotational granular temperature and append a raw row.
fn record(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
) {
    let step = run_state.total_cycle;
    if step % 2000 != 0 {
        return;
    }

    let nlocal = atoms.nlocal as usize;
    if nlocal == 0 {
        return;
    }

    // ── Translational temperature: T_trans = <m (v - <v>)^2> / (3 M) ──
    let mut m_total = 0.0_f64;
    let mut mv = [0.0_f64; 3];
    for i in 0..nlocal {
        let m = atoms.mass[i];
        m_total += m;
        for d in 0..3 {
            mv[d] += m * atoms.vel[i][d];
        }
    }
    let m_total = comm.all_reduce_sum_f64(m_total);
    let mv = [
        comm.all_reduce_sum_f64(mv[0]),
        comm.all_reduce_sum_f64(mv[1]),
        comm.all_reduce_sum_f64(mv[2]),
    ];
    let v_mean = if m_total > 0.0 {
        [mv[0] / m_total, mv[1] / m_total, mv[2] / m_total]
    } else {
        [0.0; 3]
    };

    let mut ke_trans = 0.0_f64;
    for i in 0..nlocal {
        let m = atoms.mass[i];
        let mut dv2 = 0.0;
        for d in 0..3 {
            let dv = atoms.vel[i][d] - v_mean[d];
            dv2 += dv * dv;
        }
        ke_trans += m * dv2;
    }
    let ke_trans = comm.all_reduce_sum_f64(ke_trans);

    // ── Rotational temperature: T_rot = sum(I w^2) / (3 M),  I = (2/5) m r^2 ──
    let dem = registry.get::<DemAtom>();
    let mut ke_rot = 0.0_f64;
    if let Some(ref dem) = dem {
        for i in 0..nlocal {
            let m = atoms.mass[i];
            let r = dem.radius[i];
            let inertia = 0.4 * m * r * r;
            let w = dem.omega[i];
            ke_rot += inertia * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
        }
    }
    let ke_rot = comm.all_reduce_sum_f64(ke_rot);

    let t_trans = if m_total > 0.0 { ke_trans / (3.0 * m_total) } else { 0.0 };
    let t_rot = if m_total > 0.0 { ke_rot / (3.0 * m_total) } else { 0.0 };
    let t_total = t_trans + t_rot;

    if comm.rank() != 0 {
        return;
    }

    let time = step as f64 * atoms.dt;

    let base = match input.output_dir.as_deref() {
        Some(dir) => dir.to_string(),
        None => "examples/SPH_glass_sphere_calibration/05_cooling_dissipation".to_string(),
    };
    let data_dir = format!("{}/data", base);
    let filepath = format!("{}/cooling.csv", data_dir);

    static INIT: Once = Once::new();
    INIT.call_once(|| {
        fs::create_dir_all(&data_dir).ok();
        let mut f = fs::File::create(&filepath).expect("Cannot create cooling.csv");
        writeln!(f, "step,time,T_trans,T_rot,T_total").unwrap();
    });

    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&filepath)
        .expect("Cannot open cooling.csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e}",
        step, time, t_trans, t_rot, t_total
    )
    .unwrap();
}
