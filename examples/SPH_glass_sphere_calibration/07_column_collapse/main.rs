//! sphcal_column_collapse — MACRO VALIDATION gate for the SPH glass-sphere
//! calibration: granular column-collapse runout vs aspect ratio, checked against
//! the experimental scaling laws of Lube et al. (2004) and Lajeunesse et al. (2004).
//!
//! A quasi-2D rectangular column of the canonical glass grains (initial width L0,
//! height H) is held against a vertical gate wall on a flat floor. Stage 1
//! ("settle") lets the loosely-inserted column pack down under gravity into a
//! static column. Stage 2 ("collapse") removes the gate on its first step; the
//! column collapses and spreads along +x until it comes to rest. The dimensionless
//! runout (L_f - L0)/L0 is expected to follow:
//!   (L_f - L0)/L0 ~ 1.2 a        (a <~ 2-3, linear regime)
//!   (L_f - L0)/L0 ~ 1.6 a^(2/3)  (a >~ 3,   power-law regime)
//! with a = H/L0.
//!
//! This is the end-to-end macro check the calibrated material (and later the SPH
//! closure) must reproduce — it is NOT a fitted closure parameter. It depends on
//! the calibrated rolling friction mu_r from the 03_angle_of_repose deliverable:
//! smooth spheres without rolling resistance over-run. The shipped material uses a
//! PLACEHOLDER rolling_friction = 0.05 until 03 reports its calibrated value.
//!
//! This recorder is analysis-free: it dumps the final (x, y, z, radius) of every
//! particle to `<output_dir>/data/column_collapse_results.csv`. All runout
//! extraction, regime fitting and PASS/FAIL live in `sweep.py`.
//!
//! ```bash
//! cargo run --release --example sphcal_column_collapse --no-default-features -- examples/SPH_glass_sphere_calibration/07_column_collapse/config.toml
//! ```

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

/// Name of the removable vertical gate wall (matches `name = "gate"` in config).
const GATE_NAME: &str = "gate";

/// Tracks gate release so it happens exactly once.
struct CollapseTracker {
    gate_opened: bool,
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin);

    app.add_resource(CollapseTracker { gate_opened: false });

    // Remove the gate on the first step of the "collapse" stage.
    app.add_update_system(
        open_gate.run_if(in_stage("collapse")),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();

    // Dump the final deposit once the run has finished and the bed is at rest.
    dump_deposit(&app);
}

/// Deactivate the vertical gate wall on the first "collapse" step, releasing the
/// column. Static support removal — no per-particle contact state to reset.
fn open_gate(
    mut tracker: ResMut<CollapseTracker>,
    mut walls: ResMut<Walls>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
) {
    if tracker.gate_opened {
        return;
    }
    walls.deactivate_by_name(GATE_NAME);
    tracker.gate_opened = true;
    if comm.rank() == 0 {
        println!(
            "Step {}: gate removed — column released.",
            run_state.total_cycle
        );
    }
}

/// Write the deposit profile (per-particle x, y, z, radius) so `sweep.py` can
/// extract the final runout L_f. Called after `start()`, so positions are the
/// settled rest state.
fn dump_deposit(app: &App) {
    let atoms = match app.get_resource_ref::<Atom>() {
        Some(a) => a,
        None => return,
    };
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("AtomDataRegistry must exist");
    let dem = registry.expect::<DemAtom>("dump_deposit");
    let nlocal = atoms.nlocal as usize;

    let out_dir = app
        .get_resource_ref::<Input>()
        .and_then(|i| i.output_dir.clone())
        .unwrap_or_else(|| {
            "examples/SPH_glass_sphere_calibration/07_column_collapse".to_string()
        });
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    let results_file = format!("{}/column_collapse_results.csv", data_dir);

    let mut f = fs::File::create(&results_file)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
    writeln!(f, "x,y,z,radius").unwrap();
    let mut vmax = 0.0f64;
    for i in 0..nlocal {
        writeln!(
            f,
            "{:.10e},{:.10e},{:.10e},{:.10e}",
            atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2], dem.radius[i]
        )
        .unwrap();
        let v = atoms.vel[i];
        let s = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if s > vmax {
            vmax = s;
        }
    }

    println!(
        "FINAL: {} particles dumped -> {} (max |v| = {:.3e} m/s)",
        nlocal, results_file, vmax
    );
}
