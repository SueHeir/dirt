//! bench_sliding_friction — validates Coulomb (kinetic) sliding friction and the
//! slip→roll transition against the classic rigid-body "ball thrown onto a rough
//! floor" result.
//!
//! A single sphere is launched horizontally with translational speed v0 and ZERO
//! initial spin onto a flat frictional floor under gravity. While the contact
//! point slides (v > ωR), kinetic friction decelerates the center at a constant
//! a = μ g and spins the sphere up at α = μ g / ((2/5) R). Sliding ends at
//!   t* = 2 v0 / (7 μ g)
//! after which the sphere rolls without slipping at
//!   v_final = (5/7) v0      (independent of μ).
//!
//! Floor: a real `dirt_wall` z-plane at z = 0 (normal +z). `dirt_wall` now applies
//! Mindlin tangential (sliding) friction with a Coulomb cap on plane walls, using
//! the material's `friction` coefficient. The wall therefore decelerates the
//! sliding sphere at a = μg and spins it up — no giant-sphere floor hack needed.
//! The floor is perfectly flat, so there is no curvature systematic.
//!
//! This recorder is thin: it writes the raw time series (t, vx, omega_y, contact)
//! to `<output_dir>/data/sliding_friction_results.csv` at PostFinalIntegration.
//! All theory/validation/plotting lives in `sweep.py`.
//!
//! ```bash
//! cargo run --release --example bench_sliding_friction --no-default-features -- examples/bench_sliding_friction/config.toml
//! ```

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

struct SlideTracker {
    /// Open results file handle; header written on the first recorded row.
    file: Option<fs::File>,
}

impl SlideTracker {
    fn new() -> Self {
        Self { file: None }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin) // body force m*g
        .add_plugins(WallPlugin); // frictional z-plane floor

    app.add_resource(SlideTracker::new());
    app.add_update_system(record, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn record(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<SlideTracker>,
) {
    if atoms.nlocal == 0 {
        return;
    }
    let dem = registry.expect::<DemAtom>("record");

    // Single particle at index 0.
    let p = 0usize;

    // Open the results file (with header) on the first recorded step.
    if tracker.file.is_none() {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_sliding_friction".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/sliding_friction_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(f, "t,vx,omega_y,radius,in_contact").unwrap();
        tracker.file = Some(f);
        println!("=== Sliding-friction recorder ===");
        println!("  floor: dirt_wall z-plane at z = 0 (normal +z)");
        println!("  results -> {}", results_file);
    }

    let t = run_state.total_cycle as f64 * atoms.dt;
    let vx = atoms.vel[p][0];
    let omega_y = dem.omega[p][1];
    let r = dem.radius[p];

    // Particle–wall contact: overlap with the z = 0 floor is (R - z) > 0.
    let pz = atoms.pos[p][2];
    let in_contact = pz < r;

    if let Some(f) = tracker.file.as_mut() {
        writeln!(
            f,
            "{:.10e},{:.10e},{:.10e},{:.10e},{}",
            t, vx, omega_y, r, if in_contact { 1 } else { 0 }
        )
        .unwrap();
    }
}
