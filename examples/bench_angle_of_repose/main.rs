//! Angle-of-repose benchmark — validates the bulk-friction response of the
//! granular contact model by forming a static heap and measuring its repose
//! angle θ_r against the empirical expectation θ_r(μ).
//!
//! Protocol — "lift the cylinder" (heap on a frictional floor):
//!   1. Fill: particles are inserted inside a thin z-aligned cylinder wall that
//!      sits on the floor, and settle into a packed column under gravity.
//!   2. Lift: once the column's kinetic energy drops below a threshold, the
//!      confining cylinder wall is deactivated (the "lift"). The column slumps
//!      outward across the frictional floor wall and relaxes into a conical heap.
//!   3. Record: when the heap has come to rest, every particle's final
//!      (x, y, z, radius) is written to data/repose_results.csv. The geometry
//!      fit (height-vs-radius slope → θ_r = atan(slope)) is done in sweep.py.
//!
//! The base friction that keeps the bottom layer from sliding out is supplied by
//! dirt_wall's Mindlin sliding friction on the floor plane (the material's μ via
//! friction_ij) — no frozen particle bed is needed.
//!
//! main.rs is a thin recorder: it drives the two stages and dumps raw particle
//! positions. All theory/fitting/validation lives in sweep.py.
//!
//! ```bash
//! cargo run --release --example bench_angle_of_repose --no-default-features -- examples/bench_angle_of_repose/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Two-stage protocol: confine-and-settle, then lift-and-relax.
#[derive(Clone, Debug, PartialEq, Default, StageEnum)]
enum Phase {
    #[default]
    #[stage("fill")]
    Fill,
    #[stage("lift")]
    Lift,
}

/// Tracks settling and guards against writing the results file twice.
struct ReposeTracker {
    lift_step: Option<usize>,
    written: bool,
}

impl ReposeTracker {
    fn new() -> Self {
        Self {
            lift_step: None,
            written: false,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(StatesPlugin::new(
            Phase::Fill,
            ParticleSimScheduleSet::PostFinalIntegration,
        ))
        .add_plugins(StageAdvancePlugin::<Phase>::new(
            ParticleSimScheduleSet::PostFinalIntegration,
        ));

    app.add_resource(ReposeTracker::new());

    // Stage 1: settle the column, then lift the cylinder.
    app.add_update_system(
        lift_when_settled.run_if(in_state(Phase::Fill)),
        ParticleSimScheduleSet::PostFinalIntegration,
    );
    // Stage 2: wait for the heap to come to rest, then dump positions.
    app.add_update_system(
        record_when_settled.run_if(in_state(Phase::Lift)),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();
}

/// Maximum particle speed (global, m/s). Unlike a mean, a single still-moving
/// particle keeps this above the rest threshold.
fn max_speed(atoms: &Atom, comm: &CommResource) -> f64 {
    let nlocal = atoms.nlocal as usize;
    let local_max: f64 = (0..nlocal)
        .map(|i| {
            let v = atoms.vel[i];
            (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
        })
        .fold(0.0, f64::max);
    // No all_reduce_max in the comm API; max(x) = -min(-x). Single-rank here, so
    // this is exact and also correct under MPI.
    -comm.all_reduce_min_f64(-local_max)
}

/// While filling: once the column has settled (fastest particle below the rest
/// threshold), deactivate the confining cylinder wall and advance to the lift
/// stage.
fn lift_when_settled(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    mut walls: ResMut<Walls>,
    mut tracker: ResMut<ReposeTracker>,
    mut next_state: ResMut<NextState<Phase>>,
) {
    let step = run_state.total_cycle;
    // Give the column time to settle before testing; then test periodically.
    if step < 2000 || step % 200 != 0 {
        return;
    }
    let vmax = max_speed(&atoms, &comm);
    if vmax < 2e-3 {
        walls.deactivate_by_name("cylinder");
        tracker.lift_step = Some(step);
        next_state.set(Phase::Lift);
        if comm.rank() == 0 {
            println!(
                "Step {}: max speed = {:.3e} m/s — column settled, lifting cylinder",
                step, vmax
            );
        }
    }
}

/// After the lift: once the heap has come to rest, dump every particle's final
/// (x, y, z, radius) to data/repose_results.csv exactly once.
fn record_when_settled(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    mut tracker: ResMut<ReposeTracker>,
) {
    if tracker.written {
        return;
    }
    let step = run_state.total_cycle;
    let lift_step = match tracker.lift_step {
        Some(s) => s,
        None => return,
    };
    // Let the heap slump for a bit, then poll for rest.
    if step < lift_step + 2000 || step % 200 != 0 {
        return;
    }
    // The heap geometry locks in long before the last micro-jittering particle
    // stops, so record once the fastest particle is slow (heap surface static at
    // < 1 cm/s) OR after a hard cap of 150k steps post-lift, whichever comes
    // first. This bounds the per-case runtime while the angle is unchanged.
    let vmax = max_speed(&atoms, &comm);
    let timed_out = step >= lift_step + 150_000;
    if vmax >= 1e-2 && !timed_out {
        return;
    }

    // Heap is at rest. Rank 0 writes the results file. (Single-rank example;
    // gather is not needed for the default config, but guard nlocal anyway.)
    if comm.rank() != 0 {
        tracker.written = true;
        return;
    }

    let dem = registry.expect::<DemAtom>("record_when_settled");
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_angle_of_repose".to_string());
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    let results_file = format!("{}/repose_results.csv", data_dir);
    let mut f = fs::File::create(&results_file)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
    // Single material: every particle is a heap particle (no frozen bed to
    // filter out), so sweep.py fits θ_r on all recorded positions.
    writeln!(f, "x,y,z,radius").unwrap();
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal {
        writeln!(
            f,
            "{:.8e},{:.8e},{:.8e},{:.8e}",
            atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2], dem.radius[i]
        )
        .unwrap();
    }
    tracker.written = true;

    println!(
        "Step {}: max speed = {:.3e} m/s — heap at rest, wrote {} particles -> {}",
        step, vmax, nlocal, results_file
    );
}
