//! Angle-of-repose benchmark — validates the bulk-friction response of the
//! granular contact model by forming a static heap and measuring its repose
//! angle θ_r against the empirical expectation θ_r(μ).
//!
//! Protocol — lift the cylinder (heap on a frictional floor), with the declared
//! dissipative formation restitution e=0.4:
//!   1. Fill: particles are inserted inside a thin z-aligned cylinder wall that
//!      sits on the floor and settle into a packed column under gravity.
//!   2. Lift: only after the fastest particle is below the declared fill-rest
//!      threshold is the confining cylinder wall deactivated (the "lift").
//!   3. Record: when the heap has come to rest, every particle's final
//!      (x, y, z, radius) is written to data/repose_results.csv. The geometry
//!      fit (height-vs-radius slope → θ_r = atan(slope)) is done in sweep.py.
//!
//! This executable deliberately fails rather than fabricating a snapshot when
//! either rest criterion is not reached in its configured stage budget. A
//! formation run at e=0.4 is not evidence that a value transfers to canonical
//! e=0.926 glass; that separate claim requires a protocol-matched validation.
//!
//! The base friction that keeps the bottom layer from sliding out is supplied by
//! dirt_wall's Mindlin sliding friction on the floor plane (the material's μ via
//! friction_ij) — no frozen particle bed is needed.
//!
//! Cylinder walls have no WallMotion support in dirt_wall (wall_move drives plane
//! walls only), so the "lift" is the runtime deactivation below.
//!
//! main.rs is a thin recorder: it drives the two stages and dumps raw particle
//! positions. All theory/fitting/validation lives in sweep.py.
//!
//! ```bash
//! cargo run --release --example sphcal_angle_of_repose --no-default-features --features precision-double -- examples/SPH_glass_sphere_calibration/03_angle_of_repose/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::dirt_granular::tangential::ContactHistoryStore;
use dirt_core::prelude::*;
use serde::Deserialize;
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
    fill_vmax: Option<f64>,
    fill_rest_samples: usize,
    heap_rest_samples: usize,
    written: bool,
}

/// The pose measured after release must be at least as well settled as the
/// confined state from which it was released.  A looser post-lift velocity
/// threshold can capture a transient, rate-dependent slope rather than a
/// static angle of repose.
const REST_MAX_SPEED_M_S: f64 = 2e-3;
/// Consecutive 200-step samples required before treating a state as static.
/// A one-sample crossing can be a turning point in a slowly oscillating heap.
const REST_DWELL_SAMPLES: usize = 10;

/// Optional diagnostic release condition used only for cross-code replay.
/// LAMMPS data files cannot restore DIRT's granular contact-history payload.
#[derive(Clone, Deserialize)]
#[serde(default)]
struct ReposeReplayConfig {
    zero_contact_history_at_lift: bool,
}

impl Default for ReposeReplayConfig {
    fn default() -> Self {
        Self {
            zero_contact_history_at_lift: false,
        }
    }
}

impl ReposeTracker {
    fn new() -> Self {
        Self {
            lift_step: None,
            fill_vmax: None,
            fill_rest_samples: 0,
            heap_rest_samples: 0,
            written: false,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(FixesPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(StatesPlugin::new(
            Phase::Fill,
            ParticleSimScheduleSet::PostFinalIntegration,
        ))
        .add_plugins(StageAdvancePlugin::<Phase>::new(
            ParticleSimScheduleSet::PostFinalIntegration,
        ));

    app.add_resource(ReposeTracker::new());
    Config::load::<ReposeReplayConfig>(&mut app, "repose_replay");

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
    registry: ResMut<AtomDataRegistry>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    mut walls: ResMut<Walls>,
    replay: Res<ReposeReplayConfig>,
    mut tracker: ResMut<ReposeTracker>,
    mut next_state: ResMut<NextState<Phase>>,
) {
    let step = run_state.total_cycle;
    // This protocol inserts one deterministic, overlap-checked population at
    // setup.  Still require the declared population before lifting: an
    // insertion/setup regression must not turn a quiet partial deposit into a
    // misleading angle-of-repose observation.  Do not use Atom::natoms here:
    // sum owned particles across ranks instead.
    const EXPECTED_HEAP_PARTICLES: usize = 1200;
    let global_particles = -comm.all_reduce_min_f64(-(atoms.nlocal as f64));
    if (global_particles.round() as usize) < EXPECTED_HEAP_PARTICLES {
        return;
    }
    // Give the column time to settle before testing; then test periodically.
    if step < 2000 || step % 200 != 0 {
        return;
    }
    let vmax = max_speed(&atoms, &comm);
    if vmax < REST_MAX_SPEED_M_S {
        tracker.fill_rest_samples += 1;
    } else {
        tracker.fill_rest_samples = 0;
    }
    if tracker.fill_rest_samples >= REST_DWELL_SAMPLES {
        // LAMMPS `read_data` receives particle state but not the pair Mindlin/
        // SDS or wall histories accumulated during formation.  A diagnostic
        // comparison therefore has a valid common boundary only when DIRT
        // explicitly clears those histories too.  The calibration leaves this
        // option false and retains its physical formation history.
        if replay.zero_contact_history_at_lift {
            let mut history =
                registry.expect_mut::<ContactHistoryStore>("zero-history replay lift");
            for contacts in &mut history.contacts {
                contacts.clear();
            }
            walls.tangential_springs.clear();
            walls.rolling_springs.clear();
        }
        // Persist the actual, settled initial state of the collapse before the
        // wall is removed.  A cross-code replay must start from this state;
        // recreating it with another solver's pour/inserter is a different
        // formation experiment, not a contact-model comparison.
        if comm.rank() == 0 {
            let dem = registry.expect::<DemAtom>("lift_when_settled");
            let out_dir = input.output_dir.clone().unwrap_or_else(|| {
                "examples/SPH_glass_sphere_calibration/03_angle_of_repose".to_string()
            });
            let data_dir = format!("{}/data", out_dir);
            fs::create_dir_all(&data_dir).ok();
            let prelift_file = format!("{}/repose_prelift.csv", data_dir);
            let mut f = fs::File::create(&prelift_file)
                .unwrap_or_else(|e| panic!("Cannot create {}: {}", prelift_file, e));
            // Particle state is explicit. Contact-history state is declared in
            // the qualification witness: calibration retains it; the optional
            // LAMMPS diagnostic clears it on both sides.
            writeln!(f, "x,y,z,radius,vx,vy,vz,omega_x,omega_y,omega_z").unwrap();
            for i in 0..atoms.nlocal as usize {
                writeln!(
                    f,
                    "{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
                    atoms.pos[i][0],
                    atoms.pos[i][1],
                    atoms.pos[i][2],
                    dem.radius[i],
                    atoms.vel[i][0],
                    atoms.vel[i][1],
                    atoms.vel[i][2],
                    dem.omega[i][0],
                    dem.omega[i][1],
                    dem.omega[i][2],
                )
                .unwrap();
            }
        }
        walls.deactivate_by_name("cylinder");
        tracker.lift_step = Some(step);
        tracker.fill_vmax = Some(vmax);
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
    replay: Res<ReposeReplayConfig>,
    mut tracker: ResMut<ReposeTracker>,
) {
    if tracker.written {
        return;
    }
    let step = run_state.total_cycle;
    let lift_step = match tracker.lift_step {
        Some(s) => s,
        None => {
            if comm.rank() == 0 {
                eprintln!("ERROR: fill stage ended without meeting the column-rest criterion; refusing to lift or record");
            }
            std::process::exit(2);
        }
    };
    let fill_vmax = tracker.fill_vmax.expect("lift event must record fill vmax");
    // Let the heap slump for a bit, then poll for rest.
    if step < lift_step + 2000 || step % 200 != 0 {
        return;
    }
    // A snapshot is evidence only after the declared heap-rest criterion. Do
    // not turn a runtime cap into a synthetic rest condition.
    let vmax = max_speed(&atoms, &comm);
    if vmax < REST_MAX_SPEED_M_S {
        tracker.heap_rest_samples += 1;
    } else {
        tracker.heap_rest_samples = 0;
    }
    if tracker.heap_rest_samples < REST_DWELL_SAMPLES {
        if step >= lift_step + 150_000 {
            if comm.rank() == 0 {
                eprintln!("ERROR: heap did not meet the rest criterion within 150000 post-lift steps; refusing to record");
            }
            std::process::exit(2);
        }
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
        .unwrap_or_else(|| "examples/SPH_glass_sphere_calibration/03_angle_of_repose".to_string());
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
    // Write a machine-readable witness alongside the snapshot.  The campaign
    // driver validates this record rather than inferring protocol compliance
    // from human-oriented stdout, so copied/truncated CSVs cannot masquerade as
    // qualified settled deposits.
    let qualification_file = format!("{}/repose_qualification.json", data_dir);
    let mut q = fs::File::create(&qualification_file)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", qualification_file, e));
    writeln!(
        q,
        concat!(
            "{{\n  \"schema\": 3,\n  \"history_at_lift\": \"{}\",\n  \"fill_step\": {},\n",
            "  \"fill_vmax_m_s\": {:.8e},\n  \"fill_rest_samples\": {},\n  \"lift_step\": {},\n",
            "  \"heap_step\": {},\n  \"heap_vmax_m_s\": {:.8e},\n  \"heap_rest_samples\": {},\n",
            "  \"particle_count\": {}\n}}"
        ),
        if replay.zero_contact_history_at_lift {
            "cleared"
        } else {
            "retained"
        },
        lift_step,
        fill_vmax,
        tracker.fill_rest_samples,
        lift_step,
        step,
        vmax,
        tracker.heap_rest_samples,
        nlocal
    )
    .unwrap();
    tracker.written = true;

    println!(
        "Step {}: max speed = {:.3e} m/s — heap at rest, wrote {} particles -> {} (qualification -> {})",
        step, vmax, nlocal, results_file, qualification_file
    );

    // This executable's sole product is the qualified settled snapshot above.
    // Returning from this system would let App::start continue the remaining run
    // budget, which both wastes the campaign wall time and leaves the sweep
    // process alive after its evidence has been written.  Exit only after the
    // CSV has been fully written and the single-write guard has been set.
    std::process::exit(0);
}
