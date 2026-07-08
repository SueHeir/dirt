//! Hooke wall rebound benchmark.
//!
//! Launches one sphere at a rigid `dirt_wall` plane with
//! `contact_model = "hooke"` and records the wall rebound coefficient of
//! restitution, contact duration, peak overlap, and peak wall force.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Tracks the single wall contact and writes one CSV row after separation.
struct WallReboundTracker {
    /// True once the particle has entered contact.
    was_in_contact: bool,
    /// True once contact has ended and the result has been written.
    finished: bool,
    /// Normal impact speed captured just before contact.
    v_impact: f64,
    /// Normal rebound speed captured just after separation.
    v_rebound: f64,
    /// Largest contact force magnitude seen during contact.
    max_force: f64,
    /// Largest geometric overlap with the wall.
    max_overlap: f64,
    /// Timestep index at first contact.
    step_contact_start: usize,
    /// Timestep index at separation.
    step_contact_end: usize,
    /// Previous free-flight normal velocity.
    prev_vz: f64,
}

impl WallReboundTracker {
    fn new() -> Self {
        Self {
            was_in_contact: false,
            finished: false,
            v_impact: 0.0,
            v_rebound: 0.0,
            max_force: 0.0,
            max_overlap: 0.0,
            step_contact_start: 0,
            step_contact_end: 0,
            prev_vz: 0.0,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(WallPlugin);

    app.add_resource(WallReboundTracker::new());
    app.add_update_system(track_rebound, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

fn track_rebound(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<WallReboundTracker>,
) {
    if tracker.finished || atoms.nlocal == 0 {
        return;
    }

    let dem = registry.expect::<DemAtom>("track_rebound");
    let step = run_state.total_cycle;
    let z = atoms.pos[0][2] as f64;
    let vz = atoms.vel[0][2] as f64;
    let r = dem.radius[0];
    let overlap = r - z;
    let in_contact = overlap > 0.0;
    let f = atoms.force[0];
    let fmag = ((f[0] as f64).powi(2) + (f[1] as f64).powi(2) + (f[2] as f64).powi(2)).sqrt();

    if !tracker.was_in_contact && !in_contact {
        tracker.prev_vz = vz;
    } else if !tracker.was_in_contact && in_contact {
        tracker.was_in_contact = true;
        tracker.v_impact = tracker.prev_vz.abs();
        tracker.step_contact_start = step;
        tracker.max_overlap = overlap;
        tracker.max_force = fmag;
    } else if tracker.was_in_contact && in_contact {
        tracker.max_overlap = tracker.max_overlap.max(overlap);
        tracker.max_force = tracker.max_force.max(fmag);
    } else if tracker.was_in_contact && !in_contact {
        tracker.finished = true;
        tracker.step_contact_end = step;
        tracker.v_rebound = vz.abs();

        let dt = atoms.dt;
        let contact_steps = tracker.step_contact_end - tracker.step_contact_start;
        let contact_time = contact_steps as f64 * dt;
        let cor = if tracker.v_impact > 0.0 {
            tracker.v_rebound / tracker.v_impact
        } else {
            0.0
        };

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_hooke_wall_rebound".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();

        let results_file = format!("{}/data/rebound_results.csv", out_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "v_impact,v_rebound,cor_measured,max_force,contact_time,max_overlap,dt,radius,density"
        )
        .unwrap();
        writeln!(
            f,
            "{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
            tracker.v_impact,
            tracker.v_rebound,
            cor,
            tracker.max_force,
            contact_time,
            tracker.max_overlap,
            dt,
            r,
            dem.density[0],
        )
        .unwrap();

        println!("=== Hooke Wall Rebound Results ===");
        println!("  Impact velocity:  {:.6e} m/s", tracker.v_impact);
        println!("  Rebound velocity: {:.6e} m/s", tracker.v_rebound);
        println!("  COR (measured):   {:.6}", cor);
        println!("  Max wall force:   {:.6e} N", tracker.max_force);
        println!(
            "  Contact duration: {:.6e} s ({} steps)",
            contact_time, contact_steps
        );
        println!("  Peak overlap:     {:.6e} m", tracker.max_overlap);
        println!("  Results saved to: {}", results_file);
    }

    if !tracker.finished && !in_contact {
        tracker.prev_vz = vz;
    }
}
