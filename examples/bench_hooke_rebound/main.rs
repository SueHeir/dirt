//! Hooke (linear-spring) normal-contact rebound benchmark.
//!
//! Exercises DIRT's linear spring-dashpot normal contact
//! (`contact_model = "hooke"`, per-material `kn`/`kt`) — the code path in
//! `dirt_granular/src/contact.rs` that every other benchmark leaves untouched
//! (they all use the nonlinear Hertz model).
//!
//! Two identical spheres are launched head-on along x. They collide once and
//! rebound. The binary records the relative normal **impact** and **rebound**
//! speeds (their ratio is the coefficient of restitution), the **contact
//! duration**, the **peak overlap**, and the **peak contact force**, then writes
//! a one-row CSV on separation.
//!
//! The linear spring-dashpot is the one contact law with an *exact* closed-form
//! collision (a constant-coefficient damped harmonic oscillator), so the sweep
//! driver (`sweep.py`) gates these measured quantities against analytical theory
//! — not against DIRT's own output and not against another code. See `sweep.py`
//! for the reference formulae and their derivation.
//!
//! ```bash
//! cargo run --release --example bench_hooke_rebound --no-default-features \
//!     --features precision-double -- examples/bench_hooke_rebound/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Tracks the single head-on collision and records the rebound metrics.
struct ReboundTracker {
    /// True once the pair has entered contact.
    was_in_contact: bool,
    /// True once contact has ended (measurement complete).
    finished: bool,
    /// Relative normal approach speed captured just before first contact (m/s).
    v_impact: f64,
    /// Relative normal separation speed captured at separation (m/s).
    v_rebound: f64,
    /// Peak magnitude of the contact (normal) force during contact (N).
    max_force: f64,
    /// Maximum geometric overlap during contact (m).
    max_overlap: f64,
    /// Timestep index when contact first occurs.
    step_contact_start: usize,
    /// Timestep index when contact ends (separation).
    step_contact_end: usize,
    /// Relative normal velocity at the previous pre-contact step (m/s, +approach).
    prev_v_rel: f64,
}

impl ReboundTracker {
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
            prev_v_rel: 0.0,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins);

    app.add_resource(ReboundTracker::new());

    // Read force AFTER the force phase (forces are zeroed at PostInitialIntegration
    // and freshly computed in Force, so they are valid through PostFinalIntegration).
    app.add_update_system(track_rebound, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Overlap, relative normal velocity, and contact-force magnitude for the pair.
///
/// Returns `(overlap, v_rel_normal, contact_force_mag)` where `v_rel_normal` is
/// **positive while the spheres approach** and negative while they separate.
fn contact_state(atoms: &Atom, dem: &DemAtom) -> (f64, f64, f64) {
    // Sphere–sphere: pair 0–1. The measurement is symmetric in the pair
    // (distance, |v_rel|, |force|), so it does not depend on which atom is 0/1.
    let p0 = atoms.pos[0];
    let p1 = atoms.pos[1];
    let dx = (p0[0] - p1[0]) as f64;
    let dy = (p0[1] - p1[1]) as f64;
    let dz = (p0[2] - p1[2]) as f64;
    let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-30);
    let (nx, ny, nz) = (dx / dist, dy / dist, dz / dist);
    let overlap = dem.radius[0] + dem.radius[1] - dist;
    // Relative velocity of 0 w.r.t. 1, projected on the line of centres.
    // Positive = approaching (0 moving toward 1).
    let vrx = (atoms.vel[0][0] - atoms.vel[1][0]) as f64;
    let vry = (atoms.vel[0][1] - atoms.vel[1][1]) as f64;
    let vrz = (atoms.vel[0][2] - atoms.vel[1][2]) as f64;
    let v_rel = -(vrx * nx + vry * ny + vrz * nz);
    let f = atoms.force[0];
    let fmag = ((f[0] as f64).powi(2) + (f[1] as f64).powi(2) + (f[2] as f64).powi(2)).sqrt();
    (overlap, v_rel, fmag)
}

/// Monitors the collision, capturing impact speed, rebound speed, peak contact
/// force, contact duration, and peak overlap, then writes a one-row CSV.
fn track_rebound(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<ReboundTracker>,
) {
    if tracker.finished || atoms.nlocal < 2 {
        return;
    }

    let dem = registry.expect::<DemAtom>("track_rebound");
    let step = run_state.total_cycle;
    let (overlap, v_rel, fmag) = contact_state(&atoms, &dem);
    let in_contact = overlap > 0.0;

    if !tracker.was_in_contact && !in_contact {
        // Pre-contact: remember the approach speed for the impact-velocity capture.
        tracker.prev_v_rel = v_rel;
    } else if !tracker.was_in_contact && in_contact {
        // First contact.
        tracker.was_in_contact = true;
        tracker.v_impact = tracker.prev_v_rel;
        tracker.step_contact_start = step;
        tracker.max_overlap = overlap;
        tracker.max_force = fmag;
    } else if tracker.was_in_contact && in_contact {
        // During contact: track peaks.
        if overlap > tracker.max_overlap {
            tracker.max_overlap = overlap;
        }
        if fmag > tracker.max_force {
            tracker.max_force = fmag;
        }
    } else if tracker.was_in_contact && !in_contact {
        // Separation — measurement complete. The spheres now recede at constant
        // velocity (no force), so the current relative speed is the rebound speed.
        tracker.finished = true;
        tracker.step_contact_end = step;
        tracker.v_rebound = v_rel.abs();

        let dt = atoms.dt;
        let contact_steps = tracker.step_contact_end - tracker.step_contact_start;
        let contact_time = contact_steps as f64 * dt;
        let cor = if tracker.v_impact != 0.0 {
            tracker.v_rebound / tracker.v_impact
        } else {
            0.0
        };

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_hooke_rebound".to_string());
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
            dem.radius[0],
            dem.density[0],
        )
        .unwrap();

        println!("=== Hooke (linear-spring) Rebound Results ===");
        println!(
            "  Impact velocity:  {:.6e} m/s (relative, normal)",
            tracker.v_impact
        );
        println!(
            "  Rebound velocity: {:.6e} m/s (relative, normal)",
            tracker.v_rebound
        );
        println!("  COR (measured):   {:.6}", cor);
        println!("  Max contact force:{:.6e} N", tracker.max_force);
        println!(
            "  Contact duration: {:.6e} s ({} steps)",
            contact_time, contact_steps
        );
        println!("  Peak overlap:     {:.6e} m", tracker.max_overlap);
        println!("  Timestep dt:      {:.6e} s", dt);
        println!("  Results saved to: {}", results_file);
    }

    if !tracker.finished && !in_contact {
        tracker.prev_v_rel = v_rel;
    }
}
