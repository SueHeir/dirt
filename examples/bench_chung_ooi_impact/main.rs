//! Chung & Ooi (2011) elastic normal-impact benchmark — Test Cases 1 & 2.
//!
//! Reproduces two of the standard DEM code-verification cases from
//!   L. Chung & J.Y. Ooi, "Benchmark tests for verifying discrete element
//!   modelling codes at particle impact level", Granular Matter 13:643-656 (2011).
//!
//!   * Test 1 — elastic normal impact of two identical spheres.
//!   * Test 2 — elastic normal impact of a sphere with a rigid wall.
//!
//! For both, the contact is perfectly elastic (restitution = 1, no damping), so
//! the exact reference is the Hertz analytical solution for maximum contact
//! force, contact duration, and maximum overlap. This binary runs a single
//! impact and records those three measured quantities; the sweep driver
//! (`sweep.py`) generates the per-velocity configs and gates the measured
//! values against the Hertz theory (the paper's own reference for these cases).
//!
//! ```bash
//! cargo run --release --example bench_chung_ooi_impact --no-default-features \
//!     --features precision-double -- examples/bench_chung_ooi_impact/config.toml
//! ```
//!
//! Both cases are handled by the same tracker: with two local particles it
//! measures the sphere–sphere contact (pair 0–1); with one, it measures the
//! sphere–wall contact against the floor plane at z = 0 (normal +z). In both
//! cases the net force on particle 0 (no gravity, no other contacts) is the
//! contact force, so its peak magnitude is the maximum contact force.

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

/// Tracks the single collision and records the Chung & Ooi verification metrics.
struct ImpactTracker {
    /// True once the pair (or sphere–wall) has entered contact.
    was_in_contact: bool,
    /// True once contact has ended (measurement complete).
    finished: bool,
    /// Relative normal approach speed captured just before first contact (m/s).
    v_impact: f64,
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

impl ImpactTracker {
    fn new() -> Self {
        Self {
            was_in_contact: false,
            finished: false,
            v_impact: 0.0,
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
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(WallPlugin);

    app.add_resource(ImpactTracker::new());

    // Read force AFTER the force phase (forces are zeroed at PostInitialIntegration
    // and freshly computed in Force, so they are valid through PostFinalIntegration).
    app.add_update_system(track_impact, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Overlap and relative normal approach speed for the current configuration.
///
/// Returns `(overlap, v_rel_normal, contact_force_mag)` where `v_rel_normal` is
/// positive while the bodies approach along the contact normal.
fn contact_state(atoms: &Atom, dem: &DemAtom) -> (f64, f64, f64) {
    if atoms.nlocal >= 2 {
        // Sphere–sphere (Test 1): pair 0–1.
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
    } else {
        // Sphere–wall (Test 2): floor plane z = 0, normal +z.
        let z = atoms.pos[0][2] as f64;
        let vz = atoms.vel[0][2] as f64;
        let overlap = dem.radius[0] - z;
        let v_rel = -vz; // positive while descending toward the wall
        let f = atoms.force[0];
        let fmag = ((f[0] as f64).powi(2) + (f[1] as f64).powi(2) + (f[2] as f64).powi(2)).sqrt();
        (overlap, v_rel, fmag)
    }
}

/// Monitors the collision, capturing impact speed, peak contact force, contact
/// duration, and peak overlap, then writes a one-row CSV on separation.
fn track_impact(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<ImpactTracker>,
) {
    if tracker.finished || atoms.nlocal == 0 {
        return;
    }

    let dem = registry.expect::<DemAtom>("track_impact");
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
        // Separation — measurement complete.
        tracker.finished = true;
        tracker.step_contact_end = step;

        let dt = atoms.dt;
        let contact_steps = tracker.step_contact_end - tracker.step_contact_start;
        let contact_time = contact_steps as f64 * dt;

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_chung_ooi_impact".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();

        let results_file = format!("{}/data/impact_results.csv", out_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "n_particles,v_impact,max_force,contact_time,max_overlap,dt,radius,density"
        )
        .unwrap();
        writeln!(
            f,
            "{},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
            atoms.nlocal,
            tracker.v_impact,
            tracker.max_force,
            contact_time,
            tracker.max_overlap,
            dt,
            dem.radius[0],
            dem.density[0],
        )
        .unwrap();

        println!("=== Chung & Ooi (2011) Impact Results ===");
        println!("  Particles:        {}", atoms.nlocal);
        println!("  Impact velocity:  {:.6e} m/s (relative, normal)", tracker.v_impact);
        println!("  Max contact force:{:.6e} N", tracker.max_force);
        println!("  Contact duration: {:.6e} s ({} steps)", contact_time, contact_steps);
        println!("  Peak overlap:     {:.6e} m", tracker.max_overlap);
        println!("  Timestep dt:      {:.6e} s", dt);
        println!("  Results saved to: {}", results_file);
    }

    if !tracker.finished && !in_contact {
        tracker.prev_v_rel = v_rel;
    }
}
