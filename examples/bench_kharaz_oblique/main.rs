//! bench_kharaz_oblique — replicates the oblique-impact experiment of Kharaz,
//! Gorham & Salman, "An experimental study of the elastic rebound of spheres",
//! Powder Technology 120 (2001) 281-291.
//!
//! A single alumina sphere strikes a **flat, frictional glass anvil** (a real
//! `dirt_wall` z-plane, normal +z) obliquely at a fixed impact speed and a chosen
//! angle of incidence. Unlike a frozen *sphere* partner, a flat wall keeps the
//! contact normal exactly +z throughout the collision at every incidence angle —
//! this is the geometry Kharaz used (a thick glass block) and is what lets the
//! normal restitution stay angle-independent up to grazing incidence.
//!
//! The projectile is launched with velocity (v_t, 0, -v_n); the wall is infinite
//! in x/y, so no aiming is needed. This recorder captures the pre-contact
//! (impact) and post-separation (rebound) translational velocity and the rebound
//! spin, and writes them to `<output_dir>/data/oblique_results.csv` in the SAME
//! schema as `bench_oblique_impact`, so the Kharaz sweep driver (`sweep.py`) reads
//! either example's output unchanged. All theory/validation/plotting is in
//! `sweep.py`.
//!
//! ```bash
//! cargo run --release --example bench_kharaz_oblique --no-default-features \
//!     --features precision-double -- examples/bench_kharaz_oblique/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Records the oblique impact of the single projectile against the flat wall.
/// Normal is fixed +z, tangent +x, so no impact-frame rotation is needed.
struct KharazTracker {
    was_in_contact: bool,
    finished: bool,
    vn_impact: f64,
    vt_impact: f64,
    vn_rebound: f64,
    vt_rebound: f64,
    omega_y_rebound: f64,
    step_contact_start: usize,
    step_contact_end: usize,
    max_overlap: f64,
    prev_vel: [f64; 3],
}

impl KharazTracker {
    fn new() -> Self {
        Self {
            was_in_contact: false,
            finished: false,
            vn_impact: 0.0,
            vt_impact: 0.0,
            vn_rebound: 0.0,
            vt_rebound: 0.0,
            omega_y_rebound: 0.0,
            step_contact_start: 0,
            step_contact_end: 0,
            max_overlap: 0.0,
            prev_vel: [0.0; 3],
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(WallPlugin); // frictional Mindlin z-plane anvil (no gravity)

    app.add_resource(KharazTracker::new());
    app.add_update_system(track_kharaz, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn track_kharaz(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<KharazTracker>,
) {
    if tracker.finished || atoms.nlocal < 1 {
        return;
    }
    let dem = registry.expect::<DemAtom>("track_kharaz");
    let step = run_state.total_cycle;

    // Single projectile atom; flat wall sits at z = 0 with normal +z.
    let p = 0usize;
    let vel = [
        atoms.vel[p][0] as f64,
        atoms.vel[p][1] as f64,
        atoms.vel[p][2] as f64,
    ];
    let z = atoms.pos[p][2] as f64;
    let overlap = dem.radius[p] - z; // wall at z = 0
    let in_contact = overlap > 0.0;

    if !tracker.was_in_contact && !in_contact {
        tracker.prev_vel = vel;
    } else if !tracker.was_in_contact && in_contact {
        // First contact: impact frame is fixed (n = +z, t = +x). Use the
        // pre-contact velocity for the incident components.
        tracker.was_in_contact = true;
        let v = tracker.prev_vel;
        tracker.vn_impact = -v[2]; // closing speed (positive when descending)
        tracker.vt_impact = (v[0] * v[0] + v[1] * v[1]).sqrt();
        tracker.step_contact_start = step;
        tracker.max_overlap = overlap;
    } else if tracker.was_in_contact && in_contact {
        if overlap > tracker.max_overlap {
            tracker.max_overlap = overlap;
        }
    } else if tracker.was_in_contact && !in_contact {
        // Separation: record rebound components along the fixed frame.
        tracker.finished = true;
        tracker.vn_rebound = vel[2]; // +z (moving away from the wall)
        tracker.vt_rebound = vel[0]; // +x tangent
        tracker.omega_y_rebound = dem.omega[p][1];
        tracker.step_contact_end = step;

        let dt = atoms.dt;
        let contact_steps = tracker.step_contact_end - tracker.step_contact_start;
        let contact_time = contact_steps as f64 * dt;

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_kharaz_oblique".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/oblique_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "vn_impact,vt_impact,vn_rebound,vt_rebound,omega_y_rebound,contact_time,max_overlap,dt,radius,density"
        )
        .unwrap();
        writeln!(
            f,
            "{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
            tracker.vn_impact,
            tracker.vt_impact,
            tracker.vn_rebound,
            tracker.vt_rebound,
            tracker.omega_y_rebound,
            contact_time,
            tracker.max_overlap,
            dt,
            dem.radius[p],
            dem.density[p],
        )
        .unwrap();

        println!("=== Kharaz oblique impact (flat anvil, n = +z) ===");
        println!(
            "  v_n impact:  {:.6} m/s   v_t impact: {:.6} m/s",
            tracker.vn_impact, tracker.vt_impact
        );
        println!(
            "  v_n rebound: {:.6} m/s   v_t rebound:{:.6} m/s",
            tracker.vn_rebound, tracker.vt_rebound
        );
        println!("  omega_y:     {:.6} rad/s", tracker.omega_y_rebound);
        println!(
            "  contact time:{:.6e} s ({} steps)",
            contact_time, contact_steps
        );
        println!("  results -> {}", results_file);
    }

    if !tracker.finished && !in_contact {
        tracker.prev_vel = vel;
    }
}
