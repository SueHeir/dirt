//! bench_plate_sinkage — terramechanics pressure–sinkage validation.
//!
//! A settled granular bed is pressed by a flat plate driven vertically downward
//! at a constant slow velocity. The plate is a clipped, downward-facing plane
//! wall (`normal_z = -1`, finite footprint via `bound_x_*`) moving at constant
//! velocity; `dirt_wall` already accumulates the contact (reaction) force on each
//! plane wall in `WallPlane::force_accumulator`, so the vertical load on the
//! plate is read directly — no core change is needed. The recorder streams the
//! plate depth (sinkage `z`) and the vertical reaction force `F` versus time to
//! `<output_dir>/data/plate_sinkage_results.csv`.
//!
//! The empirical reference is the Bekker pressure–sinkage relation
//!   p = (k_c/b + k_φ) · z^n ,
//! i.e. pressure under the plate grows as a power law in sinkage. Validation
//! (in `sweep.py`) checks that p(z) is monotone and well-fit by p ∝ z^n with a
//! physically sensible exponent, and that wider/deeper trends are sane.
//!
//! ```bash
//! cargo run --release --example bench_plate_sinkage --no-default-features -- examples/bench_plate_sinkage/config.toml
//! ```

use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Streams the plate's depth and the vertical reaction force on it each step,
/// starting once the plate first touches the bed (sinkage datum `z = 0`).
struct SinkageTracker {
    /// Plate plane z-coordinate at the moment of first contact (sinkage datum).
    z_contact: Option<f64>,
    /// Open results file handle (created at first contact, with a header).
    file: Option<fs::File>,
}

impl SinkageTracker {
    fn new() -> Self {
        Self { z_contact: None, file: None }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin) // bed settles under gravity before the plate presses
        .add_plugins(WallPlugin);

    app.add_resource(SinkageTracker::new());
    app.add_update_system(record_sinkage, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

/// Reads the moving plate (the single plane wall whose `velocity[2] < 0`) and
/// streams its sinkage and the vertical reaction force on it. The plate's
/// footprint is the central, downward-facing plane; the static container walls
/// (floor, sides) have zero or upward normals and are skipped.
fn record_sinkage(
    walls: Res<Walls>,
    run_state: Res<RunState>,
    input: Res<Input>,
    atoms: Res<Atom>,
    mut tracker: ResMut<SinkageTracker>,
) {
    // Identify the plate: the plane wall driven downward (velocity_z < 0).
    let plate = match walls.planes.iter().find(|w| w.velocity[2] < 0.0) {
        Some(p) => p,
        None => return,
    };

    let step = run_state.total_cycle;
    let dt = atoms.dt;
    let t = step as f64 * dt;

    // `force_accumulator` sums the scalar contact force f_net (repulsive, along
    // the plate normal -z) from all particles in the footprint. Its magnitude is
    // the vertical load the bed exerts back on the plate.
    let f_reaction = plate.force_accumulator.abs();

    // Establish the sinkage datum at first *bed* contact. A threshold well above
    // a single grazing grain (a grain's weight is ~4e-3 N here) but far below the
    // bearing load avoids latching the datum on a stray settling particle.
    const CONTACT_THRESHOLD_N: f64 = 0.05;
    let in_contact = f_reaction > CONTACT_THRESHOLD_N;
    if tracker.z_contact.is_none() {
        if !in_contact {
            return; // plate still descending through free space
        }
        tracker.z_contact = Some(plate.point_z);

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_plate_sinkage".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/plate_sinkage_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(f, "time,sinkage,force").unwrap();
        tracker.file = Some(f);
        println!("=== Plate Sinkage: contact at z = {:.6e} m, streaming results ===", plate.point_z);
    }

    let z_contact = tracker.z_contact.unwrap();
    let sinkage = (z_contact - plate.point_z).max(0.0); // downward travel since contact
    if let Some(f) = tracker.file.as_mut() {
        writeln!(f, "{:.8e},{:.8e},{:.8e}", t, sinkage, f_reaction).unwrap();
    }
}
