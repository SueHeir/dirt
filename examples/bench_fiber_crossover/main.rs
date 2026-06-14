//! bench_fiber_crossover — validates inter-fiber Coulomb friction at a single
//! crossover contact against the exact F_slide = μ·N limit.
//!
//! Two bonded-sphere fibers cross perpendicularly and touch at one crossover
//! contact (the seam between intra-fiber *bonds* and inter-fiber *contact*).
//! The lower fiber lies along y at z = 0 and is frozen. The upper fiber lies
//! along x at z ≈ 2r, is pressed down by a known normal load N (constant
//! `[[addforce]]` fz on every upper sphere), and is dragged along +x at a slow
//! constant velocity (`[[move_linear]]`). The single crossover contact must
//! resist the drag up to the Coulomb limit μ·N, then slide.
//!
//! ## Isolating the crossover contact force
//!
//! The recorder runs in the **Force** phase, *after* `hertz_mindlin_contact`
//! and `dem_bond_force` but *before* any PostForce fix (addforce / move_linear /
//! freeze). At that instant `atoms.force` holds only contact + bond
//! contributions. Summing `atoms.force` over **all** upper-fiber spheres makes
//! every intra-fiber bond force cancel (Newton's third law within the fiber),
//! leaving exactly the inter-fiber crossover contact force on the upper fiber:
//!   • Σ F_z  = upward normal reaction  → balances the applied load N
//!   • Σ F_x  = tangential reaction     → rises, then plateaus at −μ·N on slip
//! No per-contact force API is needed, and no core crate is touched.
//!
//! ```bash
//! cargo run --release --example bench_fiber_crossover --no-default-features -- examples/bench_fiber_crossover/config.toml
//! ```

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs::{self, File};
use std::io::{BufWriter, Write as IoWrite};

/// Upper-fiber spheres start at z ≈ 2 mm; the lower fiber sits at z = 0. Any
/// sphere above this threshold belongs to the (dragged) upper fiber.
const UPPER_Z_THRESHOLD: f64 = 1.0e-3;

struct Recorder {
    writer: Option<BufWriter<File>>,
    record_every: usize,
    initialized: bool,
    /// Drag speed, cached at first record so displacement = v · t.
    drag_speed: f64,
}

impl Recorder {
    fn new() -> Self {
        Recorder {
            writer: None,
            record_every: 50,
            initialized: false,
            drag_speed: 0.0,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin)
        .add_plugins(DemBondPlugin);

    app.add_resource(Recorder::new());
    // Record in the Force phase, after both force producers have run, so the
    // upper-fiber force sum isolates the crossover contact (bonds cancel and
    // no PostForce fix has fired yet).
    app.add_update_system(
        record_crossover
            .after("hertz_mindlin_contact")
            .after("dem_bond_force"),
        ParticleSimScheduleSet::Force,
    );
    app.start();
}

/// System: sum the contact force the crossover exerts on the upper fiber and
/// the upper fiber's mean drag velocity, written every `record_every` steps.
fn record_crossover(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut recorder: ResMut<Recorder>,
) {
    let nlocal = atoms.nlocal as usize;
    if nlocal == 0 {
        return;
    }
    let dem = registry.expect::<DemAtom>("record_crossover");

    if !recorder.initialized {
        // Cache the prescribed drag speed from the fastest upper-fiber sphere
        // (move_linear sets it before the first force evaluation).
        let mut vmax = 0.0_f64;
        for i in 0..nlocal {
            if atoms.pos[i][2] > UPPER_Z_THRESHOLD {
                vmax = vmax.max(atoms.vel[i][0]);
            }
        }

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_fiber_crossover".to_string());
        fs::create_dir_all(format!("{}/data", out_dir)).ok();
        let path = format!("{}/data/fiber_crossover_results.csv", out_dir);
        let mut w = BufWriter::new(
            File::create(&path).unwrap_or_else(|e| panic!("Cannot create {}: {}", path, e)),
        );
        // f_normal: Σ F_z over upper fiber (upward contact reaction = applied N)
        // f_tangential: Σ F_x over upper fiber (drag-resisting contact reaction)
        // overlap: crossover centre overlap; v_drag: mean upper-fiber +x velocity
        writeln!(
            w,
            "step,t,displacement,f_normal,f_tangential,overlap,v_drag"
        )
        .unwrap();
        recorder.writer = Some(w);
        recorder.drag_speed = vmax;
        recorder.initialized = true;

        println!("=== Fiber Crossover Friction Test ===");
        println!("  drag speed v = {:.6e} m/s along +x", vmax);
        println!("  results -> {}", path);
    }

    let step = run_state.total_cycle;
    if step % recorder.record_every != 0 {
        return;
    }

    // Sum the contact force on the upper fiber (intra-fiber bonds cancel).
    let mut fx_sum = 0.0;
    let mut fz_sum = 0.0;
    let mut vx_sum = 0.0;
    let mut n_upper = 0usize;
    // Locate the two crossover spheres (closest upper/lower pair) for overlap.
    let mut mid_upper: Option<usize> = None;
    let mut mid_lower: Option<usize> = None;
    for i in 0..nlocal {
        if atoms.pos[i][2] > UPPER_Z_THRESHOLD {
            fx_sum += atoms.force[i][0];
            fz_sum += atoms.force[i][2];
            vx_sum += atoms.vel[i][0];
            n_upper += 1;
            if atoms.pos[i][0].abs() < 1.0e-3 && atoms.pos[i][1].abs() < 1.0e-3 {
                mid_upper = Some(i);
            }
        } else if atoms.pos[i][0].abs() < 1.0e-3 && atoms.pos[i][1].abs() < 1.0e-3 {
            mid_lower = Some(i);
        }
    }
    if n_upper == 0 {
        return;
    }

    let overlap = match (mid_upper, mid_lower) {
        (Some(u), Some(l)) => {
            let dx = atoms.pos[u][0] - atoms.pos[l][0];
            let dy = atoms.pos[u][1] - atoms.pos[l][1];
            let dz = atoms.pos[u][2] - atoms.pos[l][2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            (dem.radius[u] + dem.radius[l]) - dist
        }
        _ => 0.0,
    };

    let dt = atoms.dt;
    let t = step as f64 * dt;
    let displacement = recorder.drag_speed * t;
    let v_drag = vx_sum / n_upper as f64;

    if let Some(ref mut w) = recorder.writer {
        writeln!(
            w,
            "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
            step, t, displacement, fz_sum, fx_sum, overlap, v_drag
        )
        .ok();
    }
}
