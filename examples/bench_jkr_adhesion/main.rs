//! JKR adhesion pull-off benchmark — validates the adhesive pull-off force of
//! two spheres in contact against the JKR analytical result F = (3/2)·π·w·R*.
//!
//! Two identical spheres are used. The left sphere is frozen; the right sphere
//! is launched slowly inward, makes adhesive contact, and then separates. As it
//! separates it passes through the gap (adhesion-only) regime, where the contact
//! resists with a tensile (attractive) normal force — the pull-off force — then
//! snaps to zero once the spheres clear the adhesion range. This recorder logs
//! the per-step contact normal force versus separation and captures the peak
//! tensile force in the gap regime as the measured F_pulloff.
//!
//! DIRT's contact model (see `dirt_granular::contact`) implements adhesion as a
//! *constant* attractive force: JKR uses F_adh = (3/2)·π·γ·R* and DMT uses
//! F_dmt = 2·π·γ·R*, where γ is the material `surface_energy` (the work of
//! adhesion w) and R* the effective radius. The default `adhesion_model` is
//! "jkr". In the gap regime (no geometric overlap) the net normal force is
//! exactly this constant adhesion plateau, with no Hertz spring or velocity
//! damping mixed in — which is why the measurement is restricted to that regime.
//!
//! The contact force is read from `atoms.force` of the *moving* sphere at
//! `PostFinalIntegration`. The moving sphere carries no fix, so the only force
//! accumulated on it is the contact force with its single partner (the frozen
//! sphere's own force is zeroed by the freeze fix in PostForce, but we never
//! read it). The spheres start ≥ 1.1 diameters apart to satisfy the insertion
//! overlap check, and the approach velocity is small so the ~nm-scale adhesion
//! window is sampled by many steps.
//!
//! ```bash
//! cargo run --release --example bench_jkr_adhesion --no-default-features -- examples/bench_jkr_adhesion/config.toml
//! ```

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

/// Tracks the separation/force history and the peak tensile (pull-off) force.
struct PulloffTracker {
    /// Tag of the moving (pulled) sphere — the one whose force we read.
    moving_tag: Option<u32>,
    /// Tag of the frozen sphere (its force is zeroed by the freeze fix).
    frozen_tag: Option<u32>,
    /// Most negative (tensile) contact normal force seen so far [N].
    /// Negative = attractive; the pull-off force is its magnitude.
    min_fn: f64,
    /// Separation (gap between surfaces) at which the peak tensile force occurred [m].
    sep_at_pulloff: f64,
    /// True once the spheres have separated past the adhesion range.
    finished: bool,
    /// Whether the contact has ever been engaged (force became non-trivial).
    engaged: bool,
    /// Per-step trace rows: (step, separation, fn).
    trace: Vec<(usize, f64, f64)>,
    output_dir: String,
}

impl PulloffTracker {
    fn new() -> Self {
        Self {
            moving_tag: None,
            frozen_tag: None,
            min_fn: 0.0,
            sep_at_pulloff: 0.0,
            finished: false,
            engaged: false,
            trace: Vec::new(),
            output_dir: String::new(),
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin); // [[freeze]] + [[move_linear]]

    app.add_resource(PulloffTracker::new());
    app.add_update_system(track_pulloff, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn index_of_tag(atoms: &Atom, tag: u32) -> Option<usize> {
    (0..atoms.nlocal as usize).find(|&i| atoms.tag[i] == tag)
}

fn track_pulloff(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<PulloffTracker>,
) {
    if tracker.finished || atoms.nlocal < 2 {
        return;
    }
    let dem = registry.expect::<DemAtom>("track_pulloff");
    let step = run_state.total_cycle;

    // Identify the moving sphere once: the one with non-zero velocity.
    // (The frozen sphere is held at rest by the freeze fix.)
    if tracker.moving_tag.is_none() {
        let mut moving = None;
        let mut frozen = None;
        for i in 0..atoms.nlocal as usize {
            let v = [atoms.vel[i][0] as f64, atoms.vel[i][1] as f64, atoms.vel[i][2] as f64];
            let speed2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            if speed2 > 0.0 {
                moving = Some(atoms.tag[i]);
            } else {
                frozen = Some(atoms.tag[i]);
            }
        }
        tracker.moving_tag = moving;
        tracker.frozen_tag = frozen;
        if tracker.moving_tag.is_none() || tracker.frozen_tag.is_none() {
            return;
        }
    }

    let m = match index_of_tag(&atoms, tracker.moving_tag.unwrap()) {
        Some(i) => i,
        None => return,
    };
    let f = match index_of_tag(&atoms, tracker.frozen_tag.unwrap()) {
        Some(i) => i,
        None => return,
    };

    // Line of centers (frozen -> moving) and surface separation.
    let d = [
        atoms.pos[m][0] as f64 - atoms.pos[f][0] as f64,
        atoms.pos[m][1] as f64 - atoms.pos[f][1] as f64,
        atoms.pos[m][2] as f64 - atoms.pos[f][2] as f64,
    ];
    let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if dist == 0.0 {
        return;
    }
    let n = [d[0] / dist, d[1] / dist, d[2] / dist];
    let separation = dist - (dem.radius[m] + dem.radius[f]); // <0 overlap, >0 gap

    // Contact normal force on the moving sphere. With a single partner this is
    // the contact force; project onto the line of centers. The sign convention:
    // f·n > 0 pushes the moving sphere outward (repulsion), f·n < 0 pulls it
    // back toward the frozen sphere (adhesion/tension).
    let fvec = [atoms.force[m][0] as f64, atoms.force[m][1] as f64, atoms.force[m][2] as f64];
    let f_n = fvec[0] * n[0] + fvec[1] * n[1] + fvec[2] * n[2];

    // Consider the contact "engaged" once a meaningful force is present.
    if f_n.abs() > 0.0 {
        tracker.engaged = true;
    }

    // Record the per-step trace.
    tracker.trace.push((step, separation, f_n));

    // The pull-off force is the peak tension the contact sustains as it
    // separates. DIRT's adhesion is a constant attractive force that, in the
    // gap regime (no geometric overlap, separation >= 0), is the *pure*
    // adhesion plateau −F_adh with no Hertz repulsion or velocity-dependent
    // damping contaminating it — this is exactly the JKR/DMT pull-off force.
    // We therefore take the most tensile force seen while separation >= 0.
    //
    // (In the overlap regime the net force also carries the Hertz spring and a
    // velocity-dependent damping term; those are transient and not the
    // adhesive pull-off, so they are excluded from the measurement.)
    if separation >= 0.0 && f_n < tracker.min_fn {
        tracker.min_fn = f_n;
        tracker.sep_at_pulloff = separation;
    }

    // Snap-off: once the contact has engaged and the spheres are clearly
    // separating (positive gap) with the force back to exactly zero, the contact
    // has broken — the adhesion range has been cleared. Write the result.
    //
    // Fallback: also finish on the final step so a result is always written even
    // if the spheres stick (strong adhesion can capture the sphere in the well
    // so it never escapes). The peak tensile force in the gap regime is already
    // captured on the approach, so the measurement is valid either way.
    // `cycle_remaining[i]` holds the stage's total step count; `cycle_count[i]`
    // counts steps done, so the last step is when count has reached total - 1.
    let last_step = match (run_state.cycle_count.first(), run_state.cycle_remaining.first()) {
        (Some(&done), Some(&total)) => total > 0 && done + 1 >= total,
        _ => false,
    };
    let snapped = tracker.engaged && separation > 0.0 && f_n == 0.0;
    if snapped || last_step {
        finish(&mut tracker, &atoms, &dem, m, f, &input);
    }
}

fn finish(
    tracker: &mut PulloffTracker,
    atoms: &Atom,
    dem: &DemAtom,
    m: usize,
    f: usize,
    input: &Input,
) {
    tracker.finished = true;

    let r_i = dem.radius[m];
    let r_j = dem.radius[f];
    let r_eff = (r_i * r_j) / (r_i + r_j);

    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_jkr_adhesion".to_string());
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    tracker.output_dir = out_dir.clone();

    // Force–separation trace (one row per step).
    let trace_file = format!("{}/jkr_trace.csv", data_dir);
    if let Ok(mut tf) = fs::File::create(&trace_file) {
        writeln!(tf, "step,separation,f_normal").ok();
        for (s, sep, fnv) in &tracker.trace {
            writeln!(tf, "{},{:.10e},{:.10e}", s, sep, fnv).ok();
        }
    }

    // Summary row: measured pull-off force + geometry needed for theory.
    let f_pulloff = tracker.min_fn.abs();
    let results_file = format!("{}/jkr_results.csv", data_dir);
    let mut fh = fs::File::create(&results_file)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
    writeln!(fh, "f_pulloff,sep_at_pulloff,r_eff,radius,density,dt").unwrap();
    writeln!(
        fh,
        "{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
        f_pulloff,
        tracker.sep_at_pulloff,
        r_eff,
        r_i,
        dem.density[m],
        atoms.dt,
    )
    .unwrap();

    println!("=== JKR Adhesion Pull-off Results ===");
    println!("  R*              : {:.6e} m", r_eff);
    println!("  F_pulloff (meas): {:.6e} N", f_pulloff);
    println!("  separation@peak : {:.6e} m", tracker.sep_at_pulloff);
    println!("  results -> {}", results_file);
}
