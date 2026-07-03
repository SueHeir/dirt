//! DMT / SJKR cohesion benchmark — validates the **DMT** adhesive pull-off force
//! and the **SJKR** (simplified-JKR) area-proportional cohesion against theory,
//! and exercises DIRT's adhesion-model selection. This is deliberately distinct
//! from `bench_jkr_adhesion`, which validates the JKR pull-off constant
//! `F = (3/2)·π·w·R*`; here the two physically different models are:
//!
//!   * **DMT** (`adhesion_model = "dmt"`, material `surface_energy = w`):
//!     Hertz contact plus a *constant* attractive force `F_dmt = 2·π·w·R*`.
//!     Unlike JKR, DMT has **no gap (adhesion-only) regime** — the interaction
//!     range is not extended, so the pull-off is realized *inside* geometric
//!     overlap. As the contact separates quasi-statically the net normal force
//!     is `F_n(δ) = (4/3)E*√R*·δ^{3/2} − F_dmt`, which is most tensile as the
//!     overlap `δ → 0⁺`, tending to exactly `−2·π·w·R*`. The sweep extracts that
//!     `δ→0` intercept, so the measured pull-off is `2·π·w·R*` and the
//!     DMT/JKR ratio is `4/3`.
//!
//!   * **SJKR** (material `cohesion_energy = c`, no `surface_energy`):
//!     cohesion proportional to the circular contact area `A = π·R*·δ`, giving
//!     an attractive force `F_coh(δ) = c·π·R*·δ` that is **linear in overlap**
//!     and vanishes at separation (no constant pull-off plateau). The sweep
//!     isolates it by differencing an SJKR run against a pure-Hertz baseline run
//!     at matched overlap (both run at restitution = 1 so the normal force is
//!     conservative — no velocity damping — and the shared Hertz term cancels).
//!
//! This recorder is model-agnostic: it logs the *free* sphere's contact normal
//! force versus surface separation every step. The left sphere is frozen; the
//! right sphere is launched slowly inward, makes contact, compresses, and
//! separates. The free sphere carries no fix, so at `PostFinalIntegration` its
//! `atoms.force` is exactly the single contact force (the frozen partner's force
//! is zeroed by the freeze fix, but we never read it). The full force–separation
//! trace drives both validations in `sweep.py`; the peak gap/overlap tension is
//! also captured for a quick console readout.
//!
//! ```bash
//! cargo run --release --example bench_dmt_sjkr_cohesion --no-default-features \
//!   --features precision-double -- examples/bench_dmt_sjkr_cohesion/config.toml
//! ```

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

/// Tracks the per-step force–separation history of the free sphere and the peak
/// tensile (attractive) contact force.
struct CohesionTracker {
    /// Tag of the moving (free) sphere — the one whose force we read.
    moving_tag: Option<u32>,
    /// Tag of the frozen sphere (its force is zeroed by the freeze fix).
    frozen_tag: Option<u32>,
    /// Most negative (tensile) contact normal force seen so far [N].
    min_fn: f64,
    /// Surface separation at which the peak tensile force occurred [m]
    /// (<0 overlap, >0 gap).
    sep_at_min: f64,
    /// Whether the contact has ever engaged (non-zero force seen).
    engaged: bool,
    /// True once the result has been written.
    finished: bool,
    /// Per-step trace rows: (step, separation, f_normal, v_normal).
    trace: Vec<(usize, f64, f64, f64)>,
    output_dir: String,
}

impl CohesionTracker {
    fn new() -> Self {
        Self {
            moving_tag: None,
            frozen_tag: None,
            min_fn: 0.0,
            sep_at_min: 0.0,
            engaged: false,
            finished: false,
            trace: Vec::new(),
            output_dir: String::new(),
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin); // [[freeze]]

    app.add_resource(CohesionTracker::new());
    app.add_update_system(track_cohesion, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn index_of_tag(atoms: &Atom, tag: u32) -> Option<usize> {
    (0..atoms.nlocal as usize).find(|&i| atoms.tag[i] == tag)
}

fn track_cohesion(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<CohesionTracker>,
) {
    if tracker.finished || atoms.nlocal < 2 {
        return;
    }
    let dem = registry.expect::<DemAtom>("track_cohesion");
    let step = run_state.total_cycle;

    // Identify the moving sphere once: the one with non-zero velocity.
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

    // Contact normal force on the free sphere, projected onto the line of
    // centers. f·n > 0 pushes the free sphere outward (repulsion), f·n < 0
    // pulls it toward the frozen sphere (adhesion/cohesion/tension).
    let fvec = [atoms.force[m][0] as f64, atoms.force[m][1] as f64, atoms.force[m][2] as f64];
    let f_n = fvec[0] * n[0] + fvec[1] * n[1] + fvec[2] * n[2];

    // Relative normal velocity (free minus frozen; frozen is at rest).
    let vvec = [atoms.vel[m][0] as f64, atoms.vel[m][1] as f64, atoms.vel[m][2] as f64];
    let v_n = vvec[0] * n[0] + vvec[1] * n[1] + vvec[2] * n[2];

    if f_n.abs() > 0.0 {
        tracker.engaged = true;
        // Only record steps where the contact is active, to keep traces small.
        tracker.trace.push((step, separation, f_n, v_n));
    }

    // Track the peak tension (most negative f_n) for a quick console readout.
    if f_n < tracker.min_fn {
        tracker.min_fn = f_n;
        tracker.sep_at_min = separation;
    }

    // Finish once the contact has engaged and cleanly separated (positive gap,
    // force back to zero), or on the last step as a fallback.
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
    tracker: &mut CohesionTracker,
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
        .unwrap_or_else(|| "examples/bench_dmt_sjkr_cohesion".to_string());
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    tracker.output_dir = out_dir.clone();

    // Force–separation trace (one row per active-contact step).
    let trace_file = format!("{}/cohesion_trace.csv", data_dir);
    if let Ok(mut tf) = fs::File::create(&trace_file) {
        writeln!(tf, "step,separation,f_normal,v_normal").ok();
        for (s, sep, fnv, vn) in &tracker.trace {
            writeln!(tf, "{},{:.10e},{:.10e},{:.10e}", s, sep, fnv, vn).ok();
        }
    }

    // Summary row: peak tension + geometry needed for theory.
    let f_peak_tension = tracker.min_fn.abs();
    let results_file = format!("{}/cohesion_results.csv", data_dir);
    let mut fh = fs::File::create(&results_file)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
    writeln!(fh, "f_peak_tension,sep_at_peak,r_eff,radius,density,dt").unwrap();
    writeln!(
        fh,
        "{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
        f_peak_tension,
        tracker.sep_at_min,
        r_eff,
        r_i,
        dem.density[m],
        atoms.dt,
    )
    .unwrap();

    println!("=== DMT / SJKR Cohesion Results ===");
    println!("  R*                 : {:.6e} m", r_eff);
    println!("  peak tension (meas): {:.6e} N", f_peak_tension);
    println!("  separation@peak    : {:.6e} m", tracker.sep_at_min);
    println!("  trace   -> {}", trace_file);
    println!("  results -> {}", results_file);
}
