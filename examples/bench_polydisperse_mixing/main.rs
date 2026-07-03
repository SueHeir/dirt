//! Polydisperse / multi-material mixing benchmark.
//!
//! Two free spheres of (in general) **unequal radius** and **different material**
//! collide once. The example measures the observables of that single binary
//! contact and writes them to CSV; `sweep.py` compares them to Hertz theory
//! evaluated with the *mixed* pair quantities:
//!
//!   * effective (reduced) radius   R* = r1·r2 / (r1 + r2)
//!   * effective Young's modulus    E* = 1 / ((1−ν1²)/E1 + (1−ν2²)/E2)   (= `e_eff_ij`)
//!   * effective mass               m* = m1·m2 / (m1 + m2)
//!   * pair sliding friction        μ_ij = √(μ1·μ2)                        (= `friction_ij`)
//!
//! Because Hertz peak overlap and contact duration depend on R* and E* only
//! through those mixing rules, matching them for **unequal-radius, cross-material**
//! pairs is a direct check that the code's `r_eff = r1 r2/(r1+r2)` and its
//! per-pair `e_eff_ij` mixing are correct. A second family of cases gives one
//! sphere a large tangential velocity so the contact is in **gross sliding**; the
//! ratio of tangential to normal impulse then equals μ_ij, exercising the
//! `friction_ij` mixing rule.
//!
//! ```bash
//! cargo run --release --example bench_polydisperse_mixing --no-default-features \
//!     --features precision-double -- <config.toml>
//! ```
//!
//! All observables are computed from the two particles order-independently (the
//! contact normal is taken from the line of centres at first contact), so the
//! measurement does not depend on which atom index ends up first.

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

/// Tracks the single binary collision between the two particles.
struct CollisionTracker {
    /// True once the pair has been in contact.
    was_in_contact: bool,
    /// True once the pair has separated after contact (measurement complete).
    finished: bool,
    /// Contact normal (unit, low-tag → high-tag) captured at first contact.
    n_contact: [f64; 3],
    /// Velocity of the low-tag particle just before contact.
    v1_pre: [f64; 3],
    /// Velocity of the high-tag particle just before contact.
    v2_pre: [f64; 3],
    /// Masses and radii, low-tag (particle 1) then high-tag (particle 2).
    m1: f64,
    m2: f64,
    r1: f64,
    r2: f64,
    /// Contact start/end steps for the duration.
    step_start: usize,
    step_end: usize,
    /// Maximum overlap during contact.
    max_overlap: f64,
    /// Previous-step velocities (updated while not yet in contact).
    prev_v1: [f64; 3],
    prev_v2: [f64; 3],
}

impl CollisionTracker {
    fn new() -> Self {
        Self {
            was_in_contact: false,
            finished: false,
            n_contact: [0.0; 3],
            v1_pre: [0.0; 3],
            v2_pre: [0.0; 3],
            m1: 0.0,
            m2: 0.0,
            r1: 0.0,
            r2: 0.0,
            step_start: 0,
            step_end: 0,
            max_overlap: 0.0,
            prev_v1: [0.0; 3],
            prev_v2: [0.0; 3],
        }
    }
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin); // enables [[freeze]] for the frozen-target friction cases

    app.add_resource(CollisionTracker::new());
    app.add_update_system(track_collision, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Monitors the two-particle contact and, on separation, writes the measured
/// collision observables to `<output_dir>/data/collision_results.csv`.
fn track_collision(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<CollisionTracker>,
) {
    if tracker.finished || atoms.nlocal != 2 {
        return;
    }
    let dem = registry.expect::<DemAtom>("track_collision");
    let step = run_state.total_cycle;

    // Order the two atoms by tag so "particle 1" is always the first CSV row.
    let (i1, i2) = if atoms.tag[0] <= atoms.tag[1] { (0, 1) } else { (1, 0) };

    let p1 = [atoms.pos[i1][0] as f64, atoms.pos[i1][1] as f64, atoms.pos[i1][2] as f64];
    let p2 = [atoms.pos[i2][0] as f64, atoms.pos[i2][1] as f64, atoms.pos[i2][2] as f64];
    let v1 = [atoms.vel[i1][0] as f64, atoms.vel[i1][1] as f64, atoms.vel[i1][2] as f64];
    let v2 = [atoms.vel[i2][0] as f64, atoms.vel[i2][1] as f64, atoms.vel[i2][2] as f64];

    let d = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let dist = dot(&d, &d).sqrt();
    let sum_r = dem.radius[i1] + dem.radius[i2];
    let overlap = sum_r - dist;
    let in_contact = overlap > 0.0;

    if !tracker.was_in_contact && !in_contact {
        // Pre-contact: remember this step's velocities for the impact snapshot.
        tracker.prev_v1 = v1;
        tracker.prev_v2 = v2;
    } else if !tracker.was_in_contact && in_contact {
        // First contact.
        tracker.was_in_contact = true;
        let n = [d[0] / dist, d[1] / dist, d[2] / dist];
        tracker.n_contact = n;
        tracker.v1_pre = tracker.prev_v1;
        tracker.v2_pre = tracker.prev_v2;
        tracker.m1 = atoms.mass[i1] as f64;
        tracker.m2 = atoms.mass[i2] as f64;
        tracker.r1 = dem.radius[i1];
        tracker.r2 = dem.radius[i2];
        tracker.step_start = step;
        tracker.max_overlap = overlap;
    } else if tracker.was_in_contact && in_contact {
        if overlap > tracker.max_overlap {
            tracker.max_overlap = overlap;
        }
    } else if tracker.was_in_contact && !in_contact {
        // Separation — take the post-contact velocities and write results.
        tracker.finished = true;
        tracker.step_end = step;
        let dt = atoms.dt;
        let n = tracker.n_contact;

        // Normal relative velocity (particle-2 minus particle-1) along n.
        let dv_pre = [
            tracker.v2_pre[0] - tracker.v1_pre[0],
            tracker.v2_pre[1] - tracker.v1_pre[1],
            tracker.v2_pre[2] - tracker.v1_pre[2],
        ];
        let dv_post = [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]];
        let vn_impact = dot(&dv_pre, &n); // negative: approaching
        let vn_rebound = dot(&dv_post, &n);
        let cor = if vn_impact != 0.0 { (vn_rebound / vn_impact).abs() } else { 0.0 };

        // Impulse delivered to particle 1: J = m1 (v1_post − v1_pre).
        let j1 = [
            tracker.m1 * (v1[0] - tracker.v1_pre[0]),
            tracker.m1 * (v1[1] - tracker.v1_pre[1]),
            tracker.m1 * (v1[2] - tracker.v1_pre[2]),
        ];
        let jn = dot(&j1, &n); // negative (particle 1 pushed along −n)
        let jt_vec = [j1[0] - jn * n[0], j1[1] - jn * n[1], j1[2] - jn * n[2]];
        let jt = dot(&jt_vec, &jt_vec).sqrt();
        let jn_mag = jn.abs();

        let contact_steps = tracker.step_end - tracker.step_start;
        let contact_time = contact_steps as f64 * dt;

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_polydisperse_mixing/data".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/data/collision_results.csv", out_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "v_n_impact,v_n_rebound,cor,contact_time,max_overlap,jt,jn,m1,m2,r1,r2"
        )
        .unwrap();
        writeln!(
            f,
            "{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}",
            vn_impact.abs(),
            vn_rebound.abs(),
            cor,
            contact_time,
            tracker.max_overlap,
            jt,
            jn_mag,
            tracker.m1,
            tracker.m2,
            tracker.r1,
            tracker.r2,
        )
        .unwrap();

        println!("=== Polydisperse mixing collision ===");
        println!("  r1, r2:            {:.4e}, {:.4e} m", tracker.r1, tracker.r2);
        println!("  m1, m2:            {:.4e}, {:.4e} kg", tracker.m1, tracker.m2);
        println!("  |v_n| impact:      {:.6e} m/s", vn_impact.abs());
        println!("  |v_n| rebound:     {:.6e} m/s", vn_rebound.abs());
        println!("  COR (normal):      {:.6}", cor);
        println!("  Contact duration:  {:.6e} s ({} steps)", contact_time, contact_steps);
        println!("  Peak overlap:      {:.6e} m", tracker.max_overlap);
        println!(
            "  |Jt|, |Jn|:        {:.6e}, {:.6e}  (Jt/Jn = {:.6})",
            jt,
            jn_mag,
            jt / jn_mag
        );
        println!("  Results saved to:  {}", results_file);
    }

    if !tracker.finished && !in_contact {
        tracker.prev_v1 = v1;
        tracker.prev_v2 = v2;
    }
}
