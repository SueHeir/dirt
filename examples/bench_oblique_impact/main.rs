//! Oblique-impact benchmark — validates the tangential contact model
//! (Mindlin spring + Coulomb friction cap) against the exact rigid-body
//! impulse result for a sphere striking a frictional surface in the
//! gross-sliding regime.
//!
//! A projectile sphere impacts a *frozen* target sphere (an immovable but
//! fully frictional contact partner — `dirt_wall` has no sliding friction,
//! so the target must be a real sphere). The projectile is launched from an
//! offset so it arrives nearly above the target; velocities are decomposed in
//! the actual impact frame (line-of-centers normal n̂, in-plane tangent t̂),
//! so the validation is robust to small geometric tilt and uses the *measured*
//! impact conditions.
//!
//! Gross-sliding analytical reference (mass m, radius R, I = (2/5) m R²,
//! normal restitution e, friction μ, no initial spin):
//!   v_n' =  e · v_n
//!   v_t' =  v_t − μ (1+e) v_n
//!   |ω'| =  5 μ (1+e) v_n / (2 R)
//! valid while v_t > (7/2) μ (1+e) v_n  (contact point slides throughout).
//!
//! ```bash
//! cargo run --release --example bench_oblique_impact --no-default-features -- examples/bench_oblique_impact/config.toml
//! ```

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use std::fs;
use std::io::Write as IoWrite;

struct ObliqueTracker {
    projectile_tag: Option<u32>,
    was_in_contact: bool,
    finished: bool,
    // Impact frame, captured at first contact.
    nhat: [f64; 3],
    that: [f64; 3],
    vn_impact: f64,
    vt_impact: f64,
    // Rebound, captured at separation.
    vn_rebound: f64,
    vt_rebound: f64,
    omega_y_rebound: f64,
    step_contact_start: usize,
    step_contact_end: usize,
    max_overlap: f64,
    // Previous-step projectile velocity (to capture pre-contact value).
    prev_vel: [f64; 3],
}

impl ObliqueTracker {
    fn new() -> Self {
        Self {
            projectile_tag: None,
            was_in_contact: false,
            finished: false,
            nhat: [0.0; 3],
            that: [0.0; 3],
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

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin); // for [[freeze]] on the target sphere

    app.add_resource(ObliqueTracker::new());
    app.add_update_system(track_oblique, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn index_of_tag(atoms: &Atom, tag: u32) -> Option<usize> {
    (0..atoms.nlocal as usize).find(|&i| atoms.tag[i] == tag)
}

fn track_oblique(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<ObliqueTracker>,
) {
    if tracker.finished || atoms.nlocal < 2 {
        return;
    }
    let dem = registry.expect::<DemAtom>("track_oblique");
    let step = run_state.total_cycle;

    // Identify the projectile once: the moving sphere (target is frozen).
    if tracker.projectile_tag.is_none() {
        let mut best = (0usize, 0.0f64);
        for i in 0..atoms.nlocal as usize {
            let vi = [atoms.vel[i][0] as f64, atoms.vel[i][1] as f64, atoms.vel[i][2] as f64];
            let s = dot(vi, vi);
            if s > best.1 {
                best = (i, s);
            }
        }
        tracker.projectile_tag = Some(atoms.tag[best.0]);
    }
    let p = match index_of_tag(&atoms, tracker.projectile_tag.unwrap()) {
        Some(i) => i,
        None => return,
    };
    let t = (0..atoms.nlocal as usize).find(|&i| i != p).unwrap();

    let vel = [atoms.vel[p][0] as f64, atoms.vel[p][1] as f64, atoms.vel[p][2] as f64];

    // Line-of-centers from target → projectile, and center overlap.
    let d = [
        atoms.pos[p][0] as f64 - atoms.pos[t][0] as f64,
        atoms.pos[p][1] as f64 - atoms.pos[t][1] as f64,
        atoms.pos[p][2] as f64 - atoms.pos[t][2] as f64,
    ];
    let dist = dot(d, d).sqrt();
    let overlap = (dem.radius[p] + dem.radius[t]) - dist;
    let in_contact = overlap > 0.0 && dist > 0.0;

    // Optional per-step contact trace (set DIRT_TRACE=1): logs the projectile's
    // contact force decomposed into normal/tangential, for comparison with LAMMPS.
    // atoms.force[p] holds the contact force on the projectile (single partner).
    // Written to `<output_dir>/contact_trace.csv` (columns:
    // step,overlap,fn,ft_mag,ft_signed,vt_x,omega_y).
    if in_contact && std::env::var("DIRT_TRACE").is_ok() {
        let n = [d[0] / dist, d[1] / dist, d[2] / dist];
        let f = [atoms.force[p][0] as f64, atoms.force[p][1] as f64, atoms.force[p][2] as f64];
        let fn_ = dot(f, n);
        let ftv = [f[0] - fn_ * n[0], f[1] - fn_ * n[1], f[2] - fn_ * n[2]];
        let ft_mag = dot(ftv, ftv).sqrt();
        let ft_signed = ftv[0]; // tangent ≈ +x at near-vertical impact
        let vn_p = dot(vel, n);
        let vt_x = vel[0] - vn_p * n[0];
        let trace_path = format!(
            "{}/contact_trace.csv",
            input.output_dir.clone().unwrap_or_else(|| ".".to_string())
        );
        use std::io::Write as _;
        if let Ok(mut tf) = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&trace_path)
        {
            let _ = writeln!(
                tf,
                "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
                step, overlap, fn_, ft_mag, ft_signed, vt_x, dem.omega[p][1]
            );
        }
    }

    if !tracker.was_in_contact && !in_contact {
        tracker.prev_vel = vel;
    } else if !tracker.was_in_contact && in_contact {
        // First contact: build the impact frame from the pre-contact velocity.
        tracker.was_in_contact = true;
        let n = [d[0] / dist, d[1] / dist, d[2] / dist];
        let v = tracker.prev_vel;
        let vn = -dot(v, n); // closing speed (positive when approaching)
        // tangential part of pre-contact velocity
        let vt_vec = [v[0] + vn * n[0], v[1] + vn * n[1], v[2] + vn * n[2]];
        let vt_mag = dot(vt_vec, vt_vec).sqrt();
        let that = if vt_mag > 1e-12 {
            [vt_vec[0] / vt_mag, vt_vec[1] / vt_mag, vt_vec[2] / vt_mag]
        } else {
            [0.0; 3]
        };
        tracker.nhat = n;
        tracker.that = that;
        tracker.vn_impact = vn;
        tracker.vt_impact = vt_mag;
        tracker.step_contact_start = step;
        tracker.max_overlap = overlap;
    } else if tracker.was_in_contact && in_contact {
        if overlap > tracker.max_overlap {
            tracker.max_overlap = overlap;
        }
    } else if tracker.was_in_contact && !in_contact {
        // Separation: project the rebound velocity onto the stored impact frame.
        tracker.finished = true;
        tracker.vn_rebound = dot(vel, tracker.nhat);
        tracker.vt_rebound = dot(vel, tracker.that);
        tracker.omega_y_rebound = dem.omega[p][1];
        tracker.step_contact_end = step;

        let dt = atoms.dt;
        let contact_steps = tracker.step_contact_end - tracker.step_contact_start;
        let contact_time = contact_steps as f64 * dt;

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_oblique_impact".to_string());
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

        println!("=== Oblique Impact Results (impact frame) ===");
        println!("  v_n impact:  {:.6} m/s   v_t impact: {:.6} m/s", tracker.vn_impact, tracker.vt_impact);
        println!("  v_n rebound: {:.6} m/s   v_t rebound:{:.6} m/s", tracker.vn_rebound, tracker.vt_rebound);
        println!("  omega_y:     {:.6} rad/s", tracker.omega_y_rebound);
        println!("  contact time:{:.6e} s ({} steps)", contact_time, contact_steps);
        println!("  results -> {}", results_file);
    }

    if !tracker.finished && !in_contact {
        tracker.prev_vel = vel;
    }
}
