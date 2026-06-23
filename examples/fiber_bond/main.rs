//! Comprehensive BPM fiber-bond validation harness.
//!
//! One binary, multiple configs — each TOML in this directory exercises a
//! different deformation mode and writes a CSV of bond / endpoint kinematics
//! to `<output.dir>/data/fiber_bond.csv`. The companion `validate.py` reads
//! the CSV and compares against analytical predictions from Guo et al. 2018
//! (*Chem. Eng. Sci.* **175**, 118–129).
//!
//! Modes covered today:
//!
//! * **Axial elastic** (`axial_elastic.toml`) — fixed/pulled fiber, recovers
//!   `E_b` from the σ(ε) slope. Guo Eq. 1.
//! * **Cantilever bending elastic** (`cantilever_bending.toml`) — pinned/
//!   transverse-loaded fiber, recovers `E_b·I` from the small-deformation
//!   Euler-Bernoulli relation `y_tip = F·L_c³ / (3·E_b·I)`. Guo Sec. 2.1, Fig. 3.
//!
//! Run any config:
//! ```bash
//! cargo run --release --example fiber_bond --no-default-features -- \
//!     examples/fiber_bond/axial_elastic.toml
//! ```
//! Validate any one (or all) run results:
//! ```bash
//! python3 examples/fiber_bond/validate.py examples/fiber_bond/axial_elastic/data/fiber_bond.csv
//! ```
//!
//! ## What gets recorded
//!
//! Per sampled step (`record_every` steps):
//! * Endpoint positions and velocities (left- and right-most atoms by initial x).
//! * Middle-bond geometry: rest length, current length, axial strain, bending
//!   angle magnitude.
//! * Middle-bond plastic state: `θ_p_bend`, `ε_p_axial`, `θ_max_bend`,
//!   `ε_max_axial` (zero whenever the corresponding channel is configured
//!   elastic).
//! * Bond count and `bonds_broken` (cumulative).
//!
//! The recorder intentionally does **not** capture force/moment directly —
//! those are reconstructed from the bond stiffness and the recorded
//! kinematics inside `validate.py`. That keeps the recorder cheap and
//! decoupled from whether bending or axial channels are running an elastic
//! or piecewise-plastic envelope.

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use dirt_core::dirt_bond::{BondConfig, BondHistoryStore, BondMetrics};
use dirt_core::soil_core::BondStore;
use dirt_core::dirt_fixes::FixesPlugin;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufWriter, Write as IoWrite};

/// State carried across steps by the recorder system. Set up once on the
/// first record cycle; reused thereafter.
struct Recorder {
    writer: Option<BufWriter<File>>,
    /// Output directory; cached so we can re-open `profile.csv` each sample.
    out_dir: String,
    /// **Tags** (not indices) of the leftmost and rightmost atoms by initial
    /// x. Tags are stable across engine atom-reordering; the recorder
    /// re-resolves the local index each sample.
    tag_left: Option<u32>,
    tag_right: Option<u32>,
    /// Tags of the two atoms anchoring the middle bond.
    tag_mid_a: u32,
    tag_mid_b: u32,
    /// Per-atom initial positions, indexed by atom tag. Captured at setup so
    /// `profile.csv` can give a (tag, x₀, y₀, z₀, x, y, z) row regardless of
    /// any later atom reordering inside the engine.
    initial_pos: Vec<[f64; 3]>,
    /// Initial geometry, cached at setup.
    length0: f64,
    bond_len_mid0: f64,
    r_b: f64,
    area: f64,
    iben: f64,
    k_n: f64,
    k_bend: f64,
    record_every: usize,
    initialized: bool,
}

impl Recorder {
    fn new() -> Self {
        Recorder {
            writer: None,
            out_dir: String::new(),
            tag_left: None,
            tag_right: None,
            tag_mid_a: 0,
            tag_mid_b: 0,
            initial_pos: Vec::new(),
            length0: 0.0,
            bond_len_mid0: 0.0,
            r_b: 0.0,
            area: 0.0,
            iben: 0.0,
            k_n: 0.0,
            k_bend: 0.0,
            record_every: 200,
            initialized: false,
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
    // Time-gated load schedule for the `bending_plastic_guo` scenario.
    // No-op for every other config: the system inspects the output dir name
    // and only applies the Guo three-step force schedule when matched.
    app.add_update_system(apply_three_step_load, ParticleSimScheduleSet::Force);
    app.add_update_system(record_fiber_state, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

// ── Guo 2018 three-step load schedule ───────────────────────────────────────
//
// Mimics Guo 2018 Fig. 10 (compressed to a millisecond timescale): three
// successive load cycles, each a triangular pulse with a brief flat top, then
// a settle interval at zero force. Each cycle drives the chain past the
// plastic-bending yield and the residual `θ_p` is expected to accumulate from
// one cycle to the next (Guo Fig. 12).
//
//   t ∈ [0,     T_ramp):     F ramps 0 → F_peak       (cycle 1)
//   t ∈ [T_ramp, T_hold):    F = F_peak               (cycle 1 hold)
//   t ∈ [T_hold, T_unload):  F ramps F_peak → 0       (cycle 1 unload)
//   t ∈ [T_unload, T_cycle): F = 0                    (cycle 1 settle)
//   …repeat twice more (cycles 2, 3)…
//
// Activated only for runs whose output dir contains "bending_plastic_guo" —
// every other config gets a no-op.

const F_PEAK: f64        = -0.6;     // N, downward. Anchor moment F·L = 0.024 N·m, middle moment F·L/2 = 0.012 N·m vs M_p = 0.010 N·m — first ~half of the chain yields, tip deflection stays manageable.
const T_RAMP: f64        =  5.0e-3;  // ms-scale ramp so the load varies slowly compared to T_bend ≈ 5.5 ms.
const T_HOLD_END: f64    = 15.0e-3;  // 10 ms hold — ~2 bending periods, enough to reach static deflection.
const T_UNLOAD_END: f64  = 20.0e-3;  // 5 ms unload ramp.
const T_CYCLE: f64       = 30.0e-3;  // 20 ms force pulse + 10 ms settle = 30 ms / cycle.
const N_CYCLES: usize    =  3;
const TIP_TAG: u32       = 10;       // the rightmost atom in fiber_11_spaced.csv

fn three_step_force_at(t: f64) -> f64 {
    if t < 0.0 || t >= (N_CYCLES as f64) * T_CYCLE {
        return 0.0;
    }
    let cycle_idx = (t / T_CYCLE).floor() as usize;
    let phase = t - (cycle_idx as f64) * T_CYCLE;
    if phase < T_RAMP {
        F_PEAK * (phase / T_RAMP)
    } else if phase < T_HOLD_END {
        F_PEAK
    } else if phase < T_UNLOAD_END {
        F_PEAK * (1.0 - (phase - T_HOLD_END) / (T_UNLOAD_END - T_HOLD_END))
    } else {
        0.0
    }
}

fn apply_three_step_load(
    mut atoms: ResMut<Atom>,
    input: Res<Input>,
    run_state: Res<RunState>,
) {
    let out_dir = input.output_dir.clone().unwrap_or_default();
    if !out_dir.contains("bending_plastic_guo") { return; }
    let dt = atoms.dt;
    let t = run_state.total_cycle as f64 * dt;
    let f = three_step_force_at(t);
    if f == 0.0 { return; }
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal {
        if atoms.tag[i] == TIP_TAG {
            atoms.force[i][2] += f as dirt_core::soil_core::Accum;
            break;
        }
    }
}

/// Returns the **tags** of the leftmost and rightmost atoms by initial x at
/// setup time. Tags survive engine atom-reordering (neighbour-list rebuilds),
/// whereas local indices do not — the recorder resolves index from tag each
/// sample.
fn find_endpoint_tags_by_initial_x(atoms: &Atom) -> (Option<u32>, Option<u32>) {
    let nlocal = atoms.nlocal as usize;
    if nlocal == 0 { return (None, None); }
    let mut i_left = 0usize;
    let mut i_right = 0usize;
    for i in 1..nlocal {
        if atoms.pos[i][0] < atoms.pos[i_left][0] { i_left = i; }
        if atoms.pos[i][0] > atoms.pos[i_right][0] { i_right = i; }
    }
    (Some(atoms.tag[i_left]), Some(atoms.tag[i_right]))
}

/// Find the bond closest to the centre of the fiber. Returns the lower-tag end's
/// local index and the partner tag.
fn find_middle_bond(
    atoms: &Atom,
    registry: &AtomDataRegistry,
) -> Option<(usize, u32, u32, f64)> {
    let bonds = registry.get::<BondStore>()?;
    let nlocal = atoms.nlocal as usize;
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    for i in 0..nlocal {
        x_lo = x_lo.min(atoms.pos[i][0] as f64);
        x_hi = x_hi.max(atoms.pos[i][0] as f64);
    }
    let x_centre = 0.5 * (x_lo + x_hi);
    let mut best: Option<(usize, u32, u32, f64, f64)> = None;
    for i in 0..nlocal {
        if i >= bonds.bonds.len() { break; }
        for b in &bonds.bonds[i] {
            // Process each bond once: lower tag is the owner.
            if atoms.tag[i] >= b.partner_tag { continue; }
            // Partner must also be local for the mid-bond probe — this is OK
            // for our single-rank validation configs.
            let j_opt = (0..nlocal).find(|&k| atoms.tag[k] == b.partner_tag);
            let j = match j_opt { Some(k) => k, None => continue };
            let bond_centre_x = 0.5 * (atoms.pos[i][0] as f64 + atoms.pos[j][0] as f64);
            let d = (bond_centre_x - x_centre).abs();
            let take = match &best {
                None => true,
                Some((_, _, _, _, dbest)) => d < *dbest,
            };
            if take {
                best = Some((i, atoms.tag[i], b.partner_tag, b.r0, d));
            }
        }
    }
    best.map(|(i, ta, tb, r0, _)| (i, ta, tb, r0))
}

/// Sampling system: writes one CSV row every `record_every` steps.
fn record_fiber_state(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    bond_config: Res<BondConfig>,
    bond_metrics: Res<BondMetrics>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut rec: ResMut<Recorder>,
) {
    if !rec.initialized {
        // ── one-shot setup: cache geometry & open CSV ──
        let (tag_left, tag_right) = find_endpoint_tags_by_initial_x(&atoms);
        let mid = match find_middle_bond(&atoms, &registry) {
            Some(m) => m,
            None => return, // bonds not yet populated
        };
        let (i_mid_a, tag_mid_a, tag_mid_b, r0_mid) = mid;
        // Partner local index.
        let i_mid_b = match (0..atoms.nlocal as usize).find(|&k| atoms.tag[k] == tag_mid_b) {
            Some(i) => i,
            None => return,
        };

        let dem = registry.expect::<DemAtom>("record_fiber_state");
        let r_i = dem.radius[i_mid_a];
        let r_j = dem.radius[i_mid_b];
        let r_b = bond_config.bond_radius_ratio * r_i.min(r_j);
        let area = PI * r_b * r_b;
        let iben = 0.25 * PI * r_b.powi(4);

        let k_n = match bond_config.youngs_modulus {
            Some(e) => e * area / r0_mid,
            None => bond_config.normal_stiffness,
        };
        let k_bend = match bond_config.youngs_modulus {
            Some(e) => e * iben / r0_mid,
            None => bond_config.bending_stiffness,
        };

        let nlocal_now = atoms.nlocal as usize;
        let il_idx = tag_left.and_then(|t| (0..nlocal_now).find(|&k| atoms.tag[k] == t));
        let ir_idx = tag_right.and_then(|t| (0..nlocal_now).find(|&k| atoms.tag[k] == t));
        let length0 = match (il_idx, ir_idx) {
            (Some(il), Some(ir)) => atoms.pos[ir][0] as f64 - atoms.pos[il][0] as f64,
            _ => 0.0,
        };

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/fiber_bond".to_string());
        fs::create_dir_all(format!("{}/data", out_dir)).ok();
        let path = format!("{}/data/fiber_bond.csv", out_dir);
        let mut w = BufWriter::new(
            File::create(&path).unwrap_or_else(|e| panic!("Cannot create {}: {}", path, e)),
        );
        writeln!(
            w,
            // Header — every column is dimensional (SI). Plastic-state columns
            // are zero on elastic-only runs.
            "step,t,\
             left_x,left_y,left_z,left_vx,left_vy,left_vz,\
             right_x,right_y,right_z,right_vx,right_vy,right_vz,\
             length_global,length_mid,delta_mid,strain_axial_mid,\
             delta_t_mid_mag,dth_bend_mid_mag,dth_bend_y_mid,dth_twist_mid,\
             theta_p_bend_mid_mag,theta_p_bend_y_mid,eps_p_axial_mid,\
             theta_max_bend_mid,eps_max_axial_mid,\
             r_b,area,iben,k_n,k_bend,length0,bond_len_mid0,\
             bond_count,bonds_broken"
        )
        .unwrap();

        // Capture initial positions indexed by atom tag so profile rows can
        // be matched up later regardless of any reordering inside the engine.
        let nlocal = atoms.nlocal as usize;
        let max_tag = (0..nlocal).map(|i| atoms.tag[i] as usize).max().unwrap_or(0);
        let mut initial_pos = vec![[0.0f64; 3]; max_tag + 1];
        for i in 0..nlocal {
            initial_pos[atoms.tag[i] as usize] = [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ];
        }

        rec.writer = Some(w);
        rec.out_dir = out_dir.clone();
        rec.tag_left = tag_left;
        rec.tag_right = tag_right;
        rec.tag_mid_a = tag_mid_a;
        rec.tag_mid_b = tag_mid_b;
        rec.initial_pos = initial_pos;
        rec.length0 = length0;
        rec.bond_len_mid0 = r0_mid;
        rec.r_b = r_b;
        rec.area = area;
        rec.iben = iben;
        rec.k_n = k_n;
        rec.k_bend = k_bend;
        rec.initialized = true;

        println!("=== Fiber Bond Validation Harness ===");
        println!("  output         : {}", path);
        println!("  L₀ (end-to-end): {:.6e} m", length0);
        println!("  middle bond    : tags {}↔{}, r₀ = {:.6e} m",
            tag_mid_a, tag_mid_b, r0_mid);
        println!("  bond radius    : {:.6e} m", r_b);
        println!("  area, I_ben    : {:.6e} m², {:.6e} m⁴", area, iben);
        println!("  K_n, K_bend    : {:.6e} N/m, {:.6e} N·m/rad", k_n, k_bend);

        // One-shot dump of per-bond sampled thresholds. The recorder is
        // running here for the first time, which is after `init_bond_history`
        // has populated `thresholds[0..4]` for every bond via
        // `per_bond_uniform_samples(seed, tag_a, tag_b)`. Writing the table
        // now lets the validator predict which bond will break first and at
        // what strain — without having to re-implement the SplitMix64 /
        // SmallRng / Weibull sampling chain in Python.
        let thresholds_path = format!("{}/data/bond_thresholds.csv", out_dir);
        if let Ok(f) = File::create(&thresholds_path) {
            let mut bw = BufWriter::new(f);
            writeln!(bw, "tag_a,tag_b,r0,thr0,thr1,thr2,thr3").ok();
            let hist = registry.expect::<BondHistoryStore>("threshold_dump");
            let bonds = registry.expect::<BondStore>("threshold_dump");
            let nlocal_now = atoms.nlocal as usize;
            for i in 0..nlocal_now.min(hist.history.len()) {
                let tag_a = atoms.tag[i];
                for entry in &hist.history[i] {
                    // Each bond visited once (canonical: lower-tag end owns).
                    if tag_a > entry.partner_tag { continue; }
                    let r0_b = if i < bonds.bonds.len() {
                        bonds.bonds[i]
                            .iter()
                            .find(|b| b.partner_tag == entry.partner_tag)
                            .map(|b| b.r0)
                            .unwrap_or(0.0)
                    } else { 0.0 };
                    writeln!(
                        bw,
                        "{},{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
                        tag_a, entry.partner_tag, r0_b,
                        entry.thresholds[0], entry.thresholds[1],
                        entry.thresholds[2], entry.thresholds[3],
                    ).ok();
                }
            }
            println!("  bond thresholds: {}", thresholds_path);
        }
    }

    let step = run_state.total_cycle;
    if step % rec.record_every != 0 { return; }

    let dt = atoms.dt;
    let t = step as f64 * dt;

    // Endpoint state — resolve tag → local index every sample so we are
    // robust against engine atom-reordering (neighbour-list rebuilds, MPI
    // exchange, etc.).
    let nlocal_now = atoms.nlocal as usize;
    let i_left = rec.tag_left.and_then(|t| (0..nlocal_now).find(|&k| atoms.tag[k] == t));
    let i_right = rec.tag_right.and_then(|t| (0..nlocal_now).find(|&k| atoms.tag[k] == t));
    let (lx, ly, lz, lvx, lvy, lvz) = match i_left {
        Some(i) => (atoms.pos[i][0] as f64, atoms.pos[i][1] as f64, atoms.pos[i][2] as f64,
                    atoms.vel[i][0] as f64, atoms.vel[i][1] as f64, atoms.vel[i][2] as f64),
        None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
    };
    let (rx, ry, rz, rvx, rvy, rvz) = match i_right {
        Some(i) => (atoms.pos[i][0] as f64, atoms.pos[i][1] as f64, atoms.pos[i][2] as f64,
                    atoms.vel[i][0] as f64, atoms.vel[i][1] as f64, atoms.vel[i][2] as f64),
        None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
    };
    let length_global = match (i_left, i_right) {
        (Some(_), Some(_)) => {
            let dx = rx - lx; let dy = ry - ly; let dz = rz - lz;
            (dx*dx + dy*dy + dz*dz).sqrt()
        }
        _ => f64::NAN,
    };

    // Middle-bond kinematics + plastic state.
    let nlocal = atoms.nlocal as usize;
    let i_a = (0..nlocal).find(|&k| atoms.tag[k] == rec.tag_mid_a);
    let i_b = (0..nlocal).find(|&k| atoms.tag[k] == rec.tag_mid_b);
    let (length_mid, delta_mid, strain_axial_mid,
         delta_t_mid_mag, dth_bend_mid_mag, dth_bend_y_mid, dth_twist_mid,
         theta_p_bend_mid_mag, theta_p_bend_y_mid, eps_p_axial_mid,
         theta_max_bend_mid, eps_max_axial_mid) =
        match (i_a, i_b) {
            (Some(i), Some(j)) => {
                let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
                let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
                let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
                let len = (dx*dx + dy*dy + dz*dz).sqrt();
                let delta = len - rec.bond_len_mid0;
                let strain = delta / rec.bond_len_mid0;

                let hist = registry.expect::<BondHistoryStore>("record_fiber_state");
                let entry = hist.history.get(i)
                    .and_then(|v| v.iter().find(|h| h.partner_tag == rec.tag_mid_b));
                let (dt_m, dth_bm, dth_y, dth_tw, tpb, tpb_y, epa, tmb, ema) = match entry {
                    Some(h) => {
                        let dt_m = (h.delta_t[0].powi(2) + h.delta_t[1].powi(2) + h.delta_t[2].powi(2)).sqrt();
                        let nhat = [dx / len, dy / len, dz / len];
                        let dn = h.delta_theta[0]*nhat[0]
                               + h.delta_theta[1]*nhat[1]
                               + h.delta_theta[2]*nhat[2];
                        let dth_bend = [
                            h.delta_theta[0] - dn * nhat[0],
                            h.delta_theta[1] - dn * nhat[1],
                            h.delta_theta[2] - dn * nhat[2],
                        ];
                        let dth_bm = (dth_bend[0].powi(2) + dth_bend[1].powi(2) + dth_bend[2].powi(2)).sqrt();
                        // y-component (signed) of bending vectors — for 2D
                        // bending in the xz-plane (load in −z, fiber along x),
                        // all bending lives in this single component.
                        let dth_y = dth_bend[1];
                        let tpb = (h.theta_p_bend[0].powi(2)
                                 + h.theta_p_bend[1].powi(2)
                                 + h.theta_p_bend[2].powi(2)).sqrt();
                        let tpb_y = h.theta_p_bend[1];
                        (dt_m, dth_bm, dth_y, dn, tpb, tpb_y, h.eps_p_axial, h.theta_max_bend, h.eps_max_axial)
                    }
                    None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
                };
                (len, delta, strain, dt_m, dth_bm, dth_y, dth_tw, tpb, tpb_y, epa, tmb, ema)
            }
            _ => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
        };

    let bonds_broken = bond_metrics.total_bonds_broken;
    let bond_count = bond_metrics.bond_count;

    // Snapshot cached recorder geometry so the writeln below borrows `w`
    // mutably and never re-borrows `rec`.
    let r_b = rec.r_b;
    let area = rec.area;
    let iben = rec.iben;
    let k_n = rec.k_n;
    let k_bend = rec.k_bend;
    let length0 = rec.length0;
    let bond_len_mid0 = rec.bond_len_mid0;

    if let Some(ref mut w) = rec.writer {
        writeln!(
            w,
            "{},{:.8e},\
             {:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},\
             {:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},\
             {:.8e},{:.8e},{:.8e},{:.8e},\
             {:.8e},{:.8e},{:.8e},{:.8e},\
             {:.8e},{:.8e},{:.8e},\
             {:.8e},{:.8e},\
             {:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},\
             {},{}",
            step, t,
            lx, ly, lz, lvx, lvy, lvz,
            rx, ry, rz, rvx, rvy, rvz,
            length_global, length_mid, delta_mid, strain_axial_mid,
            delta_t_mid_mag, dth_bend_mid_mag, dth_bend_y_mid, dth_twist_mid,
            theta_p_bend_mid_mag, theta_p_bend_y_mid, eps_p_axial_mid,
            theta_max_bend_mid, eps_max_axial_mid,
            r_b, area, iben, k_n, k_bend,
            length0, bond_len_mid0,
            bond_count, bonds_broken,
        ).ok();
    }

    // Profile snapshot: every sample, overwrite `profile.csv` with the
    // current per-atom state. The last-sample file is the final / steady
    // state that the validator plots against the analytical profile.
    let profile_path = format!("{}/data/profile.csv", rec.out_dir);
    if let Ok(f) = File::create(&profile_path) {
        let mut bw = BufWriter::new(f);
        writeln!(bw, "tag,x0,y0,z0,x,y,z,mass").ok();
        let nlocal = atoms.nlocal as usize;
        for i in 0..nlocal {
            let tag = atoms.tag[i] as usize;
            let p0 = rec.initial_pos.get(tag).copied().unwrap_or([0.0; 3]);
            writeln!(
                bw,
                "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
                atoms.tag[i],
                p0[0], p0[1], p0[2],
                atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2],
                atoms.mass[i],
            ).ok();
        }
    }
}
