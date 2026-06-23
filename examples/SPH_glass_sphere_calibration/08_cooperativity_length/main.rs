//! sphcal_cooperativity_length — DEM measurement of the granular
//! **cooperativity length** ξ for the MUD SPH model's nonlocal closure.
//!
//! MUD's nonlocal-granular-fluidity (NGF) branch keeps one idea from Henann–Kamrin
//! inside the granular-temperature model: a cooperativity length
//! `ξ(μ) = A·d / √|μ − μ_s|` that diverges at yield, *driving the contact branch to
//! creep below yield*. This rig measures the two closures that branch needs:
//!
//!  1. **ξ vs μ** — the amplitude `A` and the divergence as μ → μ_s.
//!  2. **the Zhang–Kamrin bridge `g ∝ √T`** — confirming the fluidity field is a
//!     measure of velocity fluctuations, so MUD can drive it from temperature.
//!
//! Same Lees–Edwards homogeneous-shear rig as `bench_lebc_shear` (so μ, p, T, Φ, I
//! are measured identically), plus a **spatial velocity-fluctuation correlation**
//! along the vorticity (z) direction:
//!   C(Δz) = ⟨δv'(0)·δv'(Δz)⟩ / ⟨|δv'|²⟩,   δv' = v − γ̇(y−y_c) x̂,
//! restricted to near-columns (|Δx|,|Δy| < d) so the Lees–Edwards xy tilt is
//! irrelevant. The correlation length ξ is the integral of C(Δz) up to its first
//! zero crossing. `sweep.py` runs a γ̇ grid, time-averages the steady window, and
//! fits `ξ(μ) = A d/√(μ−μ_s)` and `g = γ̇/μ ∝ √T`.
//!
//! ```bash
//! cargo run --release --example sphcal_cooperativity_length --no-default-features -- \
//!     examples/SPH_glass_sphere_calibration/08_cooperativity_length/config.toml
//! ```
//! Single-rank (the pair correlation is computed on local atoms with minimum image).

use std::fs;
use std::io::Write as IoWrite;
use std::sync::Once;

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;

/// Grain diameter d [m] (mean of the polydisperse pack) — sets the correlation
/// bin width and the inertial number.
const D_GRAIN: f64 = 0.5e-3;
/// Solid-grain density ρ_s [kg/m³] (for the inertial number).
const RHO_S: f64 = 2500.0;
/// Correlation bins along z: width Δz = 0.5 d, up to 6 d (< L/2 ≈ 7.1 d).
const NBINS: usize = 12;
const DR: f64 = 0.5 * D_GRAIN;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(FixesPlugin)
        .add_plugins(DeformPlugin);

    app.add_update_system(record_cooperativity, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

/// Record μ, p, T, Φ, I and the velocity-fluctuation correlation length ξ.
fn record_cooperativity(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    domain: Res<Domain>,
    virial: Option<Res<VirialStress>>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
) {
    let step = run_state.total_cycle;
    if step % 2000 != 0 {
        return;
    }
    let nlocal = atoms.nlocal as usize;
    let vol = domain.volume;
    if vol <= 0.0 {
        return;
    }

    let ly = domain.size[1];
    let gdot = if ly > 0.0 { domain.boundary_vel[0] / ly } else { 0.0 };
    let yc = 0.5 * (domain.boundaries_low[1] + domain.boundaries_high[1]);

    // ── Per-particle fluctuation velocity δv' = v − γ̇(y−y_c) x̂ ──
    let mut dv = vec![[0.0f64; 3]; nlocal];
    let mut kin = [0.0f64; 6];
    let mut ke_fluct = 0.0f64;
    let mut m_total = 0.0f64;
    let mut vol_solid = 0.0f64;
    let dem = registry.get::<DemAtom>();
    for i in 0..nlocal {
        let m = atoms.mass[i];
        m_total += m;
        let vx = atoms.vel[i][0] - gdot * (atoms.pos[i][1] - yc);
        let vy = atoms.vel[i][1];
        let vz = atoms.vel[i][2];
        dv[i] = [vx, vy, vz];
        kin[0] += m * vx * vx;
        kin[1] += m * vy * vy;
        kin[2] += m * vz * vz;
        kin[3] += m * vx * vy;
        kin[4] += m * vx * vz;
        kin[5] += m * vy * vz;
        ke_fluct += m * (vx * vx + vy * vy + vz * vz);
        if let Some(ref dem) = dem {
            let r = dem.radius[i];
            vol_solid += (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
        }
    }

    // ── Stress tensor → p, τ, μ (identical to bench_lebc_shear) ──
    let vir = match virial.as_ref() {
        Some(v) => [v.xx, v.yy, v.zz, v.xy, v.xz, v.yz],
        None => [0.0; 6],
    };
    let mut acc = [
        kin[0], kin[1], kin[2], kin[3], kin[4], kin[5],
        vir[0], vir[1], vir[2], vir[3], vir[4], vir[5],
        ke_fluct, m_total, vol_solid,
    ];
    for a in acc.iter_mut() {
        *a = comm.all_reduce_sum_f64(*a);
    }
    let kin = [acc[0], acc[1], acc[2], acc[3], acc[4], acc[5]];
    let vir = [acc[6], acc[7], acc[8], acc[9], acc[10], acc[11]];
    let ke_fluct = acc[12];
    let m_total = acc[13];
    let vol_solid = acc[14];

    let mut sig = [0.0f64; 6];
    for k in 0..6 {
        sig[k] = (kin[k] - vir[k]) / vol;
    }
    let p = (sig[0] + sig[1] + sig[2]) / 3.0;
    let dxx = sig[0] - p;
    let dyy = sig[1] - p;
    let dzz = sig[2] - p;
    let tau = (0.5 * (dxx * dxx + dyy * dyy + dzz * dzz)
        + sig[3] * sig[3] + sig[4] * sig[4] + sig[5] * sig[5])
        .sqrt();
    let t_gran = if m_total > 0.0 { ke_fluct / (3.0 * m_total) } else { 0.0 };
    let phi = vol_solid / vol;
    let mu = if p > 0.0 { tau / p } else { 0.0 };
    let inertial = if p > 0.0 { gdot * D_GRAIN / (p / RHO_S).sqrt() } else { 0.0 };
    let g_fluidity = if mu > 0.0 { gdot / mu } else { 0.0 }; // NGF fluidity g = γ̇/μ

    // ── Spatial velocity-fluctuation correlation along z (vorticity) ──
    // C(Δz) over near-columns (|Δx|,|Δy| < d, minimum image), so the LE tilt
    // (in xy) does not enter. mean_dv2 = ⟨|δv'|²⟩ normalizes to C(0)=1.
    let (lx, lz) = (domain.size[0], domain.size[2]);
    let col = D_GRAIN; // column half-width
    let mut corr = [0.0f64; NBINS];
    let mut cnt = [0u64; NBINS];
    let mut sum_dv2 = 0.0f64;
    for i in 0..nlocal {
        sum_dv2 += dv[i][0] * dv[i][0] + dv[i][1] * dv[i][1] + dv[i][2] * dv[i][2];
    }
    let mean_dv2 = if nlocal > 0 { sum_dv2 / nlocal as f64 } else { 0.0 };
    let minimg = |d: f64, l: f64| -> f64 {
        if l > 0.0 { d - l * (d / l).round() } else { d }
    };
    for i in 0..nlocal {
        for j in (i + 1)..nlocal {
            let dyij = minimg(atoms.pos[i][1] - atoms.pos[j][1], ly);
            if dyij.abs() > col {
                continue;
            }
            let dxij = minimg(atoms.pos[i][0] - atoms.pos[j][0], lx);
            if dxij.abs() > col {
                continue;
            }
            let dzij = minimg(atoms.pos[i][2] - atoms.pos[j][2], lz).abs();
            let k = (dzij / DR) as usize;
            if k < NBINS {
                corr[k] += dv[i][0] * dv[j][0] + dv[i][1] * dv[j][1] + dv[i][2] * dv[j][2];
                cnt[k] += 1;
            }
        }
    }
    // Normalize each bin by ⟨|δv'|²⟩, then the **integral correlation length**
    // ξ = (∫ C dΔz) / C(0) over the leading positive run — a decay *range*
    // independent of the (amplitude) value of C at contact. For an exponential
    // C(Δz)=C₀e^{−Δz/ξ} this returns ξ exactly.
    let mut cvals = [0.0f64; NBINS];
    let mut integral = 0.0f64;
    let mut still_positive = true;
    for k in 0..NBINS {
        let c = if cnt[k] > 0 && mean_dv2 > 0.0 {
            corr[k] / cnt[k] as f64 / mean_dv2
        } else {
            0.0
        };
        cvals[k] = c;
        if still_positive {
            if c > 0.0 {
                integral += c * DR;
            } else {
                still_positive = false;
            }
        }
    }
    let xi = if cvals[0] > 0.0 { integral / cvals[0] } else { 0.0 };

    if comm.rank() != 0 {
        return;
    }
    let time = step as f64 * atoms.dt;
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/SPH_glass_sphere_calibration/08_cooperativity_length/data".to_string());
    let path = format!("{}/cooperativity_results.csv", out_dir);
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        fs::create_dir_all(&out_dir).ok();
        let mut f = fs::File::create(&path).expect("cannot create cooperativity_results.csv");
        let cols: Vec<String> = (0..NBINS).map(|k| format!("C{k}")).collect();
        writeln!(
            f,
            "step,time,gdot,I,mu,p,tau,T,phi,g,xi,dr,{}",
            cols.join(",")
        )
        .unwrap();
    });
    let mut f = fs::OpenOptions::new().append(true).open(&path).expect("open csv");
    let cstr: Vec<String> = cvals.iter().map(|c| format!("{c:.6e}")).collect();
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{}",
        step, time, gdot, inertial, mu, p, tau, t_gran, phi, g_fluidity, xi, DR, cstr.join(",")
    )
    .unwrap();
}
