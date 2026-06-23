//! sphcal_enduring_contact — enduring-contact (rate-independent) stress branch.
//!
//! Tier-4 of the SPH glass-bead calibration. A triperiodic box of glass beads is
//! sheared at a constant rate γ̇ via the native Lees–Edwards (triclinic-box) deform
//! style (`[deform] xy = { style = "erate", rate = γ̇ }`), gravity off, fixed box →
//! homogeneous simple shear. This is the *frictional* dense branch: across a Φ range
//! we measure the full stress tensor and granular temperature, and the sweep forms
//! the enduring-contact residual
//!
//!     σ_contact(Φ) ≈ p_DEM(Φ) − p_KT(Φ, T_measured),
//!
//! the dense-regime gap between the measured DEM pressure and the collisional
//! kinetic-theory (Lun) prediction evaluated at the *measured* T and Φ. It is ≈0
//! below Φ≈0.4 and opens toward jamming — the rate-independent contact-network
//! stress the granular-temperature de-fluidization model consumes.
//!
//! The recorder streams, each thermo interval, the **full stress tensor** (the
//! Love–Weber contact virial `VirialStress / V` plus the kinetic term
//! `Σ m v'⊗v' / V`, with the streaming velocity `v̄(y)=γ̇·y` removed), the
//! pressure `p`, shear stress `σ_xy`, normal-stress differences `N₁,N₂`, the
//! shear-profile-subtracted **granular temperature** `T`, and the solid fraction
//! `Φ` to `<output_dir>/data/enduring_contact_results.csv`. `sweep.py` time-averages
//! the steady window and reports σ_contact(Φ) = p_DEM − p_KT.
//!
//! ```bash
//! cargo run --release --example sphcal_enduring_contact --no-default-features -- examples/SPH_glass_sphere_calibration/04_enduring_contact/config.toml
//! ```

use std::fs;
use std::io::Write as IoWrite;
use std::sync::Once;

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)   // gravity is set to 0 in config; kept for generality
        .add_plugins(FixesPlugin)     // viscous damping for the settle stage
        .add_plugins(DeformPlugin);   // the Lees–Edwards xy shear driver

    app.add_update_system(record_shear, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

/// Record the steady-shear stress tensor, granular temperature, and solid fraction.
fn record_shear(
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

    // Shear rate γ̇ = Δv / L_y (0 during the settle stage → no profile subtraction).
    let ly = domain.size[1];
    let gdot = if ly > 0.0 { domain.boundary_vel[0] / ly } else { 0.0 };
    let y_lo = domain.boundaries_low[1];
    let y_hi = domain.boundaries_high[1];
    let yc = 0.5 * (y_lo + y_hi);

    // ── Kinetic stress + granular temperature (streaming velocity removed) ──
    // v'(i) = v(i) − γ̇·(y_i − y_center) x̂. Kinetic stress σ_kin = Σ m v'⊗v' / V.
    let mut kin = [0.0f64; 6]; // xx, yy, zz, xy, xz, yz
    let mut ke_fluct = 0.0f64; // Σ m |v'|²
    let mut m_total = 0.0f64;
    let mut vol_solid = 0.0f64;
    let dem = registry.get::<DemAtom>();
    for i in 0..nlocal {
        let m = atoms.mass[i];
        m_total += m;
        let vx = atoms.vel[i][0] - gdot * (atoms.pos[i][1] - yc);
        let vy = atoms.vel[i][1];
        let vz = atoms.vel[i][2];
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

    // ── Contact (Love–Weber) virial: σ_contact_ij = −VirialStress_ij / V ──
    let vir = match virial.as_ref() {
        Some(v) => [v.xx, v.yy, v.zz, v.xy, v.xz, v.yz],
        None => [0.0; 6],
    };

    // Reduce across ranks.
    let mut acc = [
        kin[0], kin[1], kin[2], kin[3], kin[4], kin[5],
        vir[0], vir[1], vir[2], vir[3], vir[4], vir[5],
        ke_fluct, m_total, vol_solid,
    ];
    for a in acc.iter_mut() {
        *a = comm.all_reduce_sum_f64(*a);
    }
    let (kin, vir) = ([acc[0], acc[1], acc[2], acc[3], acc[4], acc[5]],
                      [acc[6], acc[7], acc[8], acc[9], acc[10], acc[11]]);
    let ke_fluct = acc[12];
    let m_total = acc[13];
    let vol_solid = acc[14];

    // Total stress σ = kinetic + contact. Contact part flips sign of the virial.
    let mut sig = [0.0f64; 6];
    for k in 0..6 {
        sig[k] = (kin[k] - vir[k]) / vol;
    }
    let p = (sig[0] + sig[1] + sig[2]) / 3.0;
    let n1 = sig[0] - sig[1];
    let n2 = sig[1] - sig[2];
    // von Mises shear stress τ = sqrt(½ σ':σ') from the deviator σ'.
    let dxx = sig[0] - p; let dyy = sig[1] - p; let dzz = sig[2] - p;
    let tau = (0.5 * (dxx * dxx + dyy * dyy + dzz * dzz)
        + sig[3] * sig[3] + sig[4] * sig[4] + sig[5] * sig[5]).sqrt();
    // Granular temperature T = (1/3) Σ m|v'|² / Σ m.
    let t_gran = if m_total > 0.0 { ke_fluct / (3.0 * m_total) } else { 0.0 };
    let phi = vol_solid / vol;

    if comm.rank() != 0 {
        return;
    }

    let time = step as f64 * atoms.dt;
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/SPH_glass_sphere_calibration/04_enduring_contact/data".to_string());
    let path = format!("{}/enduring_contact_results.csv", out_dir);

    static INIT: Once = Once::new();
    INIT.call_once(|| {
        fs::create_dir_all(&out_dir).ok();
        let mut f = fs::File::create(&path).expect("cannot create enduring_contact_results.csv");
        writeln!(
            f,
            "step,time,gdot,sxx,syy,szz,sxy,sxz,syz,p,tau,N1,N2,T,phi"
        )
        .unwrap();
    });

    let mut f = fs::OpenOptions::new().append(true).open(&path).expect("cannot open results csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
        step, time, gdot, sig[0], sig[1], sig[2], sig[3], sig[4], sig[5], p, tau, n1, n2, t_gran, phi
    )
    .unwrap();
}
