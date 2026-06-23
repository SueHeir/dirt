//! sphcal_compressibility — isotropic-compression bulk modulus of glass beads.
//!
//! Deliverable #2 of the SPH glass-sphere calibration campaign: the **bulk
//! compressibility closure** K(Φ) for the MUD SPH model. A triperiodic box of
//! glass beads (gravity off) is compressed quasi-statically and isotropically —
//! the `[deform] vel` driver pushes all three box faces inward at equal speed —
//! while pressure and solid fraction are streamed each thermo interval. The
//! resulting equation of state P(Φ) and its log-slope give the bulk modulus
//!
//!     K(Φ) ≈ ΔP / (ΔΦ/Φ) = dP/d(lnΦ),
//!
//! the EOS slope the SPH solver consumes to close the pressure–density relation.
//!
//! The recorder reads `Res<VirialStress>` (Love–Weber contact virial) and
//! `Res<Domain>`, and writes, each thermo interval, the **pressure**
//! `p = −trace(VirialStress)/(3·V)` plus the small kinetic term (for completeness;
//! it vanishes in the quasi-static limit) and the solid fraction
//! `Φ = Σ (π/6) d³ / V` to `<output_dir>/data/compressibility_results.csv`.
//! `sweep.py` fits P(Φ) over the dense branch and extracts K = dP/d(lnΦ).
//!
//! ```bash
//! cargo run --release --example sphcal_compressibility --no-default-features -- \
//!     examples/SPH_glass_sphere_calibration/02_compressibility/config.toml
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
        .add_plugins(GravityPlugin) // gravity is set to 0 in config; kept for generality
        .add_plugins(FixesPlugin)   // viscous damping to keep the compression quasi-static
        .add_plugins(DeformPlugin); // isotropic inward box compression (vel style on x,y,z)

    app.add_update_system(record_compression, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

/// Record the isotropic pressure and solid fraction each thermo interval.
fn record_compression(
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

    // ── Kinetic stress (diagonal only needed for the trace) and solid volume ──
    // No streaming velocity: isotropic compression has no mean shear profile.
    let mut kin_trace = 0.0f64; // Σ m (vx²+vy²+vz²)
    let mut vol_solid = 0.0f64;
    let dem = registry.get::<DemAtom>();
    for i in 0..nlocal {
        let m = atoms.mass[i];
        let v = atoms.vel[i];
        kin_trace += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if let Some(ref dem) = dem {
            let r = dem.radius[i];
            vol_solid += (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
        }
    }

    // ── Contact (Love–Weber) virial trace: σ_contact = −VirialStress / V ──
    let vir_trace = match virial.as_ref() {
        Some(v) => v.xx + v.yy + v.zz,
        None => 0.0,
    };

    // Reduce across ranks (no-op on a single process).
    let mut acc = [kin_trace, vir_trace, vol_solid];
    for a in acc.iter_mut() {
        *a = comm.all_reduce_sum_f64(*a);
    }
    let (kin_trace, vir_trace, vol_solid) = (acc[0], acc[1], acc[2]);

    // Mean pressure p = trace(σ)/3, σ = (kinetic − contact-virial) / V.
    // p = (Σ m|v|² − trace(VirialStress)) / (3 V).
    let p = (kin_trace - vir_trace) / (3.0 * vol);
    let phi = vol_solid / vol;

    if comm.rank() != 0 {
        return;
    }

    let time = step as f64 * atoms.dt;
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| {
            "examples/SPH_glass_sphere_calibration/02_compressibility/data".to_string()
        });
    let path = format!("{}/compressibility_results.csv", out_dir);

    static INIT: Once = Once::new();
    INIT.call_once(|| {
        fs::create_dir_all(&out_dir).ok();
        let mut f = fs::File::create(&path).expect("cannot create compressibility_results.csv");
        writeln!(f, "step,time,phi,p,volume").unwrap();
    });

    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("cannot open results csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e}",
        step, time, phi, p, vol
    )
    .unwrap();
}
