//! Coarse-19 yarn tensile pull — DIRT reproduction of the LAMMPS `bpm/guo` run.
//!
//! Loads the 1-ply coarse-19 yarn (5330 bonded spheres, 2 bond types) from the
//! LAMMPS data file, pulls the two z-end grips apart at constant velocity, and
//! writes the stress-strain curve and the surviving-bond count. Compare
//! `data/stress_strain.csv` against the LAMMPS
//! `..._pull_atomstress-strain.dat` (strain, stress_zz, nbond).
//!
//! ```bash
//! cargo run --release --example bench_yarn_pull_coarse19 --no-default-features \
//!     --features precision-double -- examples/bench_yarn_pull_coarse19/config.toml
//! ```

use std::fs;
use std::io::Write;
use std::sync::OnceLock;

use dirt_core::prelude::*;

// Gauge length (LAMMPS `boxmiddle` z-span) and per-step elongation, so strain is
// reported identically to the LAMMPS run (`elongation = 0.000216 * step`).
const GAUGE_LENGTH: f64 = 9000.0;
const ELONGATION_PER_STEP: f64 = 0.000216;
// Per-grip pull velocity (±vz), matching the LAMMPS `fix move linear`.
const V_PULL: f64 = 0.057197935013997156;
// Grip z-bounds (LAMMPS `boxleft` / `boxright`).
const LEFT_GRIP_ZMAX: f64 = 500.0;
const RIGHT_GRIP_ZMIN: f64 = 9500.0;

// Number of relaxation steps before the pull begins. The generated yarn is a
// compressed packing (~0.12 inter-fibre overlaps); DIRT under-damps that
// relaxation, so the grips are held fixed for this many steps to let the packing
// settle (KE → ~0) before pulling — otherwise residual relaxation energy breaks
// bonds prematurely. Overridable with YARN_RELAX_STEPS.
fn relax_steps() -> u64 {
    std::env::var("YARN_RELAX_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(150_000)
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin)
        .add_plugins(DemBondPlugin);

    // Grip control mirrors `move_linear`: set the grip velocity before position
    // integration AND zero the grip force/torque at PostForce, so neither
    // velocity-Verlet half-kick lets the restraining bond forces alter the
    // prescribed grip motion (a rigid grip, not a soft velocity constraint).
    app.add_update_system(grip_control, ParticleSimScheduleSet::PreInitialIntegration);
    app.add_update_system(grip_zero_force, ParticleSimScheduleSet::PostForce);
    app.add_update_system(measure_pull, ParticleSimScheduleSet::PostForce);

    app.start();
}

/// Zero the translational force on the grip atoms so the velocity-Verlet
/// half-kicks do not change their prescribed velocity — makes the grips rigid
/// clamps that deliver the full imposed displacement into the gauge section.
fn grip_zero_force(mut atoms: ResMut<Atom>) {
    let nlocal = atoms.nlocal() as usize;
    for i in 0..nlocal {
        let z = atoms.pos[i][2] as f64;
        if z < LEFT_GRIP_ZMAX || z > RIGHT_GRIP_ZMIN {
            atoms.force[i] = [0.0; 3];
        }
    }
}

/// Freeze the z-end grips during relaxation (`step < relax_steps`), then pull
/// them apart at constant velocity. Overrides the grip velocities each step, so
/// no `[[move_linear]]` fix is needed.
fn grip_control(mut atoms: ResMut<Atom>, run_state: Res<RunState>) {
    let pulling = run_state.total_cycle as u64 >= relax_steps();
    let nlocal = atoms.nlocal() as usize;
    for i in 0..nlocal {
        let z = atoms.pos[i][2] as f64;
        if z < LEFT_GRIP_ZMAX {
            atoms.vel[i] = [0.0, 0.0, if pulling { -V_PULL as Real } else { 0.0 }];
        } else if z > RIGHT_GRIP_ZMIN {
            atoms.vel[i] = [0.0, 0.0, if pulling { V_PULL as Real } else { 0.0 }];
        }
    }
}

/// Record strain, axial stress (Love–Weber contact+bond virial), and the number
/// of surviving bonds every `thermo` steps.
fn measure_pull(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    virial: Option<Res<VirialStress>>,
    metrics: Option<Res<dirt_core::dirt_bond::BondMetrics>>,
    domain: Res<Domain>,
    particles: ParticlesWith<'_, Read<dirt_core::soil_core::BondStore>>,
) {
    let step = run_state.total_cycle;
    if step % 20000 != 0 {
        return;
    }
    if comm.rank() != 0 {
        return;
    }

    // Strain counts from the start of the pull (after the relaxation phase).
    let relax = relax_steps();
    let pull_step = (step as u64).saturating_sub(relax);
    let strain = ELONGATION_PER_STEP * pull_step as f64 / GAUGE_LENGTH;

    // ── Diagnostics: where does the imposed strain actually go? ──
    // Gauge extension = right-grip mean z − left-grip mean z, vs its initial value.
    let nlocal = atoms.nlocal() as usize;
    let (mut lz, mut ln, mut rz, mut rn) = (0.0_f64, 0u32, 0.0_f64, 0u32);
    for i in 0..nlocal {
        let z = atoms.pos[i][2] as f64;
        if z < LEFT_GRIP_ZMAX {
            lz += z;
            ln += 1;
        } else if z > RIGHT_GRIP_ZMIN {
            rz += z;
            rn += 1;
        }
    }
    let gauge = if ln > 0 && rn > 0 {
        rz / rn as f64 - lz / ln as f64
    } else {
        0.0
    };
    // Mean/max axial bond strain (dist/r0 − 1), from BondStore + a tag→index map.
    let mut tag_idx = std::collections::HashMap::with_capacity(nlocal);
    for i in 0..nlocal {
        tag_idx.insert(atoms.tag[i], i);
    }
    let (mut sum_eps, mut n_eps, mut max_eps) = (0.0_f64, 0u64, 0.0_f64);
    particles.with(|bonds| {
        for i in 0..nlocal.min(bonds.bonds.len()) {
            let pi = atoms.pos[i];
            for b in &bonds.bonds[i] {
                if let Some(&j) = tag_idx.get(&b.partner_tag) {
                    if j <= i {
                        continue; // count each bond once
                    }
                    let pj = atoms.pos[j];
                    let dist = (((pi[0] - pj[0]).powi(2)
                        + (pi[1] - pj[1]).powi(2)
                        + (pi[2] - pj[2]).powi(2)) as f64)
                        .sqrt();
                    if b.r0 > 0.0 {
                        let eps = dist / b.r0 - 1.0;
                        sum_eps += eps;
                        n_eps += 1;
                        if eps.abs() > max_eps.abs() {
                            max_eps = eps;
                        }
                    }
                }
            }
        }
    });
    let mean_bond_eps = if n_eps > 0 { sum_eps / n_eps as f64 } else { 0.0 };

    // σ_ij = −VirialStress_ij / V (Love–Weber; pair + bond). Report both the
    // volume-normalised zz pressure and the raw summed virial (≡ LAMMPS summed
    // `compute stress/atom`).
    let vol = domain.size[0] * domain.size[1] * domain.size[2];
    let (stress_zz, virial_zz) = match virial.as_deref() {
        Some(v) => (-v.zz / vol, -v.zz),
        None => (0.0, 0.0),
    };
    let bonds = metrics.as_deref().map(|m| m.bond_count).unwrap_or(0);

    let out_dir = input
        .output_dir
        .as_deref()
        .unwrap_or("examples/bench_yarn_pull_coarse19")
        .to_string();
    let path = format!("{}/data/stress_strain.csv", out_dir);
    static INIT: OnceLock<()> = OnceLock::new();
    let _ = INIT.get_or_init(|| {
        fs::create_dir_all(format!("{}/data", out_dir)).ok();
        let mut f = fs::File::create(&path).expect("create stress_strain.csv");
        writeln!(
            f,
            "step,strain,stress_zz,virial_zz,bond_count,gauge,mean_bond_eps,max_bond_eps"
        )
        .unwrap();
    });
    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("open stress_strain.csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{},{:.6e},{:.6e},{:.6e}",
        step, strain, stress_zz, virial_zz, bonds, gauge, mean_bond_eps, max_eps
    )
    .unwrap();

    println!(
        "[yarn_pull] step={step} strain={:.3}% stress={:.3e} bonds={bonds} gauge={gauge:.1} mean_bond_eps={mean_bond_eps:.4e} max_bond_eps={max_eps:.4e}",
        strain * 100.0,
        stress_zz
    );
}
