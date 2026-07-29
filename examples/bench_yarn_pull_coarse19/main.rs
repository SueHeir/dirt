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
// Per-grip pull velocity (±vz); 2·V_PULL·dt = 0.000216 elongation/step.
const V_PULL: f64 = 0.057197935013997156;
// Grip z-bounds (LAMMPS `boxleft` / `boxright`) — used ONLY to assign the static
// grip membership once at setup (below), matching LAMMPS's `group ... region`
// which is fixed at setup. No separate relaxation phase is needed: the generated
// packing's t=0 settling transient is identical to LAMMPS's and damps on its own
// within ~20k steps with all bonds intact.
const LEFT_GRIP_ZMAX: f64 = 500.0;
const RIGHT_GRIP_ZMIN: f64 = 9500.0;

/// Static grip membership: the set of atom tags in each z-end grip, captured once
/// from the initial configuration (like LAMMPS's `group left region boxleft`,
/// assigned at setup and fixed thereafter). Re-evaluating membership by *current*
/// z each step would let interior atoms drift across the boundary and be captured
/// into the rigid grip, shortening the effective gauge.
static GRIP_TAGS: std::sync::OnceLock<(
    std::collections::HashSet<u64>,
    std::collections::HashSet<u64>,
)> = std::sync::OnceLock::new();

fn grip_tags(
    atoms: &Atom,
) -> &'static (std::collections::HashSet<u64>, std::collections::HashSet<u64>) {
    GRIP_TAGS.get_or_init(|| {
        let (mut left, mut right) = (
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        );
        for i in 0..atoms.nlocal() as usize {
            let z = atoms.pos[i][2] as f64;
            if z < LEFT_GRIP_ZMAX {
                left.insert(atoms.tag[i] as u64);
            } else if z > RIGHT_GRIP_ZMIN {
                right.insert(atoms.tag[i] as u64);
            }
        }
        (left, right)
    })
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin)
        .add_plugins(DemBondPlugin);

    // Pull the two z-end grips apart at constant velocity, mirroring LAMMPS's
    // `fix move linear` on the static groups: set the grip velocity before the
    // position update, and zero the grip force at PostForce so the velocity-Verlet
    // half-kicks can't let bond/contact forces alter the prescribed motion.
    app.add_update_system(grip_control, ParticleSimScheduleSet::PreInitialIntegration);
    app.add_update_system(grip_zero_force, ParticleSimScheduleSet::PostForce);
    app.add_update_system(measure_pull, ParticleSimScheduleSet::PostForce);

    app.start();
}

/// Set the prescribed grip velocity (±V_PULL on z) on the static grip atoms.
fn grip_control(mut atoms: ResMut<Atom>) {
    let (left, right) = grip_tags(&atoms);
    for i in 0..atoms.nlocal() as usize {
        let t = atoms.tag[i] as u64;
        if left.contains(&t) {
            atoms.vel[i] = [0.0, 0.0, -V_PULL as Real];
        } else if right.contains(&t) {
            atoms.vel[i] = [0.0, 0.0, V_PULL as Real];
        }
    }
}

/// Zero the force on the static grip atoms so neither velocity-Verlet half-kick
/// changes their prescribed velocity — a rigid clamp, not a soft constraint.
fn grip_zero_force(mut atoms: ResMut<Atom>) {
    let (left, right) = grip_tags(&atoms);
    for i in 0..atoms.nlocal() as usize {
        let t = atoms.tag[i] as u64;
        if left.contains(&t) || right.contains(&t) {
            atoms.force[i] = [0.0; 3];
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
    // NOTE: no early return on rank != 0 — every rank must reach the all_reduce
    // collectives below. Only the file-write / println at the end are rank-0.

    // Strain from the imposed grip displacement (move_linear pulls from step 0).
    let strain = ELONGATION_PER_STEP * step as f64 / GAUGE_LENGTH;

    // ── Diagnostics: where does the imposed strain actually go? ──
    // Gauge extension = right-grip mean z − left-grip mean z, vs its initial value.
    let nlocal = atoms.nlocal() as usize;
    let (left, right) = grip_tags(&atoms);
    let (mut lz, mut ln, mut rz, mut rn) = (0.0_f64, 0u32, 0.0_f64, 0u32);
    for i in 0..nlocal {
        let t = atoms.tag[i] as u64;
        if left.contains(&t) {
            lz += atoms.pos[i][2] as f64;
            ln += 1;
        } else if right.contains(&t) {
            rz += atoms.pos[i][2] as f64;
            rn += 1;
        }
    }
    // Reduce grip sums across ranks (grip atoms may live on either rank).
    let lz = comm.all_reduce_sum_f64(lz);
    let rz = comm.all_reduce_sum_f64(rz);
    let ln = comm.all_reduce_sum_f64(ln as f64);
    let rn = comm.all_reduce_sum_f64(rn as f64);
    let gauge = if ln > 0.0 && rn > 0.0 {
        rz / rn - lz / ln
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
    // Reduce bond-strain diagnostics across ranks (max via negated min).
    let sum_eps = comm.all_reduce_sum_f64(sum_eps);
    let n_eps = comm.all_reduce_sum_f64(n_eps as f64);
    let max_eps = -comm.all_reduce_min_f64(-max_eps.abs());
    let mean_bond_eps = if n_eps > 0.0 { sum_eps / n_eps } else { 0.0 };

    // Force-magnitude diagnostic (for contact-force comparison vs LAMMPS). Uses the
    // pre-zero force where possible; grip force is zeroed so interior dominates.
    {
        let (mut fmax, mut fsum, mut fn_gt1) = (0.0_f64, 0.0_f64, 0u32);
        for i in 0..nlocal {
            let f = atoms.force[i];
            let fm = ((f[0] * f[0] + f[1] * f[1] + f[2] * f[2]) as f64).sqrt();
            if fm > fmax {
                fmax = fm;
            }
            fsum += fm;
            if fm > 1.0 {
                fn_gt1 += 1;
            }
        }
        let fmax = -comm.all_reduce_min_f64(-fmax);
        let fsum = comm.all_reduce_sum_f64(fsum);
        let fn_gt1 = comm.all_reduce_sum_f64(fn_gt1 as f64) as u32;
        let ntot = comm.all_reduce_sum_f64(nlocal as f64);
        if comm.rank() == 0 {
            println!(
                "[force] step={step} max|f|={fmax:.4} mean|f|={:.4} n(|f|>1)={fn_gt1}",
                fsum / ntot
            );
        }
    }

    // σ_ij = −VirialStress_ij / V (Love–Weber; pair + bond). Report both the
    // volume-normalised zz pressure and the raw summed virial (≡ LAMMPS summed
    // `compute stress/atom`).
    let vol = domain.size[0] * domain.size[1] * domain.size[2];
    // VirialStress and BondMetrics.bond_count are per-rank local — reduce them.
    let local_vzz = virial.as_deref().map(|v| v.zz).unwrap_or(0.0);
    let virial_zz = -comm.all_reduce_sum_f64(local_vzz);
    let stress_zz = virial_zz / vol;
    let local_bonds = metrics.as_deref().map(|m| m.bond_count as f64).unwrap_or(0.0);
    let bonds = comm.all_reduce_sum_f64(local_bonds) as usize;

    if comm.rank() != 0 {
        return;
    }

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
