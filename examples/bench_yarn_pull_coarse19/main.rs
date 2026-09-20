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
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use dirt_core::dirt_atom::DemAtom;
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

/// Optional relaxation phase (env `YARN_RELAX_STEPS`): grips are freed and the
/// whole system is viscously damped so an un-relaxed packing (e.g. the assembled
/// 4-ply, whose inter-ply contacts overlap) settles before the pull. 0 = none
/// (the pre-relaxed 1-ply needs no relax). Pair with `BOND_BREAK_MIN_STEP` = the
/// same value so the settling transient can't break bonds.
fn relax_steps() -> u64 {
    std::env::var("YARN_RELAX_STEPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}
fn relax_damping() -> f64 {
    std::env::var("YARN_RELAX_DAMPING")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(30.0)
}
/// OVITO-dump cadence (env `DUMP_EVERY`; 0 = off). Under MPI each rank writes its
/// OWNED atoms to `traj_r{rank}.lammpstrj`; merge_traj.py stitches them per frame.
fn dump_every() -> u64 {
    std::env::var("DUMP_EVERY")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

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
) -> &'static (
    std::collections::HashSet<u64>,
    std::collections::HashSet<u64>,
) {
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
    app.add_update_system(relax_damping_system, ParticleSimScheduleSet::PostForce);
    app.add_update_system(grip_zero_force, ParticleSimScheduleSet::PostForce);
    app.add_update_system(measure_pull, ParticleSimScheduleSet::PostForce);
    app.add_update_system(dump_traj, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Set the prescribed grip velocity (±V_PULL on z) on the static grip atoms.
/// During the relax phase the grips are freed (they settle with the rest).
fn grip_control(mut atoms: ResMut<Atom>, run_state: Res<RunState>) {
    if (run_state.total_cycle as u64) < relax_steps() {
        return;
    }
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
/// Skipped during relax (grips are free then).
fn grip_zero_force(mut atoms: ResMut<Atom>, run_state: Res<RunState>) {
    if (run_state.total_cycle as u64) < relax_steps() {
        return;
    }
    let (left, right) = grip_tags(&atoms);
    for i in 0..atoms.nlocal() as usize {
        let t = atoms.tag[i] as u64;
        if left.contains(&t) || right.contains(&t) {
            atoms.force[i] = [0.0; 3];
        }
    }
}

/// Global viscous damping applied to every atom DURING the relax phase only, so
/// an un-relaxed packing settles to equilibrium (critically-damped) before the
/// pull begins. `force -= gamma * v`.
fn relax_damping_system(mut atoms: ResMut<Atom>, run_state: Res<RunState>) {
    if (run_state.total_cycle as u64) >= relax_steps() {
        return;
    }
    let g = relax_damping() as Real;
    for i in 0..atoms.nlocal() as usize {
        atoms.force[i][0] -= g * atoms.vel[i][0];
        atoms.force[i][1] -= g * atoms.vel[i][1];
        atoms.force[i][2] -= g * atoms.vel[i][2];
    }
}

/// Per-rank OVITO/LAMMPS-dump: each rank writes its OWNED atoms (id mol type x y z
/// radius) to `traj_r{rank}.lammpstrj` every `DUMP_EVERY` steps. Column `mol` is
/// the fibre id (mol // 19 = ply) so OVITO can colour by fibre or ply.
/// `merge_traj.py` stitches the per-rank files into one frame-consistent file.
fn dump_traj(
    atoms: Res<Atom>,
    particles: ParticlesWith<'_, Read<DemAtom>>,
    run_state: Res<RunState>,
    input: Res<Input>,
    comm: Res<CommResource>,
) {
    let every = dump_every();
    if every == 0 || (run_state.total_cycle as u64) % every != 0 {
        return;
    }
    let step = run_state.total_cycle;
    let nlocal = atoms.nlocal() as usize;
    let out_dir = input
        .output_dir
        .as_deref()
        .unwrap_or("examples/bench_yarn_pull_coarse19")
        .to_string();
    let path = format!("{}/data/traj_r{}.lammpstrj", out_dir, comm.rank());
    static INIT: OnceLock<()> = OnceLock::new();
    let _ = INIT.get_or_init(|| {
        fs::create_dir_all(format!("{}/data", out_dir)).ok();
        fs::File::create(&path).ok(); // truncate at first frame (per rank/process)
    });
    let mut s = String::with_capacity(nlocal * 60);
    s.push_str(&format!("ITEM: TIMESTEP\n{step}\n"));
    s.push_str(&format!("ITEM: NUMBER OF ATOMS\n{nlocal}\n"));
    s.push_str("ITEM: BOX BOUNDS ff ff ff\n-500 500\n-500 500\n-1500 11500\n");
    s.push_str("ITEM: ATOMS id mol type x y z radius\n");
    particles.with(|dem| {
        for i in 0..nlocal {
            let p = atoms.pos[i];
            let mol = dem.molecule_id[i] as i64;
            s.push_str(&format!(
                "{} {} 1 {:.3} {:.3} {:.3} {:.2}\n",
                atoms.tag[i], mol, p[0], p[1], p[2], dem.radius[i]
            ));
        }
    });
    let mut f = fs::OpenOptions::new().append(true).open(&path).unwrap();
    f.write_all(s.as_bytes()).unwrap();
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
    // Output cadence: coarse (every 20k) until the first bond breaks, then fine
    // (every 1k) to resolve the failure. Reduce the global bond count on EVERY step
    // (collective-safe, no-op in serial) so all ranks agree on the cadence.
    let local_bonds_now = metrics
        .as_deref()
        .map(|m| m.bond_count as f64)
        .unwrap_or(0.0);
    let bonds_now = comm.all_reduce_sum_f64(local_bonds_now) as u64;
    static INIT_BONDS: AtomicU64 = AtomicU64::new(u64::MAX);
    if bonds_now > 0 {
        let _ =
            INIT_BONDS.compare_exchange(u64::MAX, bonds_now, Ordering::Relaxed, Ordering::Relaxed);
    }
    let init = INIT_BONDS.load(Ordering::Relaxed);
    let broken = init != u64::MAX && bonds_now < init;
    // Fine enough to render smooth stress-strain even in the elastic ramp; tighten
    // further once bonds start breaking. (Overridable via YARN_THERMO.)
    let coarse = std::env::var("YARN_THERMO")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(2000);
    let cadence = if broken { coarse.min(1000) } else { coarse };
    if step % cadence != 0 {
        return;
    }
    // NOTE: no early return on rank != 0 — every rank must reach the all_reduce
    // collectives below. Only the file-write / println at the end are rank-0.

    // Strain from the imposed grip displacement; counts from the end of relax.
    let pull_step = (step as u64).saturating_sub(relax_steps());
    let strain = ELONGATION_PER_STEP * pull_step as f64 / GAUGE_LENGTH;

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
    let bonds = bonds_now as usize; // already reduced above

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
