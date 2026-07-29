//! Coarse-19 yarn tensile pull, sphero-cylinder version. Each fibre (a chain of
//! overlapping spheres in `bench_yarn_pull_coarse19`) is replaced by a chain of
//! single-body capsules bonded COM-to-COM. Geometry from `make_capsule_yarn.py`
//! (data/capsule_yarn_k{K}.txt). Material, bond law, CoR, neighbour cadence, dt
//! and grips match the sphere yarn so the stress-strain curves compare directly.
//!
//! ```bash
//! python3 examples/bench_capsule_yarn_pull/make_capsule_yarn.py 4
//! cargo run --release --example bench_capsule_yarn_pull --no-default-features \
//!     --features precision-double -- examples/bench_capsule_yarn_pull/config.toml
//! ```

use std::fs;
use std::io::Write;
use std::sync::OnceLock;

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use dirt_core::soil_core::{BondEntry, BondStore};

const GAUGE_LENGTH: f64 = 9000.0;
const ELONGATION_PER_STEP: f64 = 0.000216;
const V_PULL: f64 = 0.057197935013997156;
const LEFT_GRIP_ZMAX: f64 = 500.0;
const RIGHT_GRIP_ZMIN: f64 = 9500.0;
// Capsule density = the sphere yarn's `set density` value (mass only affects the
// dynamic transient; the quasi-static curve is rate-controlled, density-independent).
const DENSITY: f64 = 9.93245e-07;

fn spec_path() -> String {
    std::env::var("CAPSULE_YARN_SPEC").unwrap_or_else(|_| {
        "examples/bench_capsule_yarn_pull/data/capsule_yarn_k4.txt".to_string()
    })
}

/// Optional: limit to the first N fibres (env FIBRE_LIMIT) for a fast stability check.
fn fibre_limit() -> u32 {
    std::env::var("FIBRE_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(u32::MAX)
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(dirt_core::soil_verlet::VelocityVerletPlugin::new())
        .add_plugins(CapsuleContactPlugin)
        .add_plugins(CapsulePlugin)
        .add_plugins(DemBondPlugin);

    // Build the capsule yarn (insert capsules + COM-to-COM bonds) at the insertion
    // stage, before the bond history / ghost-cutoff setup runs at PostSetup.
    app.add_setup_system(build_capsule_yarn, ScheduleSetupSet::Setup);

    // Rigid z-end grips + stress-strain measurement (same as the sphere yarn).
    app.add_update_system(grip_control, ParticleSimScheduleSet::PreInitialIntegration);
    app.add_update_system(grip_zero_force, ParticleSimScheduleSet::PostForce);
    app.add_update_system(measure_pull, ParticleSimScheduleSet::PostForce);

    app.start();
}

/// One capsule read from the spec file.
struct SpecCapsule {
    com: [f64; 3],
    quat: [f64; 4],
    half_length: f64,
    molecule: u32,
}

/// Insert the capsule yarn from the plain-text spec and bond consecutive capsules
/// within each fibre (COM-to-COM, r0 = the COM distance).
fn build_capsule_yarn(
    comm: Res<CommResource>,
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    scheduler_manager: Res<SchedulerManager>,
) {
    if scheduler_manager.index != 0 || comm.rank() != 0 {
        return;
    }
    let text = fs::read_to_string(spec_path()).expect("read capsule yarn spec");
    let mut lines = text.lines();
    let hdr: Vec<f64> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    let (radius, n_cap, n_bond) = (hdr[0], hdr[2] as usize, hdr[3] as usize);

    let limit = fibre_limit();
    let mut caps: Vec<SpecCapsule> = Vec::with_capacity(n_cap);
    for _ in 0..n_cap {
        let f: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        let mol: u32 = f[8].parse().unwrap();
        if mol > limit {
            continue;
        }
        caps.push(SpecCapsule {
            com: [f[0].parse().unwrap(), f[1].parse().unwrap(), f[2].parse().unwrap()],
            quat: [
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
            ],
            half_length: f[7].parse().unwrap(),
            molecule: mol,
        });
    }
    // Bonds reference original spec indices; keep a map old_index -> new_index so
    // the FIBRE_LIMIT subset stays consistent.
    let mut keep: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut raw = text.lines();
    raw.next(); // header
    let mut new_i = 0usize;
    for old_i in 0..n_cap {
        let mol: u32 = raw.next().unwrap().split_whitespace().nth(8).unwrap().parse().unwrap();
        if mol <= limit {
            keep.insert(old_i, new_i);
            new_i += 1;
        }
    }

    // ── insert capsules ──
    for c in &caps {
        let def = CapsuleDef {
            name: "yarn_capsule".to_string(),
            radius,
            length: 2.0 * c.half_length,
        };
        // atom_type 0 = the single [[dem.materials]] entry (material table is 0-indexed).
        insert_capsule(&mut atoms, &registry, &def, c.com, [0.0; 3], c.quat, DENSITY, 0);
    }
    // molecule id (intra-fibre contact exclusion, mirrors the sphere yarn)
    {
        let mut dem = registry.expect_mut::<DemAtom>("build_capsule_yarn");
        for (i, c) in caps.iter().enumerate() {
            dem.molecule_id[i] = c.molecule as f64;
        }
    }
    // bonds (COM-to-COM); r0 from the inserted positions
    {
        let mut bonds = registry.expect_mut::<BondStore>("build_capsule_yarn");
        while bonds.bonds.len() < atoms.len() {
            bonds.bonds.push(Vec::new());
        }
        for _ in 0..n_bond {
            let f: Vec<usize> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|s| s.parse::<i64>().unwrap() as usize)
                .collect();
            let (oi, oj, bt) = (f[0], f[1], f[2] as u32);
            let (i, j) = match (keep.get(&oi), keep.get(&oj)) {
                (Some(&i), Some(&j)) => (i, j),
                _ => continue,
            };
            let pi = atoms.pos[i];
            let pj = atoms.pos[j];
            let r0 = (((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2))
                as f64)
                .sqrt();
            let (ti, tj) = (atoms.tag[i], atoms.tag[j]);
            bonds.bonds[i].push(BondEntry { partner_tag: tj, bond_type: bt, r0 });
            bonds.bonds[j].push(BondEntry { partner_tag: ti, bond_type: bt, r0 });
        }
    }
    println!(
        "build_capsule_yarn: inserted {} capsules, {} bonds (radius {:.2})",
        caps.len(),
        n_bond,
        radius
    );
}

/// Static grip membership captured once from initial z (like the sphere yarn).
static GRIP_TAGS: OnceLock<(std::collections::HashSet<u64>, std::collections::HashSet<u64>)> =
    OnceLock::new();

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

fn grip_zero_force(mut atoms: ResMut<Atom>) {
    let (left, right) = grip_tags(&atoms);
    for i in 0..atoms.nlocal() as usize {
        let t = atoms.tag[i] as u64;
        if left.contains(&t) || right.contains(&t) {
            atoms.force[i] = [0.0; 3];
        }
    }
}

/// Record strain, axial stress (virial), and surviving-bond count every `thermo`.
fn measure_pull(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    virial: Option<Res<VirialStress>>,
    metrics: Option<Res<dirt_core::dirt_bond::BondMetrics>>,
    domain: Res<Domain>,
) {
    let step = run_state.total_cycle;
    if step % 20000 != 0 || comm.rank() != 0 {
        return;
    }
    let strain = ELONGATION_PER_STEP * step as f64 / GAUGE_LENGTH;

    let nlocal = atoms.nlocal() as usize;
    let (left, right) = grip_tags(&atoms);
    let (mut lz, mut ln, mut rz, mut rn) = (0.0_f64, 0u32, 0.0_f64, 0u32);
    let (mut fmax, mut fgt1) = (0.0_f64, 0u32);
    for i in 0..nlocal {
        let t = atoms.tag[i] as u64;
        let z = atoms.pos[i][2] as f64;
        if left.contains(&t) {
            lz += z;
            ln += 1;
        } else if right.contains(&t) {
            rz += z;
            rn += 1;
        }
        let f = atoms.force[i];
        let fm = ((f[0] * f[0] + f[1] * f[1] + f[2] * f[2]) as f64).sqrt();
        if fm > fmax {
            fmax = fm;
        }
        if fm > 1.0 {
            fgt1 += 1;
        }
    }
    let gauge = if ln > 0 && rn > 0 {
        rz / rn as f64 - lz / ln as f64
    } else {
        0.0
    };

    let vol = domain.size[0] * domain.size[1] * domain.size[2];
    let stress_zz = match virial.as_deref() {
        Some(v) => -v.zz / vol,
        None => 0.0,
    };
    let bonds = metrics.as_deref().map(|m| m.bond_count).unwrap_or(0);

    let out_dir = input
        .output_dir
        .as_deref()
        .unwrap_or("examples/bench_capsule_yarn_pull")
        .to_string();
    let path = format!("{}/data/stress_strain.csv", out_dir);
    static INIT: OnceLock<()> = OnceLock::new();
    let _ = INIT.get_or_init(|| {
        fs::create_dir_all(format!("{}/data", out_dir)).ok();
        let mut f = fs::File::create(&path).expect("create stress_strain.csv");
        writeln!(f, "step,strain,stress_zz,bond_count,gauge").unwrap();
    });
    let mut f = fs::OpenOptions::new().append(true).open(&path).unwrap();
    writeln!(f, "{},{:.8e},{:.8e},{},{:.6e}", step, strain, stress_zz, bonds, gauge).unwrap();

    println!(
        "[capsule_yarn] step={step} strain={:.3}% stress={stress_zz:.3e} bonds={bonds} \
         gauge={gauge:.1} max|f|={fmax:.4} n(|f|>1)={fgt1}",
        strain * 100.0
    );
}
