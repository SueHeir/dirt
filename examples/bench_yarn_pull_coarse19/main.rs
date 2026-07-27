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

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin)
        .add_plugins(DemBondPlugin);

    app.add_update_system(measure_pull, ParticleSimScheduleSet::PostForce);

    app.start();
}

/// Record strain, axial stress (Love–Weber contact+bond virial), and the number
/// of surviving bonds every `thermo` steps.
fn measure_pull(
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    virial: Option<Res<VirialStress>>,
    metrics: Option<Res<dirt_core::dirt_bond::BondMetrics>>,
    domain: Res<Domain>,
) {
    let step = run_state.total_cycle;
    if step % 20000 != 0 {
        return;
    }
    if comm.rank() != 0 {
        return;
    }

    let strain = ELONGATION_PER_STEP * step as f64 / GAUGE_LENGTH;

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
        writeln!(f, "step,strain,stress_zz,virial_zz,bond_count").unwrap();
    });
    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("open stress_strain.csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{}",
        step, strain, stress_zz, virial_zz, bonds
    )
    .unwrap();

    println!(
        "[yarn_pull] step={step} strain={:.4}% stress_zz={:.4e} bonds={bonds}",
        strain * 100.0,
        stress_zz
    );
}
