//! bench_cundall_damping — validates the Cundall **non-viscous** global damping
//! fix (`[[cundall]]`, `dirt_fixes`) against its exact, analytically-known
//! effect on a single free particle, and cross-checks the linear part against
//! LAMMPS `fix damping/cundall`.
//!
//! # What is tested
//!
//! A single sphere lives in an empty box (no walls, no contacts). Two
//! independent, sign-flipping motions run simultaneously so that BOTH branches
//! of the Cundall sign function are exercised for BOTH the linear and angular
//! coefficients:
//!
//! * **Linear (γ_l):** the sphere is launched straight up (+z) under gravity.
//!   Gravity is applied in the Force phase; the `[[cundall]]` fix (PostForce)
//!   then damps the *total* force component-by-component exactly as LAMMPS does:
//!
//!   ```text
//!   F_z <- F_z * (1 - γ_l * sign(F_z · v_z))
//!   ```
//!
//!   While rising, F_z (= −mg) opposes v_z (> 0): sign < 0, so the force is
//!   *amplified* to −mg(1+γ_l) ⇒ constant acceleration a_up = −g(1+γ_l).
//!   While falling, F_z acts along v_z (< 0): sign > 0, so the force is *reduced*
//!   to −mg(1−γ_l) ⇒ a_down = −g(1−γ_l). Both are exact, mass-independent, and
//!   piecewise constant — a clean analytical target.
//!
//! * **Angular (γ_a):** the sphere is given an initial spin (+z) and a constant
//!   applied torque `TORQUE_Z` about z (opposing the spin). The `[[cundall]]`
//!   fix damps the torque by `1 − γ_a·sign(T_z·ω_z)`, giving spin-down
//!   acceleration α_down = T_z·(1+γ_a)/I while ω_z > 0 and α_up = T_z·(1−γ_a)/I
//!   after ω_z crosses zero. (There is no "angular gravity" fix, so the constant
//!   torque is applied here in the example — it stays example-local.)
//!
//! The recorder logs `t, v_z, omega_z` each step and writes a self-describing
//! CSV whose header carries the *actual* run parameters (g, I⁻¹, applied torque,
//! and the γ_l/γ_a read back from the live `[[cundall]]` fix), so `sweep.py`
//! validates against theory with zero duplicated constants.
//!
//! ```bash
//! cargo run --release --example bench_cundall_damping --no-default-features \
//!   --features precision-double -- examples/bench_cundall_damping/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Initial angular velocity about +z (rad/s), set once at step 0.
const OMEGA0_Z: f64 = 5.0;
/// Constant applied torque about z (N·m), opposing the initial spin. Applied
/// every step in the Force phase so the `[[cundall]]` fix (PostForce) damps it.
const TORQUE_Z: f64 = -1.0e-6;

/// Records the free particle's linear/angular state each step and writes a CSV.
struct CundallTracker {
    tag: Option<u32>,
    spin_set: bool,
    rows: Vec<(f64, f64, f64)>, // (t, v_z, omega_z)
    written: bool,
}

impl CundallTracker {
    fn new() -> Self {
        Self { tag: None, spin_set: false, rows: Vec::new(), written: false }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins) // DEM atom + rotational Verlet
        .add_plugins(GravityPlugin) // body force F = m g (linear driver)
        .add_plugins(FixesPlugin); // provides the [[cundall]] fix under test

    app.add_resource(CundallTracker::new());

    // Set the initial spin once, and apply the constant torque every step,
    // both before the Cundall fix runs (PostForce).
    app.add_update_system(set_initial_spin, ParticleSimScheduleSet::PreForce);
    app.add_update_system(apply_constant_torque, ParticleSimScheduleSet::Force);
    app.add_update_system(record, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// One-shot: give the single sphere its initial angular velocity about +z.
fn set_initial_spin(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    mut tracker: ResMut<CundallTracker>,
) {
    if tracker.spin_set || atoms.nlocal < 1 {
        return;
    }
    let mut dem = registry.expect_mut::<DemAtom>("set_initial_spin");
    dem.omega[0] = [0.0, 0.0, OMEGA0_Z];
    tracker.tag = Some(atoms.tag[0]);
    tracker.spin_set = true;
}

/// Apply the constant driving torque about z every step (Force phase), so the
/// Cundall fix damps a non-zero torque and the angular decay is well-defined.
fn apply_constant_torque(atoms: Res<Atom>, registry: Res<AtomDataRegistry>) {
    if atoms.nlocal < 1 {
        return;
    }
    let mut dem = registry.expect_mut::<DemAtom>("apply_constant_torque");
    dem.torque[0][2] += TORQUE_Z;
}

/// Record (t, v_z, omega_z) each step; flush a self-describing CSV at the end.
fn record(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    input: Res<Input>,
    gravity: Res<GravityConfig>,
    fixes: Res<FixesRegistry>,
    mut tracker: ResMut<CundallTracker>,
) {
    if tracker.written || atoms.nlocal == 0 {
        return;
    }
    let tag = match tracker.tag {
        Some(t) => t,
        None => return,
    };
    let dem = registry.expect::<DemAtom>("record");
    let i = match (0..atoms.nlocal as usize).find(|&k| atoms.tag[k] == tag) {
        Some(k) => k,
        None => return,
    };

    let dt = atoms.dt;
    let step = run_state.total_cycle;
    let t = step as f64 * dt;
    tracker.rows.push((t, atoms.vel[i][2] as f64, dem.omega[i][2]));

    let last_step = run_config.current_stage(0).steps as usize;
    if step + 1 >= last_step {
        // Read the live Cundall coefficients back from the fix registry so the
        // CSV is fully self-describing (no constants duplicated in sweep.py).
        let (gamma_l, gamma_a) = fixes
            .cundall
            .first()
            .map(|c| (c.gamma_l, c.gamma_a))
            .unwrap_or((0.0, 0.0));
        let inv_inertia = dem.inv_inertia[i];

        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_cundall_damping".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/cundall_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "# dt={:.10e} g={:.10e} inv_inertia={:.10e} torque_z={:.10e} gamma_l={:.10e} gamma_a={:.10e}",
            dt, -gravity.gz, inv_inertia, TORQUE_Z, gamma_l, gamma_a
        )
        .unwrap();
        writeln!(f, "t,vz,omega_z").unwrap();
        for (t, vz, wz) in &tracker.rows {
            writeln!(f, "{:.10e},{:.10e},{:.10e}", t, vz, wz).unwrap();
        }
        tracker.written = true;
        println!("=== Cundall Damping Results ===");
        println!("  rows recorded:    {}", tracker.rows.len());
        println!("  gamma_l={:.4}  gamma_a={:.4}", gamma_l, gamma_a);
        println!("  final v_z:        {:.6e} m/s", atoms.vel[i][2] as f64);
        println!("  final omega_z:    {:.6e} rad/s", dem.omega[i][2]);
        println!("  results saved to: {}", results_file);
    }
}
