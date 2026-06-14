//! bench_rolling_decay — validates the constant-torque rolling-resistance
//! model against the exact pure-rolling deceleration for that model.
//!
//! A single sphere is launched horizontally with a matching spin (v = ωR, pure
//! rolling) on a flat frictional floor `[[wall]]` (z = 0, normal +z), under
//! gravity. With rolling friction μ_r the wall contact applies a decelerating
//! rolling-resistance couple τ_r = μ_r·F_n·r_eff that opposes the spin, while
//! Mindlin static (sliding) friction enforces the rolling constraint. Because a
//! flat wall has r_eff = R (no curvature correction), the sphere decelerates at
//! the EXACT constant rate
//!
//!     a = (5/7) · μ_r · g
//!
//! (derived in README.md: rolling-resistance couple τ_r = μ_r·F_n·R + Mindlin
//! static friction enforcing v = ωR, with I = 2/5 m R²; the inertia enters as
//! I/R + mR = 7/5 mR, hence 5/7). The recorder logs t, x, vx, omega so
//! `sweep.py` can fit the slope and check that the sphere stays in pure rolling
//! (vx ≈ ωR) throughout.
//!
//! `dirt_wall` now carries the full friction trio on every wall type: normal,
//! Mindlin sliding friction (material `friction`/`friction_ij`), and rolling
//! resistance (material `rolling_friction`/`rolling_friction_ij`, `constant`
//! and `sds` models per `rolling_model`). So the floor is a real wall plane —
//! no giant frozen sphere is needed.
//!
//! The moving sphere has no initial-spin TOML knob (`[[particles.insert]]`
//! exposes velocity but not omega), so the spin is set once at the first step
//! by `init_pure_rolling` below — entirely within this example.
//!
//! ```bash
//! cargo run --release --example bench_rolling_decay --no-default-features -- examples/bench_rolling_decay/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Records the rolling sphere's state each step and writes a CSV at the end.
struct RollingTracker {
    /// Tag of the single rolling sphere.
    mover_tag: Option<u32>,
    /// Spin set yet? (one-shot pure-rolling initialization)
    initialized: bool,
    /// Buffered rows: (t, x, vx, omega) for the mover.
    rows: Vec<(f64, f64, f64, f64)>,
    /// Particle radius, captured for the header.
    radius: f64,
    written: bool,
}

impl RollingTracker {
    fn new() -> Self {
        Self {
            mover_tag: None,
            initialized: false,
            rows: Vec::new(),
            radius: 0.0,
            written: false,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin) // body force F = m g
        .add_plugins(WallPlugin); // frictional + rolling floor wall (z = 0)

    app.add_resource(RollingTracker::new());

    // Set the initial spin (pure rolling) before forces are computed.
    app.add_update_system(init_pure_rolling, ParticleSimScheduleSet::PreForce);
    // Record state after the integration completes.
    app.add_update_system(record, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// One-shot setup, run before the first force pass. Seats the single sphere on
/// the floor wall (z = 0, normal +z) with a hair of overlap, and gives it the
/// matching spin ω = (0, +v0/R, 0) so it starts in pure rolling (v = ωR,
/// contact point at rest).
fn init_pure_rolling(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    mut tracker: ResMut<RollingTracker>,
) {
    if tracker.initialized || atoms.nlocal < 1 {
        return;
    }
    let mut dem = registry.expect_mut::<DemAtom>("init_pure_rolling");

    let i = 0usize;
    let r = dem.radius[i];

    // Desired horizontal launch speed: take it from the inserted velocity.
    let v0 = atoms.vel[i][0];

    // Seat the sphere on the floor wall (z = 0) with a small overlap so normal
    // contact is live and near force balance (≈ the Hertz static overlap that
    // supports the weight) from step 0.
    let overlap0 = 2.0e-6;
    let cx = atoms.pos[i][0];
    let cy = atoms.pos[i][1];
    atoms.pos[i] = [cx, cy, r - overlap0];
    atoms.vel[i] = [v0, 0.0, 0.0];

    // Pure rolling in +x on a floor BELOW the sphere (contact point at the
    // bottom, lever r·n̂ = −R ẑ) ⇒ ω = (0, +v0/R, 0) makes the contact-point
    // velocity vanish: v + ω×(−Rẑ) = v − ω_y R x̂ = 0.
    dem.omega[i] = [0.0, v0 / r, 0.0];

    tracker.mover_tag = Some(atoms.tag[i]);
    tracker.radius = r;
    tracker.initialized = true;
}

/// Record (t, x, vx, omega_y) for the rolling sphere each step, and flush to CSV
/// once the sphere has effectively stopped (or at the last step).
fn record(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    input: Res<Input>,
    mut tracker: ResMut<RollingTracker>,
) {
    if tracker.written || atoms.nlocal == 0 {
        return;
    }
    let tag = match tracker.mover_tag {
        Some(t) => t,
        None => return,
    };
    let dem = registry.expect::<DemAtom>("record");

    let i = match (0..atoms.nlocal as usize).find(|&k| atoms.tag[k] == tag) {
        Some(k) => k,
        None => return,
    };

    let step = run_state.total_cycle;
    let dt = atoms.dt;
    let t = step as f64 * dt;
    let x = atoms.pos[i][0];
    let vx = atoms.vel[i][0];
    // Forward-rolling spin is ω_y = +vx/R; report ω_y directly so that, in pure
    // rolling, omega·R = vx and the two decay together.
    let omega = dem.omega[i][1];

    tracker.rows.push((t, x, vx, omega));

    // Stop logging once it has essentially halted, or at the final step.
    let last_step = run_config.current_stage(0).steps as usize;
    // Only treat near-zero speed as "halted" once the sphere has actually been
    // rolling (guards against the first step before the spin/velocity settle).
    let moved = tracker.rows.iter().any(|r| r.2 > 1e-3);
    let halted = moved && vx <= 1e-4;
    if halted || step + 1 >= last_step {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_rolling_decay".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/rolling_decay_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(f, "# radius={:.10e} dt={:.10e}", tracker.radius, dt).unwrap();
        writeln!(f, "t,x,vx,omega").unwrap();
        for (t, x, vx, w) in &tracker.rows {
            writeln!(f, "{:.10e},{:.10e},{:.10e},{:.10e}", t, x, vx, w).unwrap();
        }
        tracker.written = true;
        println!("=== Rolling Decay Results ===");
        println!("  rows recorded:    {}", tracker.rows.len());
        println!("  final vx:         {:.6e} m/s", vx);
        println!("  results saved to: {}", results_file);
    }
}
