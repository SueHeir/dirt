//! bench_sds_rolling — validates the SDS (spring-dashpot-slider) rolling-
//! resistance model (`rolling_model = "sds"`, contact.rs) against ITS OWN exact
//! analytical behaviour in two regimes:
//!
//!   1. **Elastic spring-dashpot decay** (cap disengaged): the rolling spin
//!      relaxes as the exact damped linear oscillator `I δ̈ + γ_r δ̇ + k_r δ = 0`.
//!      In the over-damped regime this is a clean (near-)exponential ω decay set
//!      by the larger eigenvalue |s₁| of that ODE.
//!   2. **Coulomb cap** (slider saturated): under a large sustained spin the
//!      spring+dashpot torque exceeds the cap `τ_max = μ_r·F_n·r_eff`, so the
//!      slider holds `τ = τ_max` and the spin decays at the exact constant rate
//!      `α = (5/4)·μ_r·g/R` (same closed form as the twisting/constant benches:
//!      F_n = m g, equal-sphere r_eff = R/2, I = (2/5) m R²).
//!
//! ## Setup (identical geometry to bench_twisting_friction, different spin axis)
//! Two identical spheres are stacked along +z. The LOWER sphere is a frozen
//! anchor (`[[freeze]]`, so it can neither translate nor spin). The UPPER sphere
//! rests on it under gravity, seated at the static Hertz overlap so `F_n = m g`
//! from step 0 with no settling ring-down, and is given a pure ROLLING spin
//! `ω = (ω0, 0, 0)` about a horizontal axis (⊥ to the contact normal n̂ = +ẑ).
//!
//! `ω` about x̂ is a pure rolling relative rotation (its projection onto n̂ is
//! zero → no twisting), so the rolling-resistance couple is the model under test.
//! Sliding friction is turned OFF in the config (`friction = 0`) so the tangential
//! surface slip that a horizontal spin would otherwise drive cannot add a couple:
//! the ONLY torque perpendicular to n̂ is the SDS rolling resistance. The rolling
//! couple is a pure torque (no force), so the sphere never translates — ω_x simply
//! spins down while ω_y, ω_z and the horizontal position stay at zero.
//!
//! The upper sphere has no initial-spin TOML knob (`[[particles.insert]]` exposes
//! velocity but not omega), so the exact seating and the spin are applied once at
//! the first step by `init_roll` below — entirely within this example. The spin
//! magnitude ω0 is read from `SDS_OMEGA0` (default `OMEGA0`) so sweep.py can drive
//! the same binary at different spins without editing the source.
//!
//! The recorder logs `t, omega_roll, omega_perp, drift` so sweep.py can fit the
//! elastic decay rate / the saturated slope and confirm the roll stayed pure
//! (ω_perp ≈ 0, drift ≈ 0).
//!
//! ```bash
//! cargo run --release --example bench_sds_rolling --no-default-features --features precision-double -- examples/bench_sds_rolling/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Default initial rolling spin ω0 (rad/s) about +x̂. Overridable via the
/// `SDS_OMEGA0` env var so one binary serves both the small-spin elastic case and
/// the large-spin saturated case without a source edit.
const OMEGA0: f64 = 8.0;

/// Static Hertz overlap (m) that makes the normal force balance the upper
/// sphere's weight, so `F_n = m g` from step 0 with no settling ring-down.
///
/// From `F_n = (4/3) E* √(r_eff) · δ^{3/2} = m g` with R = 5 mm, ρ = 2500 kg/m³,
/// E = 1e8 Pa, ν = 0.3 (⇒ E* = E/(2(1−ν²)), r_eff = R/2): δ ≈ 2.307e-6 m. These
/// material/geometry values are held fixed across the sweep (only the rolling
/// model parameters and ω0 change), so this seating stays exact for every case.
const SEAT_OVERLAP: f64 = 2.307e-6;

/// Reads the requested initial spin: `SDS_OMEGA0` if set & parseable, else OMEGA0.
fn omega0() -> f64 {
    std::env::var("SDS_OMEGA0")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|v| v.is_finite() && *v > 0.0)
        .unwrap_or(OMEGA0)
}

/// Records the rolling sphere's state each step and writes a CSV at the end.
struct RollTracker {
    /// Tag of the spinning (upper) sphere.
    mover_tag: Option<u32>,
    /// Seated + spun yet? (one-shot initialization)
    initialized: bool,
    /// Buffered rows: (t, omega_roll, omega_perp, drift) for the mover.
    rows: Vec<(f64, f64, f64, f64)>,
    /// Particle radius, captured for the header.
    radius: f64,
    /// Applied initial spin, captured for the header.
    omega0: f64,
    /// Horizontal start position (x, y) of the mover, for the drift check.
    start_xy: (f64, f64),
    written: bool,
}

impl RollTracker {
    fn new() -> Self {
        Self {
            mover_tag: None,
            initialized: false,
            rows: Vec::new(),
            radius: 0.0,
            omega0: omega0(),
            start_xy: (0.0, 0.0),
            written: false,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin) // body force F = m g on the free (upper) sphere
        .add_plugins(FixesPlugin); // [[freeze]] anchors the lower sphere

    app.add_resource(RollTracker::new());

    // Seat + spin the upper sphere before forces are computed.
    app.add_update_system(init_roll, ParticleSimScheduleSet::PreForce);
    // Record state after the integration completes.
    app.add_update_system(record, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// One-shot setup, run before the first force pass. Picks the UPPER sphere as the
/// spinning mover (the lower one is the frozen anchor), seats it on the anchor at
/// the static Hertz overlap (co-axial in x,y so the contact normal is exactly +ẑ),
/// and gives it a pure rolling spin ω = (ω0, 0, 0) about x̂.
fn init_roll(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    mut tracker: ResMut<RollTracker>,
) {
    if tracker.initialized || atoms.nlocal < 2 {
        return;
    }
    let mut dem = registry.expect_mut::<DemAtom>("init_roll");

    // Identify the two spheres by height: mover = upper, anchor = lower.
    let (mut mover, mut anchor) = (0usize, 1usize);
    if atoms.pos[0][2] < atoms.pos[1][2] {
        mover = 1;
        anchor = 0;
    }

    let r = dem.radius[mover];
    let ax = atoms.pos[anchor][0] as f64;
    let ay = atoms.pos[anchor][1] as f64;
    let az = atoms.pos[anchor][2] as f64;

    // Seat the mover co-axially above the anchor with the static overlap so the
    // contact normal is +ẑ and F_n = m g from step 0.
    let mz = az + 2.0 * (r as f64) - SEAT_OVERLAP;
    atoms.pos[mover] = [
        ax as dirt_core::soil_core::Real,
        ay as dirt_core::soil_core::Real,
        mz as dirt_core::soil_core::Real,
    ];
    atoms.vel[mover] = [0.0, 0.0, 0.0];

    // Pure rolling spin about x̂ (⊥ to n̂ = +ẑ): projection onto n̂ is zero, so this
    // is rolling, not twisting.
    dem.omega[mover] = [tracker.omega0, 0.0, 0.0];

    tracker.mover_tag = Some(atoms.tag[mover]);
    tracker.radius = r as f64;
    tracker.start_xy = (ax, ay);
    tracker.initialized = true;
}

/// Record (t, omega_roll, omega_perp, drift) for the rolling sphere each step, and
/// flush to CSV once the spin has effectively stopped (or at the last step).
fn record(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    input: Res<Input>,
    mut tracker: ResMut<RollTracker>,
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
    let ox = dem.omega[i][0]; // rolling spin (about x̂)
    let oy = dem.omega[i][1];
    let oz = dem.omega[i][2];
    let omega_perp = (oy * oy + oz * oz).sqrt(); // off-axis spin leak (should stay ~0)
    let dx = atoms.pos[i][0] as f64 - tracker.start_xy.0;
    let dy = atoms.pos[i][1] as f64 - tracker.start_xy.1;
    let drift = (dx * dx + dy * dy).sqrt();

    tracker.rows.push((t, ox, omega_perp, drift));

    // Stop logging once the spin has essentially halted, or at the final step.
    let last_step = run_config.current_stage(0).steps as usize;
    let moved = tracker
        .rows
        .iter()
        .any(|r| r.1.abs() > 1e-2 * tracker.omega0);
    let halted = moved && ox.abs() <= 1e-3 * tracker.omega0;
    if halted || step + 1 >= last_step {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_sds_rolling".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/sds_rolling_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "# radius={:.10e} dt={:.10e} omega0={:.10e}",
            tracker.radius, dt, tracker.omega0
        )
        .unwrap();
        writeln!(f, "t,omega_roll,omega_perp,drift").unwrap();
        for (t, ox, op, dr) in &tracker.rows {
            writeln!(f, "{:.10e},{:.10e},{:.10e},{:.10e}", t, ox, op, dr).unwrap();
        }
        tracker.written = true;
        println!("=== SDS Rolling Results ===");
        println!("  rows recorded:    {}", tracker.rows.len());
        println!("  omega0:           {:.6e} rad/s", tracker.omega0);
        println!("  final omega_roll: {:.6e} rad/s", ox);
        println!("  results saved to: {}", results_file);
    }
}
