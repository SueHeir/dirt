//! bench_marshall_twisting — validates the **Marshall** twisting model, whose
//! twisting stiffness/damping/friction coefficients are DERIVED from the active
//! tangential (Mindlin) model, against the exact analytical spin-down of a
//! sphere twisting on an enduring contact.
//!
//! ## Setup (identical seating to bench_twisting_friction)
//! Two identical spheres are stacked along +z. The LOWER sphere is a frozen
//! anchor (`[[freeze]]`, so it cannot translate OR spin — a true immovable
//! contact partner). The UPPER sphere rests on it under gravity, seated at the
//! static Hertz overlap so the normal contact force is `F_n = m g` with no
//! settling transient, and is given a pure twisting spin `ω = (0, 0, ω0)` about
//! the contact normal n̂ = +ẑ.
//!
//! A pure twist about n̂ produces **zero** relative surface velocity at the
//! contact point (it lies on the spin axis), so there is no sliding and no
//! rolling — the ONLY contact torque is the twisting couple, applied purely
//! along n̂. Crucially, even though the tangential (sliding) friction μ_t is now
//! nonzero (Marshall needs it — it sets μ_twist), the tangential spring stays at
//! zero because the contact-point slip is zero, so no tangential force or torque
//! is produced; the spin decays purely under the Marshall twisting couple.
//!
//! ## Analytical spin-down (Marshall derived coefficients)
//! Per LAMMPS `pair_granular` `twisting marshall` (Marshall 2009, eqs 32–33) the
//! twisting friction coefficient is `μ_twist = (2/3) a μ_t`, with `a = √(R* δ)`
//! the Hertz contact radius and μ_t the tangential friction. The saturated
//! (sliding) twisting couple is therefore
//!
//!     τ_tw = μ_twist · F_n = (2/3) a μ_t · F_n.
//!
//! With `F_n = m g`, equal-sphere `r_eff = R/2` (so `a = √(R δ / 2)`), and
//! `I = (2/5) m R²`, the spin decays at the exact constant rate
//!
//!     α = dω/dt = − τ_tw / I = − (5/3) · a μ_t g / R²      (exact, saturated).
//!
//! The Marshall model winds a derived torsional spring `k_twist = ½ k_t a²` up to
//! the cap `τ_max = μ_twist F_n` within a few steps and then holds it (the twist
//! keeps growing, so the slider stays saturated), giving the constant
//! deceleration above; `sweep.py` fits the slope past the brief wind-up.
//!
//! The recorder logs `t, omega_z, omega_perp, drift` so `sweep.py` can fit the
//! spin-down slope, check it against α above, and confirm the twist stayed pure
//! (ω_perp ≈ 0, drift ≈ 0).
//!
//! The upper sphere has no initial-spin TOML knob (`[[particles.insert]]`
//! exposes velocity but not omega), so the spin and the exact seating are set
//! once at the first step by `init_twist` below — entirely within this example.
//!
//! ```bash
//! cargo run --release --example bench_marshall_twisting --no-default-features --features precision-double -- examples/bench_marshall_twisting/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Initial twisting spin ω0 (rad/s) given to the upper sphere about +ẑ. There is
/// no omega knob in `[[particles.insert]]`, so it is applied here. Chosen so the
/// slowest-μ_t case still spins down within the configured step budget.
const OMEGA0: f64 = 8.0;

/// Static Hertz overlap (m) that makes the normal force balance the upper
/// sphere's weight, so `F_n = m g` from step 0 with no settling ring-down.
///
/// From `F_n = (4/3) E* √(r_eff) · δ^{3/2} = m g` with the config's fixed
/// R = 5 mm, ρ = 2500 kg/m³, E = 1e8 Pa, ν = 0.3 (⇒ E* = E/(2(1−ν²)),
/// r_eff = R/2): δ ≈ 2.307e-6 m. These material/geometry values are held fixed
/// across the sweep (only μ_t changes), so this seating stays exact for every
/// case — and, because the Marshall contact radius `a = √(r_eff δ)` depends on
/// this overlap, the derived twisting coefficients are likewise exact.
const SEAT_OVERLAP: f64 = 2.307e-6;

/// Records the twisting sphere's state each step and writes a CSV at the end.
struct TwistTracker {
    /// Tag of the spinning (upper) sphere.
    mover_tag: Option<u32>,
    /// Seated + spun yet? (one-shot initialization)
    initialized: bool,
    /// Buffered rows: (t, omega_z, omega_perp, drift) for the mover.
    rows: Vec<(f64, f64, f64, f64)>,
    /// Particle radius, captured for the header.
    radius: f64,
    /// Horizontal start position (x, y) of the mover, for the drift check.
    start_xy: (f64, f64),
    written: bool,
}

impl TwistTracker {
    fn new() -> Self {
        Self {
            mover_tag: None,
            initialized: false,
            rows: Vec::new(),
            radius: 0.0,
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

    app.add_resource(TwistTracker::new());

    // Seat + spin the upper sphere before forces are computed.
    app.add_update_system(init_twist, ParticleSimScheduleSet::PreForce);
    // Record state after the integration completes.
    app.add_update_system(record, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// One-shot setup, run before the first force pass. Picks the UPPER sphere as
/// the spinning mover (the lower one is the frozen anchor), seats it on the
/// anchor at the static Hertz overlap (co-axial in x,y so the contact normal is
/// exactly +ẑ), and gives it a pure twisting spin ω = (0, 0, ω0).
fn init_twist(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    mut tracker: ResMut<TwistTracker>,
) {
    if tracker.initialized || atoms.nlocal < 2 {
        return;
    }
    let mut dem = registry.expect_mut::<DemAtom>("init_twist");

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

    // Pure twisting spin about the contact normal (+ẑ): no sliding, no rolling.
    dem.omega[mover] = [0.0, 0.0, OMEGA0];

    tracker.mover_tag = Some(atoms.tag[mover]);
    tracker.radius = r as f64;
    tracker.start_xy = (ax, ay);
    tracker.initialized = true;
}

/// Record (t, omega_z, omega_perp, drift) for the twisting sphere each step, and
/// flush to CSV once the spin has effectively stopped (or at the last step).
fn record(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    input: Res<Input>,
    mut tracker: ResMut<TwistTracker>,
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
    let ox = dem.omega[i][0];
    let oy = dem.omega[i][1];
    let oz = dem.omega[i][2];
    let omega_perp = (ox * ox + oy * oy).sqrt();
    let dx = atoms.pos[i][0] as f64 - tracker.start_xy.0;
    let dy = atoms.pos[i][1] as f64 - tracker.start_xy.1;
    let drift = (dx * dx + dy * dy).sqrt();

    tracker.rows.push((t, oz, omega_perp, drift));

    // Stop logging once the spin has essentially halted, or at the final step.
    let last_step = run_config.current_stage(0).steps as usize;
    let moved = tracker.rows.iter().any(|r| r.1 > 1e-2);
    let halted = moved && oz <= 1e-3;
    if halted || step + 1 >= last_step {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_marshall_twisting".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let results_file = format!("{}/twisting_results.csv", data_dir);
        let mut f = fs::File::create(&results_file)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
        writeln!(
            f,
            "# radius={:.10e} dt={:.10e} omega0={:.10e}",
            tracker.radius, dt, OMEGA0
        )
        .unwrap();
        writeln!(f, "t,omega_z,omega_perp,drift").unwrap();
        for (t, oz, op, dr) in &tracker.rows {
            writeln!(f, "{:.10e},{:.10e},{:.10e},{:.10e}", t, oz, op, dr).unwrap();
        }
        tracker.written = true;
        println!("=== Marshall Twisting Results ===");
        println!("  rows recorded:    {}", tracker.rows.len());
        println!("  final omega_z:    {:.6e} rad/s", oz);
        println!("  results saved to: {}", results_file);
    }
}
