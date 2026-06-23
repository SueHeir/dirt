//! Full resident GPU granular sim: gravity + planar walls (floor + box) +
//! particle friction. A loose block of a few thousand grains is dropped under
//! gravity into a square box and settles into a pile, entirely on-device.
//!
//! Run with: `cargo run -p dirt_gpu --example pile --release`
//!
//! Demonstrates the resident stepper via the Force-hook stack: soil's
//! `GpuState` (resident pos/vel/force + cell list + velocity-Verlet + rotation
//! aux-DOF) with dirt's `GranularForce` (Hertz-Mindlin contacts) and `WallForce`
//! (planar walls) registered as Force hooks. set_state/set_params once, run_steps
//! in chunks, occasional downloads. Prints sanity metrics over time:
//!   - total kinetic energy (rises on impact, then DECAYS toward ~0 as friction
//!     and damping dissipate it),
//!   - minimum particle z (stays >= ~-small: grains rest ON the floor, never
//!     tunnel through),
//!   - final pile height.

use dirt_gpu::{GranularConfig, GranularForce, WallForce};
use soil_gpu::{Boundary, GpuContext, GpuState, Grid, Plane};

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("No GPU adapter available; cannot run pile example.");
        return;
    };

    // ── Particle block ────────────────────────────────────────────────────
    let r = 0.05f32; // grain radius (m)
    let side = 14usize; // 14^3 = 2744 grains
    let spacing = 2.05 * r; // just touching -> tiny gaps
    let mut pos = Vec::new();
    for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                let f = (ix + iy * side + iz * side * side) as f64;
                pos.push([
                    1.5 * r + ix as f32 * spacing + (0.13 * f).sin() as f32 * 0.03 * r,
                    1.5 * r + iy as f32 * spacing + (0.27 * f).cos() as f32 * 0.03 * r,
                    1.5 * r + iz as f32 * spacing,
                ]);
            }
        }
    }
    let n = pos.len();
    let vel = vec![[0.0f32; 3]; n];
    let om = vec![[0.0f32; 3]; n];
    let radius = vec![r; n];

    let mass = 1.0e-3f32; // light grains -> gravity vs. stiffness well-conditioned
    let inv_mass = vec![1.0 / mass; n];
    let inv_inertia = vec![1.0 / (0.4 * mass * r * r); n]; // solid sphere

    // ── Material params ───────────────────────────────────────────────────
    let e_eff = 1.0e5f32; // effective Young's modulus
    let beta = 0.5f32; // normal/tangential damping coefficient
    let g_eff = 4.0e4f32; // effective shear modulus
    let mu = 0.5f32; // sliding friction (>0 -> real pile)
    let gravity = [0.0f32, 0.0, -9.81];

    // Stable dt from the Hertz contact period tc = 2*pi*sqrt(m/k_n), evaluated
    // at a modest overlap; use dt ~ tc/40 (conservative).
    let delta_est = 0.05 * r;
    let k_n = (4.0 / 3.0) * e_eff * (delta_est * r).sqrt();
    let tc = 2.0 * std::f32::consts::PI * (mass / k_n).sqrt();
    let dt = tc / 40.0;

    // ── Walls: floor + four sides (a square box) ──────────────────────────
    let box_w = side as f32 * spacing + r;
    let mut boundary = Boundary::new();
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));    // floor z=0
    boundary.push(Plane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));    // x=0
    boundary.push(Plane::new([box_w, 0.0, 0.0], [-1.0, 0.0, 0.0])); // x=box_w
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]));    // y=0
    boundary.push(Plane::new([0.0, box_w, 0.0], [0.0, -1.0, 0.0])); // y=box_w

    // soil's Grid takes the literal contact cutoff (sum of radii = 2*r).
    let grid = Grid::from_positions(&pos, 2.0 * r);

    println!("GPU adapter: {}", ctx.adapter_info);
    println!(
        "pile: n={n} grains  r={r}  box_w={box_w:.3}  k_n~{k_n:.0}  tc~{tc:.2e}s  dt={dt:.2e}s",
    );
    println!("grid: cells={:?} total={} binsize={:.3}", grid.n, grid.total_cells, grid.bin_size);

    // ── Resident setup via the Force-hook stack (uploaded once) ───────────
    let mut gpu = GpuState::new(ctx, n, grid.total_cells);
    gpu.set_params(dt, gravity);
    gpu.set_state(&pos, &vel, &inv_mass, grid);
    let omega = gpu.add_aux_dof(); // rotation: state=ω, rate=τ, inv_coeff=1/inertia
    gpu.set_aux_inv_coeff(omega, &inv_inertia);
    gpu.set_aux_state(omega, &om);
    let cfg = GranularConfig { e_eff, beta, g_eff, mu, dt };
    gpu.add_force_hook(Box::new(GranularForce::new(&gpu, &grid, omega, &radius, cfg)));
    gpu.add_force_hook(Box::new(WallForce::new(&gpu, omega, &radius, &boundary, e_eff, beta, g_eff, mu, dt)));

    let i_solid = (0.4 * mass * r * r) as f64;
    let kinetic = |v: &[[f32; 3]], w: &[[f32; 3]]| -> f64 {
        let mut e = 0.0;
        for i in 0..n {
            let vs = (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) as f64;
            let ws = (w[i][0] * w[i][0] + w[i][1] * w[i][1] + w[i][2] * w[i][2]) as f64;
            e += 0.5 * mass as f64 * vs + 0.5 * i_solid * ws;
        }
        e
    };

    // ── Run resident on-device, downloading only occasionally ─────────────
    let chunk = 500;
    let chunks = 80; // 40,000 steps total
    println!("\n  step      KE          min z     max z");
    let mut ke_peak = 0.0f64;
    let start = std::time::Instant::now();
    for c in 0..=chunks {
        if c > 0 {
            gpu.run_steps(chunk);
            gpu.wait();
        }
        let p = gpu.download_pos();
        let v = gpu.download_vel();
        let w = gpu.download_aux_state(omega);
        let ke = kinetic(&v, &w);
        ke_peak = ke_peak.max(ke);
        let mut min_z = f32::MAX;
        let mut max_z = f32::MIN;
        for pi in &p {
            min_z = min_z.min(pi[2]);
            max_z = max_z.max(pi[2]);
        }
        if c % 5 == 0 || c == chunks {
            println!(
                "  {:>5}  {:>10.3e}  {:>8.4}  {:>8.4}",
                c * chunk, ke, min_z, max_z
            );
        }
    }
    let elapsed = start.elapsed();

    let p = gpu.download_pos();
        let v = gpu.download_vel();
        let w = gpu.download_aux_state(omega);
    let ke_final = kinetic(&v, &w);
    let mut min_z = f32::MAX;
    let mut max_z = f32::MIN;
    for pi in &p {
        min_z = min_z.min(pi[2]);
        max_z = max_z.max(pi[2]);
    }
    let pile_height = max_z + r;

    println!("\n── settled ──");
    println!(
        "KE: peak={ke_peak:.3e}  final={ke_final:.3e}  (decayed to {:.2e} of peak)",
        ke_final / ke_peak.max(1e-30)
    );
    println!("floor: min particle z = {min_z:.4}  (>= ~0 -> grains rest on the floor, no tunneling)");
    println!("pile height = {pile_height:.4} m");
    println!(
        "{} steps of {n} grains in {:.2}s ({:.2} Msteps*grains/s)",
        (chunks * chunk),
        elapsed.as_secs_f64(),
        (chunks * chunk * n) as f64 / elapsed.as_secs_f64() / 1e6
    );
}
