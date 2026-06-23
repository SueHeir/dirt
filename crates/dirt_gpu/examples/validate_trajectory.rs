//! Tier 2: full GPU *trajectory* validation against the recorded CPU baseline.
//!
//! Runs baseline example scenarios END TO END on the resident GPU path and diffs
//! the measured physics metrics against the CPU-single baseline (the right
//! reference: the GPU is f32 = single). CPU-double is shown for context.
//!
//!   - hertz_rebound  : sphere on a wall, normal Hertz → COR, contact_time, overlap
//!   - sliding_friction: sphere sliding on a wall under gravity, Coulomb friction
//!     spins it up to rolling-without-slipping → final vx, omega_y
//!
//! Effective params come from dirt's own `MaterialTable` (same E/ν→e_eff and
//! restitution→β as the CPU), so the only difference is f32-vs-f64 arithmetic.
//!
//! ```text
//! cargo run --release -p dirt_gpu --example validate_trajectory \
//!     --no-default-features --features precision-double
//! ```

use dirt_atom::MaterialTable;
use dirt_gpu::{GpuContext, GranularConfig, GranularForce, WallForce};
use soil_gpu::{Boundary, GpuState, Grid, Plane};

const DENSITY: f64 = 2500.0;

fn mass_of(r: f64) -> f64 {
    DENSITY * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3)
}

/// Effective (e_eff, beta, g_eff, mu) from dirt's MaterialTable — identical to CPU.
fn matched(young: f64, poisson: f64, e: f64, fric: f64) -> (f32, f32, f32, f32) {
    let mut mt = MaterialTable::new();
    mt.add_material("m", young, poisson, e, fric, 0.0, 0.0);
    mt.build_pair_tables();
    (
        mt.e_eff_ij[0][0] as f32,
        mt.beta_ij[0][0] as f32,
        mt.g_eff_ij[0][0] as f32,
        mt.friction_ij[0][0] as f32,
    )
}

fn relrow(name: &str, gpu: f64, single: f64, double: f64) {
    let rs = if single != 0.0 { (gpu - single).abs() / single.abs() } else { f64::NAN };
    let rd = if double != 0.0 { (gpu - double).abs() / double.abs() } else { f64::NAN };
    println!("  {name:<12} {gpu:>15.8e} {single:>15.8e} {rs:>11.2e} {double:>15.8e} {rd:>11.2e}");
}

/// Floor-wall GpuState builder shared by both scenarios (single particle).
fn build(ctx: GpuContext, pos: [f32; 3], vel: [f32; 3], r: f32, dt: f32, gravity: [f32; 3],
         params: (f32, f32, f32, f32)) -> (GpuState, usize) {
    let (e_eff, beta, g_eff, mu) = params;
    let m = mass_of(r as f64) as f32;
    let posv = vec![pos];
    let velv = vec![vel];
    let radius = vec![r];
    let inv_mass = vec![1.0f32 / m];
    let inv_inertia = vec![1.0f32 / (0.4 * m * r * r)];
    let mut boundary = Boundary::new();
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])); // floor z=0
    let grid = Grid::from_positions(&posv, 2.0 * r);
    let grid2 = Grid::from_positions(&posv, 2.0 * r);
    let mut gpu = GpuState::new(ctx, 1, grid.total_cells);
    gpu.set_params(dt, gravity);
    gpu.set_state(&posv, &velv, &inv_mass, grid);
    let omega = gpu.add_aux_dof();
    gpu.set_aux_inv_coeff(omega, &inv_inertia);
    gpu.set_aux_state(omega, &[[0.0f32; 3]]);
    // Mirror pile.rs: register GranularForce too (single particle ⇒ no pairs, but
    // it initialises the same resident contact-history substrate WallForce uses).
    let cfg = GranularConfig { e_eff, beta, g_eff, mu, dt };
    gpu.add_force_hook(Box::new(GranularForce::new(&gpu, &grid2, omega, &radius, cfg)));
    gpu.add_force_hook(Box::new(WallForce::new(
        &gpu, omega, &radius, &boundary, e_eff, beta, g_eff, mu, dt,
    )));
    (gpu, omega)
}

fn hertz(ctx: GpuContext) {
    // glass E=70e9 ν=0.22 e=0.9 μ=0; dt from CPU baseline; no gravity.
    let r = 0.005f64;
    let dt = 7.6225511569e-7f32;
    let params = matched(70.0e9, 0.22, 0.9, 0.0);
    let (mut gpu, _om) = build(ctx, [0.0, 0.0, r as f32 + 0.0005], [0.0, 0.0, -1.0],
                               r as f32, dt, [0.0; 3], params);

    let (mut was, mut vimp, mut vreb, mut prev) = (false, 0.0f64, 0.0f64, -1.0f64);
    let (mut max_ov, mut s0, mut s1) = (0.0f64, 0usize, 0usize);
    for step in 0..20000 {
        gpu.run_steps(1);
        gpu.wait();
        let z = gpu.download_pos()[0][2] as f64;
        let vz = gpu.download_vel()[0][2] as f64;
        let ov = r - z;
        let cont = ov > 0.0;
        if !was && cont { was = true; vimp = prev; s0 = step; max_ov = ov; }
        else if was && cont { if ov > max_ov { max_ov = ov; } }
        else if was && !cont { vreb = vz; s1 = step; break; }
        if !cont { prev = vz; }
    }
    let ct = (s1 - s0) as f64 * dt as f64;
    let cor = (vreb / vimp).abs();
    println!("\n=== hertz_rebound (normal Hertz, wall) ===");
    println!("  metric            GPU(f32)        CPU-single   relΔ(s)        CPU-double   relΔ(d)");
    relrow("COR", cor, 9.0164524317e-1, 9.0166213568e-1);
    relrow("contact_t", ct, 3.5063735322e-5, 3.5063735322e-5);
    relrow("max_overlap", max_ov, 1.1249557137e-5, 1.1249383514e-5);
}

fn sliding(ctx: GpuContext) {
    // glass E=70e9 ν=0.22 e=0.3 μ=0.5; dt=2e-6; gravity on; vx0=1, zero spin.
    let r = 0.005f64;
    let dt = 2.0e-6f32;
    let params = matched(70.0e9, 0.22, 0.3, 0.5);
    // Start just above the floor (overlap ~0) so it settles gently then slides.
    let (mut gpu, omega) = build(ctx, [0.0, 0.0, r as f32 + 1.0e-6], [1.0, 0.0, 0.0],
                                 r as f32, dt, [0.0, 0.0, -9.81], params);
    println!("  [probe] chunk   z(rel r)      vz         vx        omega_y   contact");
    for c in 0..45 {
        gpu.run_steps(1000);
        if c % 5 == 0 || c == 44 {
            gpu.wait();
            let z = gpu.download_pos()[0][2] as f64;
            let v = gpu.download_vel()[0];
            let oy = gpu.download_aux_state(omega)[0][1] as f64;
            println!("  [probe] {:>5}  {:>+10.3e}  {:>+9.2e}  {:>8.5}  {:>8.3}  {}",
                     (c + 1) * 1000, z - r, v[2] as f64, v[0] as f64, oy, (z < r) as u8);
        }
    }
    gpu.wait();
    let vx = gpu.download_vel()[0][0] as f64;
    let omega_y = gpu.download_aux_state(omega)[0][1] as f64;
    println!("\n=== sliding_friction (tangential Coulomb, wall, gravity) ===");
    println!("  metric            GPU(f32)        CPU-single   relΔ(s)        CPU-double   relΔ(d)");
    relrow("vx_final", vx, 7.1499371529e-1, 7.1428571428e-1);
    relrow("omega_y_final", omega_y, 1.4299874306e2, 1.4285714286e2);
}

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("No GPU adapter available.");
        return;
    };
    println!("GPU adapter: {}", ctx.adapter_info);
    hertz(ctx.clone());
    sliding(ctx);
}
