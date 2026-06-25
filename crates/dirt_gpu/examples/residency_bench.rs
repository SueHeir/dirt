//! Step-1 residency benchmark: quantify the win from keeping the granular sim
//! RESIDENT on the GPU vs the host-authoritative per-step-sync model.
//!
//! Run: `cargo run -p dirt_gpu --example residency_bench --release`
//!
//! Two ways to advance the SAME granular scene K steps, both using soil's
//! `GpuState` + dirt's `GranularForce`/`WallForce` hooks + gravity:
//!
//!   RESIDENT       — `run_steps(K)`: one command encoder, one submit. Cell-list
//!                    rebuild + force + velocity-Verlet integrate all on-device
//!                    for all K steps; host touches nothing until a final
//!                    download. This is the `pile`/`validate_trajectory` path.
//!
//!   PER-STEP SYNC  — the host-authoritative schedule plugin's cost model
//!                    (`dirt_granular::gpu_granular_force`): every step re-upload
//!                    pos/vel (+ cell-list rebuild) via `set_state`, `eval_force_once`,
//!                    download force + torque, then integrate on the host.
//!
//! Correctness of the resident integrator itself is already established
//! (validate_trajectory: resident run_steps matches the CPU baseline to ~1e-4;
//! the host-authoritative plugin matches CPU hertz to 1.2e-8). This bench adds
//! (a) a resident-batched-vs-stepwise identity check and (b) the throughput
//! delta that motivates moving the schedule path onto the resident model.

use std::time::Instant;

use dirt_gpu::{GranularConfig, GranularForce, WallForce};
use soil_gpu::{Boundary, GpuContext, GpuState, Grid, Plane};

struct Scene {
    pos: Vec<[f32; 3]>,
    vel: Vec<[f32; 3]>,
    om: Vec<[f32; 3]>,
    radius: Vec<f32>,
    inv_mass: Vec<f32>,
    inv_inertia: Vec<f32>,
    cfg: GranularConfig,
    boundary: Boundary,
    grid: Grid,
    gravity: [f32; 3],
    dt: f32,
    n: usize,
}

fn build_scene(side: usize) -> Scene {
    let r = 0.05f32;
    let spacing = 2.05 * r;
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
    let mass = 1.0e-3f32;
    let e_eff = 1.0e5f32;
    let beta = 0.5f32;
    let g_eff = 4.0e4f32;
    let mu = 0.5f32;
    let gravity = [0.0f32, 0.0, -9.81];
    let delta_est = 0.05 * r;
    let k_n = (4.0 / 3.0) * e_eff * (delta_est * r).sqrt();
    let tc = 2.0 * std::f32::consts::PI * (mass / k_n).sqrt();
    let dt = tc / 40.0;

    let box_w = side as f32 * spacing + r;
    let mut boundary = Boundary::new();
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
    boundary.push(Plane::new([box_w, 0.0, 0.0], [-1.0, 0.0, 0.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]));
    boundary.push(Plane::new([0.0, box_w, 0.0], [0.0, -1.0, 0.0]));

    let grid = Grid::from_positions(&pos, 2.0 * r);
    Scene {
        pos,
        vel: vec![[0.0; 3]; n],
        om: vec![[0.0; 3]; n],
        radius: vec![r; n],
        inv_mass: vec![1.0 / mass; n],
        inv_inertia: vec![1.0 / (0.4 * mass * r * r); n],
        cfg: GranularConfig::new(e_eff, beta, g_eff, mu, dt),
        boundary,
        grid,
        gravity,
        dt,
        n,
    }
}

/// Fresh resident GpuState with granular + wall hooks and (optionally) gravity.
fn make_state(ctx: &GpuContext, s: &Scene, gravity_on: bool) -> (GpuState, usize) {
    let mut gpu = GpuState::new(ctx.clone(), s.n, s.grid.total_cells);
    gpu.set_params(s.dt, if gravity_on { s.gravity } else { [0.0; 3] });
    gpu.set_state(&s.pos, &s.vel, &s.inv_mass, s.grid);
    let omega = gpu.add_aux_dof();
    gpu.set_aux_inv_coeff(omega, &s.inv_inertia);
    gpu.set_aux_state(omega, &s.om);
    gpu.add_force_hook(Box::new(GranularForce::new(&gpu, &s.grid, omega, &s.radius, s.cfg)));
    gpu.add_force_hook(Box::new(WallForce::new(
        &gpu, omega, &s.radius, &s.boundary, s.cfg.e_eff, s.cfg.beta, s.cfg.g_eff, s.cfg.mu, s.dt,
    )));
    (gpu, omega)
}

fn max_pos_diff(a: &[[f32; 3]], b: &[[f32; 3]]) -> f32 {
    a.iter()
        .zip(b)
        .map(|(p, q)| {
            (0..3)
                .map(|d| (p[d] - q[d]).abs())
                .fold(0.0f32, f32::max)
        })
        .fold(0.0f32, f32::max)
}

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("No GPU adapter; cannot run residency_bench.");
        return;
    };
    println!("GPU adapter: {}", ctx.adapter_info);
    let k = 500usize; // steps per neighbour-rebuild window (resident chunk)

    println!("\n  N        K     resident(ms)  per-step-sync(ms)  speedup   batch-vs-stepwise");
    for &side in &[20usize, 40] {
        let s = build_scene(side);
        let n = s.n;

        // ── Correctness: run_steps(K) [1 submit] vs K×run_steps(1) [K submits] ──
        let (a, _) = make_state(&ctx, &s, true);
        a.run_steps(k);
        a.wait();
        let pa = a.download_pos();
        let (b, _) = make_state(&ctx, &s, true);
        for _ in 0..k {
            b.run_steps(1);
        }
        b.wait();
        let pb = b.download_pos();
        let ident = max_pos_diff(&pa, &pb);

        // ── Timing: resident run_steps(K) ──
        let (res, _) = make_state(&ctx, &s, true);
        let t0 = Instant::now();
        res.run_steps(k);
        res.wait();
        let _ = res.download_pos();
        let t_res = t0.elapsed().as_secs_f64() * 1e3;

        // ── Timing: per-step-sync (host-authoritative plugin cost model) ──
        let (hs, omega) = make_state(&ctx, &s, true);
        let mut pos = s.pos.clone();
        let mut vel = s.vel.clone();
        let mut om = s.om.clone();
        let dt = s.dt;
        let t0 = Instant::now();
        for _ in 0..k {
            // re-upload moved state + cell-list rebuild (set_state), force-only eval, download
            hs.set_state(&pos, &vel, &s.inv_mass, s.grid);
            hs.set_aux_state(omega, &om);
            hs.eval_force_once();
            let f = hs.download_force();
            let tq = hs.download_aux_rate(omega);
            // host integrate (semi-implicit Euler — representative per-step host work;
            // the exact integrator is irrelevant to the sync-cost measurement)
            for i in 0..n {
                for d in 0..3 {
                    vel[i][d] += dt * f[i][d] * s.inv_mass[i];
                    om[i][d] += dt * tq[i][d] * s.inv_inertia[i];
                    pos[i][d] += dt * vel[i][d];
                }
            }
        }
        let t_hs = t0.elapsed().as_secs_f64() * 1e3;

        println!(
            "  {:<7} {:<5} {:>11.2}  {:>16.2}  {:>6.1}x   {:>.2e}",
            n,
            k,
            t_res,
            t_hs,
            t_hs / t_res,
            ident
        );
    }
    println!(
        "\nident = max|pos| between batched run_steps(K) and stepwise K×run_steps(1)\n\
         (≈0 ⇒ residency batching is exact). speedup = per-step-sync / resident."
    );
}
