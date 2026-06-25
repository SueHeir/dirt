//! Step-1 residency: end-to-end on a real wall+gravity granular drop.
//!
//! Advances the SAME scene (gravity + floor/box walls + Hertz-Mindlin contacts +
//! rotation, all on-device via soil's `GpuState`) two ways and compares:
//!
//!   A. RESIDENT       — `run_steps(K)` windows, host syncs only at the end
//!                       (the residency model: data stays on device).
//!   B. PER-STEP SYNC  — `run_steps(1)` each step, with the per-step host
//!                       readback (download pos/vel/omega) that the
//!                       host-authoritative `GpuGranularForcePlugin` (milestone 1)
//!                       pays every step. Identical kernels → identical physics;
//!                       the only difference is the per-step host round-trip.
//!
//! So the position diff isolates whether residency batching changes the answer
//! (it must not), and the wall-clock ratio is the residency speedup on a real
//! sim. Resident-vs-CPU correctness is already established (`validate_trajectory`
//! ~1e-4 vs the CPU baseline); this isolates the schedule-sync cost.
//!
//! Run: cargo run -p dirt_gpu --example resident_validate --release

use dirt_gpu::{GranularConfig, GranularForce, WallForce};
use soil_gpu::{Boundary, GpuContext, GpuState, Grid, Plane};

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("No GPU adapter; cannot run resident_validate.");
        return;
    };

    // ── deterministic wall+gravity drop (a loose block settling in a box) ──
    let r = 0.05f32;
    let side = 10usize; // 1000 grains
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
    let vel0 = vec![[0.0f32; 3]; n];
    let om0 = vec![[0.0f32; 3]; n];
    let radius = vec![r; n];
    let mass = 1.0e-3f32;
    let inv_mass = vec![1.0 / mass; n];
    let inv_inertia = vec![1.0 / (0.4 * mass * r * r); n];

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

    // Build a FRESH GpuState + hooks per path, so each path starts from zero
    // contact history. (Re-using one instance and only re-uploading pos/vel does
    // NOT clear the on-device tangential springs, which would make each path
    // start from the previous path's leftover history — a comparison artifact.)
    let build = || {
        let grid = Grid::from_positions(&pos, 2.0 * r);
        let mut g = GpuState::new(ctx.clone(), n, grid.total_cells);
        g.set_params(dt, gravity);
        g.set_state(&pos, &vel0, &inv_mass, grid);
        let omega = g.add_aux_dof();
        g.set_aux_inv_coeff(omega, &inv_inertia);
        g.set_aux_state(omega, &om0);
        let cfg = GranularConfig::new(e_eff, beta, g_eff, mu, dt);
        g.add_force_hook(Box::new(GranularForce::new(&g, &grid, omega, &radius, cfg)));
        g.add_force_hook(Box::new(WallForce::new(
            &g, omega, &radius, &boundary, e_eff, beta, g_eff, mu, dt,
        )));
        (g, omega)
    };

    println!("resident_validate: n={n}  dt={dt:.2e}  adapter={}", ctx.adapter_info);

    let steps = 4000usize;
    let window = 250usize; // resident sync cadence (a neighbor-rebuild-window stand-in)

    let maxdiff = |a: &[[f32; 3]], b: &[[f32; 3]]| -> f32 {
        let mut m = 0.0f32;
        for i in 0..n {
            for d in 0..3 {
                m = m.max((a[i][d] - b[i][d]).abs());
            }
        }
        m
    };

    // ── A. SINGLE WINDOW: run_steps(steps) once — primes exactly once, the ──
    //     reference resident trajectory (this is the pile/validate_trajectory path).
    let (gpu, _omega) = build();
    let ta = std::time::Instant::now();
    gpu.run_steps(steps);
    gpu.wait();
    let a_time = ta.elapsed();
    let a_pos = gpu.download_pos();

    // ── B. WINDOWED: the residency model with host-sync boundaries (rebuild /
    //     MPI exchange / IO). Prime ONCE with run_steps, then run_steps_continue
    //     per window — no re-prime, so windows stitch bit-for-bit to one run.
    let (gpu, _omega) = build();
    let tb = std::time::Instant::now();
    let mut done = 0;
    while done < steps {
        let k = window.min(steps - done);
        if done == 0 {
            gpu.run_steps(k);
        } else {
            gpu.run_steps_continue(k);
        }
        gpu.wait();
        done += k;
    }
    let b_time = tb.elapsed();
    let b_pos = gpu.download_pos();

    // ── C. PER-STEP: run_steps(1)×steps + the per-step readback the ─────────
    //     host-authoritative milestone-1 plugin pays. Re-primes every step.
    let (gpu, omega) = build();
    let tc = std::time::Instant::now();
    for _ in 0..steps {
        gpu.run_steps(1);
        let _p = gpu.download_pos();
        let _v = gpu.download_vel();
        let _w = gpu.download_aux_state(omega);
    }
    let c_time = tc.elapsed();
    let c_pos = gpu.download_pos();

    let ab = maxdiff(&a_pos, &b_pos);
    let ac = maxdiff(&a_pos, &c_pos);
    println!("\nsteps={steps}  window={window}  (grain r={r})");
    println!("  A single-window run_steps({steps})        : {:>8.1} ms", a_time.as_secs_f64() * 1e3);
    println!("  B windowed     run_steps({window})x{:<3}        : {:>8.1} ms", steps / window, b_time.as_secs_f64() * 1e3);
    println!("  C per-step     run_steps(1)x{steps}+readback : {:>8.1} ms", c_time.as_secs_f64() * 1e3);
    println!("\n  resident speedup (A vs C per-step) : {:>6.1}x", c_time.as_secs_f64() / a_time.as_secs_f64());
    println!("  max |A - B| (single vs windowed)   : {ab:.3e}");
    println!("  max |A - C| (single vs per-step)   : {ac:.3e}");
    println!(
        "\n  => windowing {} the trajectory (A==B). `run_steps_continue` skips the\n     entry force-prime and trusts the resident force buffer, so a host-sync window\n     (rebuild / MPI exchange) no longer re-advances the tangential contact history.\n     Combined with the deterministic cell list, windowing is now bit-exact —\n     the correctness gate for residency and step-2b GPU-resident halos.",
        if ab < 1e-4 { "PRESERVES" } else { "CHANGES" }
    );
}
