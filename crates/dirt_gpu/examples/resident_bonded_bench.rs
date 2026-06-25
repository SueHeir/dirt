//! resident_bonded_bench — GPU resident throughput for a bonded periodic system
//! (the BPM-under-LEBC workload). A cubic lattice of grains in a periodic box,
//! each bonded to its +x/+y/+z neighbours (bonds wrap across the boundary, so the
//! periodic bond minimum-image is exercised), advanced fully resident on the GPU
//! with the contact hook + the beam bond hook composing. Reports particle-steps/s.
//!
//! Run: PERF_SIDE=12 PERF_STEPS=2000 cargo run --release -p dirt_gpu \
//!        --example resident_bonded_bench --no-default-features --features precision-double
//!
//! Compare the throughput to the CPU baseline (bench_lebc_shear: ~8 M part-steps/s
//! at 1634 grains, contact-only — bonds add CPU cost, so this is conservative).

use dirt_gpu::{
    BeamBondConfig, BeamBondForce, BondTopology, GpuContext, GpuState, GranularConfig,
    GranularForce, Grid,
};

fn main() {
    let side: usize = std::env::var("PERF_SIDE").ok().and_then(|s| s.parse().ok()).unwrap_or(12);
    let steps: usize = std::env::var("PERF_STEPS").ok().and_then(|s| s.parse().ok()).unwrap_or(2000);
    let n = side * side * side;

    let Some(ctx) = GpuContext::new() else {
        eprintln!("no GPU adapter; skipping");
        return;
    };

    let r = 0.1f32;
    let spacing = 2.05 * r; // grains nearly touching
    let l = side as f32 * spacing; // periodic box length
    let mass = (4.0 / 3.0) * std::f32::consts::PI * r * r * r;
    let inertia = 0.4 * mass * r * r;

    // Cubic lattice positions in [0, l)^3, with a tiny deterministic jitter.
    let mut pos = Vec::with_capacity(n);
    let idx = |ix: usize, iy: usize, iz: usize| (ix * side + iy) * side + iz;
    for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                let f = idx(ix, iy, iz) as f32;
                pos.push([
                    (ix as f32 + 0.5) * spacing + 0.01 * r * (0.13 * f).sin(),
                    (iy as f32 + 0.5) * spacing + 0.01 * r * (0.27 * f).cos(),
                    (iz as f32 + 0.5) * spacing,
                ]);
            }
        }
    }

    // Bond each grain to its +x/+y/+z neighbour (periodic), double-stored on both
    // endpoints → 6 bonds/grain. CSR topology.
    let mut offsets = vec![0u32; n + 1];
    let mut counts = vec![0u32; n];
    let neighbours = |ix: usize, iy: usize, iz: usize| {
        [
            idx((ix + 1) % side, iy, iz),
            idx(ix, (iy + 1) % side, iz),
            idx(ix, iy, (iz + 1) % side),
            idx((ix + side - 1) % side, iy, iz),
            idx(ix, (iy + side - 1) % side, iz),
            idx(ix, iy, (iz + side - 1) % side),
        ]
    };
    for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                counts[idx(ix, iy, iz)] = 6;
            }
        }
    }
    for i in 0..n {
        offsets[i + 1] = offsets[i] + counts[i];
    }
    let total = offsets[n] as usize;
    let mut partner = vec![0u32; total];
    let mut r0 = vec![spacing; total];
    for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                let i = idx(ix, iy, iz);
                let base = offsets[i] as usize;
                for (k, p) in neighbours(ix, iy, iz).iter().enumerate() {
                    partner[base + k] = *p as u32;
                }
            }
        }
    }
    let topo = BondTopology { offsets, partner, r0: std::mem::take(&mut r0) };

    // Periodic box grid (cells tile [0,l) exactly).
    let nc = ((l / (2.0 * r)).floor() as i32).max(3);
    let grid = Grid { n: [nc, nc, nc], origin: [0.0; 3], bin_size: l / nc as f32, total_cells: (nc * nc * nc) as usize };

    let dt = 1.0e-6f32;
    let mut gs = GpuState::new(ctx.clone(), n, grid.total_cells);
    gs.set_params(dt, [0.0; 3]);
    gs.set_box([l, l, l], [0.0; 3], 0.0, 0.0); // periodic (orthogonal; LE tilt=0 for the throughput measure)
    let vel = vec![[0.0f32; 3]; n];
    let inv_mass = vec![1.0 / mass; n];
    gs.set_state(&pos, &vel, &inv_mass, grid.clone());
    let omega = gs.add_aux_dof();
    gs.set_aux_inv_coeff(omega, &vec![1.0 / inertia; n]);
    gs.set_aux_state(omega, &vec![[0.0f32; 3]; n]);

    let radius = vec![r; n];
    let mut cc = GranularConfig::new(1.0e6, 0.2, 4.0e5, 0.3, dt);
    cc.lx = l; cc.ly = l; cc.lz = l;
    gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega, &radius, cc)));
    let bcfg = BeamBondConfig {
        bond_radius_ratio: 0.5, youngs_modulus: 1.0e7, shear_modulus: 4.0e6,
        beta_normal: 0.1, beta_shear: 0.1, beta_twist: 0.1, beta_bending: 0.1,
        sigma_max: 1.0e30, tau_max: 1.0e30, dt,
        lx: l, ly: l, lz: l, tilt_xy: 0.0, accumulate_torque: true,
    };
    gs.add_force_hook(Box::new(BeamBondForce::new(&gs, omega, &radius, &topo, bcfg)));

    // Warm up (build pipelines, prime), then time the resident window.
    gs.run_steps(50);
    gs.wait();
    let t0 = std::time::Instant::now();
    gs.run_steps(steps);
    gs.wait();
    let secs = t0.elapsed().as_secs_f64();

    let part_steps = (n as f64) * (steps as f64);
    let throughput = part_steps / secs;
    println!(
        "resident_bonded_bench: n={n} ({side}^3) bonds={total} steps={steps} box={l:.4} adapter={}",
        ctx.adapter_info
    );
    println!(
        "  wall={secs:.3}s  throughput={:.2} M part-steps/s",
        throughput / 1e6
    );
    println!("  (CPU bench_lebc_shear baseline ≈ 8 M part-steps/s at 1634 grains, contact-only)");
}
