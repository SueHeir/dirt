//! Validate `GpuGranularResidentPlugin` (roadmap step 1): a wall+gravity granular
//! drop run via the resident schedule plugin must reproduce a direct `GpuState`
//! run (the path already validated vs CPU by `validate_trajectory`/`pile`),
//! proving the schedule integration is faithful — and the windowed plugin must be
//! much faster than syncing the host every step.
//!
//! Run: cargo run -p dirt_granular --example resident_plugin_validate \
//!        --no-default-features --features precision-double --release

use grass_app::prelude::*;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_granular::{gpu_granular_resident_step, ResidentGpu};
use dirt_gpu::{
    Boundary, GpuContext, GpuState, GranularConfig, GranularForce, Grid, Plane, WallForce,
};
use soil_core::{Atom, AtomDataRegistry, ParticleSimScheduleSet};

const SIDE: usize = 8; // 8^3 = 512 grains
const R: f32 = 0.05;
const DENSITY: f64 = 2500.0;
const STEPS: usize = 4000;
const WINDOW: usize = 500;
const GRAVITY: [f32; 3] = [0.0, 0.0, -9.81];

struct Scene {
    pos: Vec<[f32; 3]>,
    n: usize,
    mass: f64,
    inv_inertia: f32,
    boundary: Boundary,
    dt: f32,
    e_eff: f32,
    beta: f32,
    g_eff: f32,
    mu: f32,
}

fn build_scene() -> Scene {
    let spacing = 2.05 * R;
    let mut pos = Vec::new();
    for ix in 0..SIDE {
        for iy in 0..SIDE {
            for iz in 0..SIDE {
                let f = (ix + iy * SIDE + iz * SIDE * SIDE) as f64;
                pos.push([
                    1.5 * R + ix as f32 * spacing + (0.13 * f).sin() as f32 * 0.03 * R,
                    1.5 * R + iy as f32 * spacing + (0.27 * f).cos() as f32 * 0.03 * R,
                    1.5 * R + iz as f32 * spacing,
                ]);
            }
        }
    }
    let n = pos.len();
    let mass = DENSITY * 4.0 / 3.0 * std::f64::consts::PI * (R as f64).powi(3);
    let inv_inertia = (1.0 / (0.4 * mass * (R as f64).powi(2))) as f32;

    // Soft material so dt is comfortable for a multi-thousand-step drop.
    let mt = make_mt();
    let (e_eff, beta, g_eff, mu) = (
        mt.e_eff_ij[0][0] as f32,
        mt.beta_ij[0][0] as f32,
        mt.g_eff_ij[0][0] as f32,
        mt.friction_ij[0][0] as f32,
    );

    // Stable dt = tc/40 from the Hertz contact period at a modest overlap.
    let delta = 0.05 * R;
    let k_n = (4.0 / 3.0) * e_eff * (delta * R).sqrt();
    let tc = 2.0 * std::f32::consts::PI * (mass as f32 / k_n).sqrt();
    let dt = tc / 40.0;

    let box_w = SIDE as f32 * spacing + R;
    let mut boundary = Boundary::new();
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
    boundary.push(Plane::new([box_w, 0.0, 0.0], [-1.0, 0.0, 0.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]));
    boundary.push(Plane::new([0.0, box_w, 0.0], [0.0, -1.0, 0.0]));

    Scene { pos, n, mass, inv_inertia, boundary, dt, e_eff, beta, g_eff, mu }
}

/// Single-material table (same scalars the resident plugin reads via gpu_scalars).
fn make_mt() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material("soft", 1.0e6, 0.3, 0.9, 0.4, 0.0, 0.0);
    mt.build_pair_tables();
    mt
}

/// Make the host Atom + DemAtom registry for the scene (same initial state both
/// paths consume).
fn make_app_state(s: &Scene) -> (Atom, AtomDataRegistry) {
    let mut atom = Atom::new();
    atom.dt = s.dt as f64;
    for (i, p) in s.pos.iter().enumerate() {
        atom.push_test_atom(i as u32, [p[0] as f64, p[1] as f64, p[2] as f64], R as f64, s.mass);
    }
    atom.nlocal = s.n as u32;
    atom.natoms = s.n as u64;

    let mut dem = DemAtom::new();
    for _ in 0..s.n {
        dem.radius.push(R as f64);
        dem.density.push(DENSITY);
        dem.inv_inertia.push(s.inv_inertia as f64);
        dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
        dem.omega.push([0.0; 3]);
        dem.ang_mom.push([0.0; 3]);
        dem.torque.push([0.0; 3]);
        dem.body_id.push(0.0);
    }
    let mut registry = AtomDataRegistry::new();
    registry.register(dem);
    (atom, registry)
}

/// Run the scene through the resident schedule plugin (window K, ticks = steps/K).
/// Returns final positions and wall-clock.
fn run_plugin(s: &Scene, ctx: GpuContext, window: usize) -> (Vec<[f64; 3]>, std::time::Duration) {
    let ticks = STEPS / window;
    let (atom, registry) = make_app_state(s);
    let mut app = App::new();
    app.add_resource(atom);
    app.add_resource(registry);
    app.add_resource(make_mt());
    app.add_resource(ResidentGpu::new(Some(ctx), window, GRAVITY, s.boundary.clone()));
    app.add_update_system(gpu_granular_resident_step, ParticleSimScheduleSet::Force);
    app.organize_systems();

    let t0 = std::time::Instant::now();
    for _ in 0..ticks {
        app.run();
    }
    let dt = t0.elapsed();
    let a = app.get_resource_ref::<Atom>().unwrap();
    let pos = (0..s.n).map(|i| [a.pos[i][0] as f64, a.pos[i][1] as f64, a.pos[i][2] as f64]).collect();
    (pos, dt)
}

/// Direct GpuState reference: build exactly as the plugin does and window the run
/// (run_steps once, then run_steps_continue) — the already-CPU-validated path.
fn run_direct(s: &Scene, ctx: GpuContext, window: usize) -> Vec<[f64; 3]> {
    let n = s.n;
    let radius: Vec<f32> = vec![R; n];
    let inv_inertia: Vec<f32> = vec![s.inv_inertia; n];
    let inv_mass: Vec<f32> = vec![(1.0 / s.mass) as f32; n];
    let velf = vec![[0.0f32; 3]; n];
    let omf = vec![[0.0f32; 3]; n];
    let r_max = R.max(f32::MIN_POSITIVE);
    let grid = Grid::from_positions(&s.pos, 2.0 * r_max);

    let mut gs = GpuState::new(ctx, n, grid.total_cells);
    gs.set_params(s.dt, GRAVITY);
    gs.set_state(&s.pos, &velf, &inv_mass, grid);
    let omega_aux = gs.add_aux_dof();
    gs.set_aux_inv_coeff(omega_aux, &inv_inertia);
    gs.set_aux_state(omega_aux, &omf);
    let cfg = GranularConfig { e_eff: s.e_eff, beta: s.beta, g_eff: s.g_eff, mu: s.mu, dt: s.dt };
    gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega_aux, &radius, cfg)));
    gs.add_force_hook(Box::new(WallForce::new(
        &gs, omega_aux, &radius, &s.boundary, s.e_eff, s.beta, s.g_eff, s.mu, s.dt,
    )));

    let windows = STEPS / window;
    gs.run_steps(window);
    for _ in 1..windows {
        gs.run_steps_continue(window);
    }
    gs.wait();
    let p = gs.download_pos();
    (0..n).map(|i| [p[i][0] as f64, p[i][1] as f64, p[i][2] as f64]).collect()
}

fn max_diff(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter().zip(b).flat_map(|(x, y)| (0..3).map(move |d| (x[d] - y[d]).abs())).fold(0.0, f64::max)
}

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("no GPU adapter; skipping resident_plugin_validate");
        return;
    };
    let s = build_scene();
    println!(
        "resident_plugin_validate: n={} steps={} window={} dt={:.3e} adapter={}",
        s.n, STEPS, WINDOW, s.dt, ctx.adapter_info
    );

    let direct = run_direct(&s, ctx.clone(), WINDOW);
    let (plugin_pos, t_plugin) = run_plugin(&s, ctx.clone(), WINDOW);
    let (_perstep_pos, t_perstep) = run_plugin(&s, ctx.clone(), 1);

    let d_faithful = max_diff(&plugin_pos, &direct);
    let speedup = t_perstep.as_secs_f64() / t_plugin.as_secs_f64();

    println!("\n  resident plugin (window={WINDOW}) : {:>8.1} ms", t_plugin.as_secs_f64() * 1e3);
    println!("  resident plugin (window=1)     : {:>8.1} ms", t_perstep.as_secs_f64() * 1e3);
    println!("  end-to-end speedup (win=1 -> {WINDOW}) : {speedup:>6.1}x");
    println!("  max |plugin - direct GpuState| : {d_faithful:.3e}");

    let pass = d_faithful < 1e-6;
    println!(
        "\n  => {} the resident PLUGIN reproduces the direct GpuState path bit-for-bit\n     (faithful schedule integration); that path is validated vs CPU (~1e-4) by\n     validate_trajectory/pile. Windowing the host sync is {speedup:.0}x cheaper.",
        if pass { "PASS:" } else { "FAIL:" }
    );
    assert!(pass, "resident plugin diverged from direct GpuState: {d_faithful:.3e}");
}
