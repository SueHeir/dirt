//! Windowed-resident GPU granular stepper (roadmap step 1).
//!
//! [`GpuGranularResidentPlugin`] is the device-resident counterpart to the
//! host-authoritative [`GpuGranularForcePlugin`](crate::gpu::GpuGranularForcePlugin)
//! (milestone 1). Instead of a host↔device round-trip every step (upload pos/vel,
//! eval force, download force, host-integrate), it keeps positions/velocities/
//! forces/omega **resident on the device** and advances the whole velocity-Verlet
//! loop (integrate + Hertz-Mindlin contact + planar walls + gravity) on the GPU
//! for a *window* of `K` steps per schedule tick via soil's `run_steps` /
//! `run_steps_continue`. The host `Atom`/`DemAtom` are synced back each window for
//! I/O and diagnostics; the device is the source of truth between syncs.
//!
//! ## Scope (single-rank, step 1)
//! - The resident `GpuState` rebuilds its own cell list on-device each step, so it
//!   needs **no host neighbor list**. Walls are planar (passed as a [`Boundary`]),
//!   gravity is a constant vector. No MPI (the window boundary is the future
//!   halo-exchange point — step 2b).
//! - This plugin **replaces** the host Verlet integrate, host neighbor build, host
//!   contact force, and host rotation for the particles it owns — so add it on its
//!   own (do NOT also add `VelocityVerletPlugin`/`RotationalDynamicsPlugin`/host
//!   force). Milestone 1's host-authoritative path is untouched and stays the
//!   default for MPI / mixed-physics configs.
//! - Bit-exactness across window boundaries depends on soil's deterministic cell
//!   list + `run_steps_continue` (the force-prime/history gate, already landed).

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_gpu::{
    Boundary, GpuContext, GpuState, GranularConfig, GranularForce, Grid, WallForce,
};
use soil_core::{Atom, AtomDataRegistry, ParticleSimScheduleSet, Real};

/// Plain Hertz-Mindlin scalars the GPU kernel uses (single material).
fn gpu_scalars(mt: &MaterialTable) -> (f32, f32, f32, f32) {
    (
        mt.e_eff_ij[0][0] as f32,
        mt.beta_ij[0][0] as f32,
        mt.g_eff_ij[0][0] as f32,
        mt.friction_ij[0][0] as f32,
    )
}

/// Resident GPU state: built lazily on the first step from the host `Atom` set,
/// then owns the trajectory on-device. `primed` distinguishes the first window
/// (`run_steps`, which evaluates the seed force) from later windows
/// (`run_steps_continue`, which trusts the resident force buffer — the bit-exact
/// windowing path).
pub struct ResidentGpu {
    ctx: Option<GpuContext>,
    gpu: Option<GpuState>,
    omega_aux: usize,
    primed: bool,
    window: usize,
    gravity: [f32; 3],
    boundary: Boundary,
    n: usize,
}

impl ResidentGpu {
    /// Construct an unbuilt resident-GPU resource (state is built lazily on the
    /// first step from the host `Atom` set). Exposed so examples/tests can wire
    /// the resident system into a minimal `App` without the full plugin group.
    pub fn new(ctx: Option<GpuContext>, window: usize, gravity: [f32; 3], boundary: Boundary) -> Self {
        Self { ctx, gpu: None, omega_aux: 0, primed: false, window: window.max(1), gravity, boundary, n: 0 }
    }

    fn build(&mut self, atoms: &Atom, registry: &AtomDataRegistry, mt: &MaterialTable) {
        let Some(ctx) = self.ctx.clone() else { return };
        let n = atoms.nlocal as usize;
        let dem = registry.expect::<DemAtom>("ResidentGpu::build");

        let radius: Vec<f32> = (0..n).map(|i| dem.radius[i] as f32).collect();
        let inv_inertia: Vec<f32> = (0..n).map(|i| dem.inv_inertia[i] as f32).collect();
        let inv_mass: Vec<f32> = (0..n).map(|i| 1.0 / atoms.mass[i] as f32).collect();
        let posf: Vec<[f32; 3]> = (0..n)
            .map(|i| [atoms.pos[i][0] as f32, atoms.pos[i][1] as f32, atoms.pos[i][2] as f32])
            .collect();
        let velf: Vec<[f32; 3]> = (0..n)
            .map(|i| [atoms.vel[i][0] as f32, atoms.vel[i][1] as f32, atoms.vel[i][2] as f32])
            .collect();
        let omf: Vec<[f32; 3]> = (0..n)
            .map(|i| [dem.omega[i][0] as f32, dem.omega[i][1] as f32, dem.omega[i][2] as f32])
            .collect();

        let r_max = radius.iter().copied().fold(0.0f32, f32::max).max(f32::MIN_POSITIVE);
        let grid = Grid::from_positions(&posf, 2.0 * r_max);
        let dt = atoms.dt as f32;

        let mut gs = GpuState::new(ctx, n, grid.total_cells);
        gs.set_params(dt, self.gravity);
        gs.set_state(&posf, &velf, &inv_mass, grid);
        let omega_aux = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega_aux, &inv_inertia);
        gs.set_aux_state(omega_aux, &omf);

        let (e_eff, beta, g_eff, mu) = gpu_scalars(mt);
        let cfg = GranularConfig { e_eff, beta, g_eff, mu, dt };
        gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega_aux, &radius, cfg)));
        gs.add_force_hook(Box::new(WallForce::new(
            &gs, omega_aux, &radius, &self.boundary, e_eff, beta, g_eff, mu, dt,
        )));

        self.gpu = Some(gs);
        self.omega_aux = omega_aux;
        self.n = n;
        self.primed = false;
    }
}

/// Resident step system: advance the device by one window of `K` steps, then sync
/// host `Atom`/`DemAtom` so I/O and diagnostics see the current state. The device
/// state is authoritative across calls — host arrays are NOT re-uploaded (that
/// would re-prime the force and corrupt the contact history).
pub fn gpu_granular_resident_step(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    mut res: ResMut<ResidentGpu>,
) {
    let n = atoms.nlocal as usize;
    if n == 0 || res.ctx.is_none() {
        return;
    }
    if res.gpu.is_none() || res.n != n {
        res.build(&atoms, &registry, &material_table);
    }
    let window = res.window;
    let omega_aux = res.omega_aux;
    let primed = res.primed;

    let (p, v, w) = {
        let gs = res.gpu.as_ref().unwrap();
        if primed {
            gs.run_steps_continue(window);
        } else {
            gs.run_steps(window);
        }
        gs.wait();
        (gs.download_pos(), gs.download_vel(), gs.download_aux_state(omega_aux))
    };
    res.primed = true;

    for i in 0..n {
        atoms.pos[i] = [p[i][0] as Real, p[i][1] as Real, p[i][2] as Real];
        atoms.vel[i] = [v[i][0] as Real, v[i][1] as Real, v[i][2] as Real];
    }
    let mut dem = registry.expect_mut::<DemAtom>("gpu_granular_resident_step");
    for i in 0..n {
        dem.omega[i] = [w[i][0] as f64, w[i][1] as f64, w[i][2] as f64];
    }
}

/// Device-resident granular stepper (roadmap step 1). Add this INSTEAD of
/// `VelocityVerletPlugin` + host force + `RotationalDynamicsPlugin`: it owns the
/// whole velocity-Verlet loop (integrate + contact + walls + gravity + rotation)
/// on the GPU for `window` steps per schedule tick. Falls back to a no-op if no
/// GPU adapter is present (use the CPU/host plugins on such machines).
pub struct GpuGranularResidentPlugin {
    /// Planar walls evaluated on-device (floor / box faces).
    pub boundary: Boundary,
    /// Constant body acceleration applied on-device (e.g. `[0,0,-9.81]`).
    pub gravity: [f32; 3],
    /// Steps advanced on-device per schedule tick. Larger = fewer host syncs =
    /// faster, but coarser host-visible output cadence.
    pub window: usize,
}

impl Plugin for GpuGranularResidentPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces", "integration"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles"]
    }

    fn build(&self, app: &mut App) {
        let ctx = GpuContext::new();
        if let Some(ref c) = ctx {
            println!("GpuGranularResident: resident granular stepper on GPU adapter: {}", c.adapter_info);
        } else {
            eprintln!("GpuGranularResident: no GPU adapter — resident step is a no-op");
        }
        app.add_resource(ResidentGpu {
            ctx,
            gpu: None,
            omega_aux: 0,
            primed: false,
            window: self.window.max(1),
            gravity: self.gravity,
            boundary: self.boundary.clone(),
            n: 0,
        });
        // One system that advances the device a whole window and writes back.
        app.add_update_system(
            gpu_granular_resident_step.label("gpu_granular_resident"),
            ParticleSimScheduleSet::Force,
        );
    }
}
