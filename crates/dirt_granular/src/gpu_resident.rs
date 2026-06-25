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

use std::any::TypeId;

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
        let cfg = GranularConfig::new(e_eff, beta, g_eff, mu, dt);
        gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega_aux, &radius, cfg)));
        gs.add_force_hook(Box::new(WallForce::new(
            &gs, omega_aux, &radius, &self.boundary, e_eff, beta, g_eff, mu, dt,
        )));

        self.gpu = Some(gs);
        self.omega_aux = omega_aux;
        self.n = n;
        self.primed = false;
    }

    /// Push the current host `Atom`/`DemAtom` state for the local bulk back onto the
    /// device, overwriting the resident buffers. Used by the coherence path when a
    /// host system mutated the state between windows (policy A): the caller then
    /// clears `primed` so the next `run_steps` re-primes the force (contact history
    /// resets — a small physics discontinuity, see coherence_plan.md). Positions are
    /// assumed to stay within the existing grid bounds (true for velocity/force
    /// edits; a teleporting write would also need a `set_grid`).
    fn reupload_locals(&self, atoms: &Atom, registry: &AtomDataRegistry) {
        let Some(gs) = self.gpu.as_ref() else { return };
        let n = self.n;
        let posf: Vec<[f32; 3]> = (0..n)
            .map(|i| [atoms.pos[i][0] as f32, atoms.pos[i][1] as f32, atoms.pos[i][2] as f32])
            .collect();
        let velf: Vec<[f32; 3]> = (0..n)
            .map(|i| [atoms.vel[i][0] as f32, atoms.vel[i][1] as f32, atoms.vel[i][2] as f32])
            .collect();
        let omf: Vec<[f32; 3]> = {
            let dem = registry.expect::<DemAtom>("ResidentGpu::reupload_locals");
            (0..n)
                .map(|i| [dem.omega[i][0] as f32, dem.omega[i][1] as f32, dem.omega[i][2] as f32])
                .collect()
        };
        gs.write_pos_slice(0, &posf);
        gs.write_vel_slice(0, &velf);
        gs.write_aux_slice(self.omega_aux, 0, &omf);
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
    mut coherence: Option<ResMut<CoherenceRegistry>>,
) {
    let n = atoms.nlocal as usize;
    if n == 0 || res.ctx.is_none() {
        return;
    }
    if res.gpu.is_none() || res.n != n {
        res.build(&atoms, &registry, &material_table);
    }

    // Coherence path (policy A): if a host system wrote `Atom` since the last
    // window, push it to the device and force a re-prime (contact history resets).
    // `take_host_dirty` is None when coherence is off → eager path below, unchanged.
    if let Some(c) = coherence.as_mut() {
        if c.take_host_dirty(TypeId::of::<Atom>()) {
            res.reupload_locals(&atoms, &registry);
            res.primed = false;
            if std::env::var("SIM_SUPPRESS_WARNINGS").is_err() {
                eprintln!(
                    "[coherence] resident GPU re-primed from host-modified state (contact history reset)"
                );
            }
        }
    }

    let window = res.window;
    let omega_aux = res.omega_aux;
    let primed = res.primed;

    {
        let gs = res.gpu.as_ref().unwrap();
        if primed {
            gs.run_steps_continue(window);
        } else {
            gs.run_steps(window);
        }
    }
    res.primed = true;

    if let Some(c) = coherence.as_mut() {
        // Lazy sync: the device is now authoritative. The host `Atom`/`DemAtom`
        // mirror is NOT downloaded here — the next host consumer triggers the pull
        // via ResidentMirrorBridge (scheduler-mediated). Saves the per-window
        // download whenever no host system reads the state this window.
        c.mark_device_dirty(TypeId::of::<Atom>());
    } else {
        // Eager path (coherence off): wait, download, and write the host mirror
        // back every window — the original step-1 behaviour.
        let (p, v, w) = {
            let gs = res.gpu.as_ref().unwrap();
            gs.wait();
            (gs.download_pos(), gs.download_vel(), gs.download_aux_state(omega_aux))
        };
        for i in 0..n {
            atoms.pos[i] = [p[i][0] as Real, p[i][1] as Real, p[i][2] as Real];
            atoms.vel[i] = [v[i][0] as Real, v[i][1] as Real, v[i][2] as Real];
        }
        let mut dem = registry.expect_mut::<DemAtom>("gpu_granular_resident_step");
        for i in 0..n {
            dem.omega[i] = [w[i][0] as f64, w[i][1] as f64, w[i][2] as f64];
        }
    }
}

/// Bridge that pulls the resident device state (pos/vel + omega) back into the
/// host `Atom`/`DemAtom` when a host consumer reads `Atom` while the mirror is
/// `DeviceDirty` (coherence_plan.md Phase 3). Borrows the resource cells it needs
/// by index; `resolve` fills those in once the resource table is final.
#[cfg(feature = "gpu_coherence")]
struct ResidentMirrorBridge {
    atom_idx: usize,
    registry_idx: usize,
    gpu_idx: usize,
}

#[cfg(feature = "gpu_coherence")]
impl ResidentMirrorBridge {
    fn unresolved() -> Self {
        ResidentMirrorBridge { atom_idx: usize::MAX, registry_idx: usize::MAX, gpu_idx: usize::MAX }
    }
}

/// Build a [`CoherenceRegistry`] preloaded with the resident GPU's `Atom` mirror
/// (pos/vel + omega bridge). Used by [`GpuGranularResidentPlugin`] and exposed so
/// examples/tests that wire the resident step manually can opt into coherence.
#[cfg(feature = "gpu_coherence")]
pub fn resident_coherence_registry() -> CoherenceRegistry {
    let mut reg = CoherenceRegistry::new();
    reg.register(TypeId::of::<Atom>(), Box::new(ResidentMirrorBridge::unresolved()));
    reg
}

#[cfg(feature = "gpu_coherence")]
impl MirrorBridge for ResidentMirrorBridge {
    fn resolve(&mut self, index: &std::collections::HashMap<TypeId, usize>) {
        self.atom_idx = index[&TypeId::of::<Atom>()];
        self.registry_idx = index[&TypeId::of::<AtomDataRegistry>()];
        self.gpu_idx = index[&TypeId::of::<ResidentGpu>()];
    }

    fn download(&self, resources: &[std::cell::RefCell<Box<dyn std::any::Any>>]) {
        let gpu_cell = resources[self.gpu_idx].borrow();
        let res = gpu_cell
            .downcast_ref::<ResidentGpu>()
            .expect("ResidentMirrorBridge: gpu_idx is not a ResidentGpu");
        let Some(gs) = res.gpu.as_ref() else { return };
        let n = res.n;
        let omega_aux = res.omega_aux;
        gs.wait();
        let p = gs.download_pos();
        let v = gs.download_vel();
        let w = gs.download_aux_state(omega_aux);

        {
            let mut atom_cell = resources[self.atom_idx].borrow_mut();
            let atoms = atom_cell
                .downcast_mut::<Atom>()
                .expect("ResidentMirrorBridge: atom_idx is not an Atom");
            for i in 0..n {
                atoms.pos[i] = [p[i][0] as Real, p[i][1] as Real, p[i][2] as Real];
                atoms.vel[i] = [v[i][0] as Real, v[i][1] as Real, v[i][2] as Real];
            }
        }
        let reg_cell = resources[self.registry_idx].borrow();
        let registry = reg_cell
            .downcast_ref::<AtomDataRegistry>()
            .expect("ResidentMirrorBridge: registry_idx is not an AtomDataRegistry");
        let mut dem = registry.expect_mut::<DemAtom>("ResidentMirrorBridge::download");
        for i in 0..n {
            dem.omega[i] = [w[i][0] as f64, w[i][1] as f64, w[i][2] as f64];
        }
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

        // Coherence (coherence_plan.md Phase 3): register the Atom mirror so host
        // systems added to this resident config sync transparently instead of
        // silently dropping. Off by default — the resident step then uses the eager
        // per-window download and this registry is never created.
        #[cfg(feature = "gpu_coherence")]
        {
            app.add_resource(resident_coherence_registry());
            println!("GpuGranularResident: host<->device coherence enabled (Atom mirror)");
        }
    }
}
