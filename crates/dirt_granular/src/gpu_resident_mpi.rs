//! Ghost-aware **windowed-resident** GPU stepper for MPI (roadmap step 1 × MPI).
//!
//! Fuses milestone-1's ghost-awareness (`gpu.rs`: upload local+ghost, keep locals)
//! with step-1 residency (`gpu_resident.rs`: advance `window` Verlet steps on-device
//! per schedule tick). At each tick the host `Atom` already holds fresh ghosts (the
//! schedule's `forward_comm`/`borders` ran in the Exchange/Neighbor phases before
//! `Force`), so this system uploads **local+ghost**, advances `window` steps on the
//! device, and writes **locals** back.
//!
//! ## Correctness boundary (the step-2b seam)
//!
//! The schedule forward-comms ghost positions **every step** (`CommState::
//! CommunicateOnly`). This system only refreshes ghosts once per *tick* (= once per
//! `window` device steps). So:
//! - `window == 1` → ghosts are fresh every device step → correct (matches CPU /
//!   1-rank within the f32 band), but there is no residency win (a host round-trip
//!   every step, like milestone-1).
//! - `window > 1` → ghosts are frozen for `window-1` steps while their owner rank
//!   advances them → boundary-local contacts see a stale ghost trajectory → the
//!   result diverges, growing with `window`.
//!
//! Realising residency *under MPI* therefore requires GPU-resident halos (step 2b):
//! pack/unpack `forward_comm` from device buffers so ghost positions refresh
//! on-device each step without a host round-trip. This module exists to demonstrate
//! that boundary; for single-rank residency use `GpuGranularResidentPlugin`.

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_gpu::{Boundary, GpuContext, GpuState, GranularConfig, GranularForce, Grid, WallForce};
use soil_core::{Atom, AtomDataRegistry, CommResource, ParticleSimScheduleSet, Real};

use crate::gpu::gpu_scalars;

/// Resident GPU state for the MPI path. Rebuilt when the local+ghost count changes
/// (every neighbour rebuild — `forward_comm` keeps the count fixed between them).
pub struct ResidentMpiGpu {
    ctx: Option<GpuContext>,
    gpu: Option<GpuState>,
    omega_aux: usize,
    window: usize,
    gravity: [f32; 3],
    boundary: Boundary,
    /// local+ghost count the device buffers are sized for.
    nall: usize,
}

impl ResidentMpiGpu {
    /// (Re)build device state sized for the current local+ghost set, with the
    /// contact + wall force hooks. Uploads local+ghost positions/velocities.
    fn build(&mut self, atoms: &Atom, registry: &AtomDataRegistry, mt: &MaterialTable) {
        let Some(ctx) = self.ctx.clone() else { return };
        let nall = atoms.len();
        let dem = registry.expect::<DemAtom>("ResidentMpiGpu::build");

        let radius: Vec<f32> = (0..nall).map(|i| dem.radius[i] as f32).collect();
        let inv_inertia: Vec<f32> = (0..nall).map(|i| dem.inv_inertia[i] as f32).collect();
        let inv_mass: Vec<f32> = (0..nall).map(|i| 1.0 / atoms.mass[i] as f32).collect();
        let posf = upload_vec3(&atoms.pos, nall);
        let velf = upload_vec3(&atoms.vel, nall);
        let omf = upload_vec3(&dem.omega, nall);

        let r_max = radius.iter().copied().fold(0.0f32, f32::max).max(f32::MIN_POSITIVE);
        let grid = Grid::from_positions(&posf, 2.0 * r_max);
        let dt = atoms.dt as f32;

        let mut gs = GpuState::new(ctx, nall, grid.total_cells);
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
        self.nall = nall;
    }
}

fn upload_vec3<T: Copy + Into<f64>>(src: &[[T; 3]], n: usize) -> Vec<[f32; 3]> {
    (0..n)
        .map(|i| [src[i][0].into() as f32, src[i][1].into() as f32, src[i][2].into() as f32])
        .collect()
}

/// Resident-MPI step (Force phase): upload local+ghost (fresh from this tick's
/// host forward_comm), advance `window` steps on-device, write LOCALS back.
/// Re-primes every tick because the ghost set / positions changed.
pub fn gpu_resident_mpi_step(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    mut res: ResMut<ResidentMpiGpu>,
) {
    let nlocal = atoms.nlocal as usize;
    let nall = atoms.len();
    if nlocal == 0 || res.ctx.is_none() {
        return;
    }
    // Rebuild device buffers when the local+ghost count changes (a rebuild step);
    // otherwise re-upload the fresh state into the existing buffers.
    if res.gpu.is_none() || res.nall != nall {
        res.build(&atoms, &registry, &material_table);
    } else {
        // Step 2b (GPU-resident halos): the local bulk stays RESIDENT on-device and
        // is NEVER re-uploaded — it's bit-identical to the host copy captured by last
        // tick's download. Only the ghost slice [nlocal..nall], refreshed by the
        // host's forward_comm this tick, is written to the device. The grid is
        // computed CPU-side from the host positions (current locals + fresh ghosts);
        // Grid::from_positions does no GPU upload.
        let posf = upload_vec3(&atoms.pos, nall);
        let r_max = (0..nall)
            .map(|i| registry.expect::<DemAtom>("step").radius[i] as f32)
            .fold(0.0f32, f32::max)
            .max(f32::MIN_POSITIVE);
        let grid = Grid::from_positions(&posf, 2.0 * r_max);
        let omega_aux = res.omega_aux;
        let gs = res.gpu.as_ref().unwrap();
        let nghost = nall - nlocal;
        if nghost > 0 {
            let gpos = upload_vec3(&atoms.pos[nlocal..nall], nghost);
            let gvel = upload_vec3(&atoms.vel[nlocal..nall], nghost);
            let gom = {
                let dem = registry.expect::<DemAtom>("step");
                upload_vec3(&dem.omega[nlocal..nall], nghost)
            };
            gs.write_pos_slice(nlocal, &gpos);
            gs.write_vel_slice(nlocal, &gvel);
            gs.write_aux_slice(omega_aux, nlocal, &gom);
        }
        gs.set_grid(grid);
    }

    let window = res.window;
    let omega_aux = res.omega_aux;
    let (p, v, w) = {
        let gs = res.gpu.as_ref().unwrap();
        gs.run_steps(window); // re-prime + window steps (ghosts changed this tick)
        gs.wait();
        (gs.download_pos(), gs.download_vel(), gs.download_aux_state(omega_aux))
    };

    // Write LOCALS back; ghosts are re-derived by the host's forward_comm next tick.
    for i in 0..nlocal {
        atoms.pos[i] = [p[i][0] as Real, p[i][1] as Real, p[i][2] as Real];
        atoms.vel[i] = [v[i][0] as Real, v[i][1] as Real, v[i][2] as Real];
    }
    let mut dem = registry.expect_mut::<DemAtom>("gpu_resident_mpi_step");
    for i in 0..nlocal {
        dem.omega[i] = [w[i][0] as f64, w[i][1] as f64, w[i][2] as f64];
    }
}

/// Ghost-aware windowed-resident GPU stepper for MPI. Add INSTEAD of host Verlet +
/// force + rotation (it owns the loop), like `GpuGranularResidentPlugin`, but it
/// uploads local+ghost each tick so cross-rank contacts compose. See module docs
/// for the `window` correctness boundary (window=1 correct; window>1 needs step-2b
/// GPU-resident halos).
pub struct GpuGranularResidentMpiPlugin {
    pub boundary: Boundary,
    pub gravity: [f32; 3],
    pub window: usize,
}

impl Plugin for GpuGranularResidentMpiPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }
    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces", "integration"]
    }
    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles", "neighbor_list"]
    }
    fn build(&self, app: &mut App) {
        // Step 3: bind one GPU per rank — rank r uses adapter r % num_adapters.
        // No-op on a single-GPU machine (every rank gets device 0).
        let local_rank = app
            .get_resource_ref::<CommResource>()
            .map(|c| c.rank() as usize)
            .unwrap_or(0);
        let ctx = GpuContext::new_for_rank(local_rank);
        if let Some(ref c) = ctx {
            println!("GpuGranularResidentMpi: rank {local_rank} window={} on GPU {}", self.window, c.adapter_info);
        } else {
            eprintln!("GpuGranularResidentMpi: no GPU adapter — resident step is a no-op");
        }
        app.add_resource(ResidentMpiGpu {
            ctx,
            gpu: None,
            omega_aux: 0,
            window: self.window.max(1),
            gravity: self.gravity,
            boundary: self.boundary.clone(),
            nall: 0,
        });
        app.add_update_system(
            gpu_resident_mpi_step.label("gpu_resident_mpi"),
            ParticleSimScheduleSet::Force,
        );
    }
}
