//! GPU offload of the particle-particle contact force (Force phase).
//!
//! [`GpuGranularForcePlugin`] is a drop-in replacement for
//! [`HertzMindlinContactPlugin`](crate::HertzMindlinContactPlugin) that computes
//! the Hertz-normal + Mindlin-tangential contact force (plus rotational torque
//! and body-force gravity) on the GPU via soil's resident loop
//! (`dirt_gpu`/`soil_gpu`), instead of the CPU `hertz_mindlin_contact_force`.
//! CPU integration and quaternion rotation (`VelocityVerletPlugin`,
//! `RotationalDynamicsPlugin`) are unchanged — this is a force-phase offload
//! (the LAMMPS GPU-package model), not full device residency.
//!
//! ## Scope & composition
//!
//! - Physics: **plain Hertz-Mindlin only** — no rolling/twisting friction,
//!   cohesion, or surface energy (the GPU kernel omits them). Use the CPU plugin
//!   for those material models.
//! - The GPU computes only the particle-particle **contact** force/torque and
//!   **accumulates** (`+=`) into `Atom::force` / `DemAtom::torque`, exactly like
//!   the CPU `hertz_mindlin_contact_force` it replaces (forces are zeroed each
//!   step by soil_core's `zero_all_forces` at `PostInitialIntegration`). So it
//!   **composes** with every other force contributor — gravity, `WallPlugin`,
//!   and any custom CPU force fix all just add their force on top. Because force
//!   is synced host↔device each step in this offload model, no `DualBuffer`
//!   coherence is needed for correctness — that is only a later optimization for
//!   collapsing the per-step transfer (full device residency).
//! - Contact (tangential spring) history lives on-device and persists across
//!   steps; it resets only when the local atom count changes (e.g. insertion).

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_gpu::{
    GpuContext, GpuState, GranularConfig, GranularForce, Grid,
};
use soil_core::{Atom, AtomDataRegistry, ParticleSimScheduleSet};

/// Resource holding the GPU device and the resident granular force state. Built
/// lazily on the first force step (once atoms exist), rebuilt if the local atom
/// count changes.
pub struct GpuGranular {
    ctx: Option<GpuContext>,
    state: Option<GpuState>,
    omega_aux: usize,
    grid: Grid,
    n: usize,
    /// Cached per-atom inverse mass (static across steps; re-uploaded each step).
    inv_mass: Vec<f32>,
}

impl GpuGranular {
    fn new(ctx: Option<GpuContext>) -> Self {
        Self {
            ctx,
            state: None,
            omega_aux: 0,
            grid: Grid { n: [1; 3], origin: [0.0; 3], bin_size: 1.0, total_cells: 1 },
            n: 0,
            inv_mass: Vec::new(),
        }
    }

    /// (Re)build the resident GPU state + contact hook for the current atom set.
    /// Resets device contact history (acceptable on an atom-count change).
    fn build(&mut self, atoms: &Atom, registry: &AtomDataRegistry, mt: &MaterialTable) {
        let Some(ctx) = self.ctx.clone() else { return };
        let n = atoms.nlocal as usize;
        let dem = registry.expect::<DemAtom>("GpuGranular::build");

        let radius: Vec<f32> = (0..n).map(|i| dem.radius[i] as f32).collect();
        let inv_inertia: Vec<f32> = (0..n).map(|i| dem.inv_inertia[i] as f32).collect();
        let inv_mass: Vec<f32> = (0..n).map(|i| 1.0 / atoms.mass[i] as f32).collect();
        let posf: Vec<[f32; 3]> =
            (0..n).map(|i| [atoms.pos[i][0] as f32, atoms.pos[i][1] as f32, atoms.pos[i][2] as f32]).collect();
        let velf: Vec<[f32; 3]> =
            (0..n).map(|i| [atoms.vel[i][0] as f32, atoms.vel[i][1] as f32, atoms.vel[i][2] as f32]).collect();
        let omf: Vec<[f32; 3]> =
            (0..n).map(|i| [dem.omega[i][0] as f32, dem.omega[i][1] as f32, dem.omega[i][2] as f32]).collect();

        let r_max = radius.iter().copied().fold(0.0f32, f32::max).max(f32::MIN_POSITIVE);
        // soil's Grid takes the literal contact cutoff (sum of radii = 2*r_max).
        let grid = Grid::from_positions(&posf, 2.0 * r_max);
        let dt = atoms.dt as f32;

        let mut gs = GpuState::new(ctx, n, grid.total_cells);
        // Gravity off: this hook contributes only the contact force, which is
        // accumulated onto the host force alongside gravity / walls / other fixes.
        gs.set_params(dt, [0.0, 0.0, 0.0]);
        gs.set_state(&posf, &velf, &inv_mass, grid);
        let omega_aux = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega_aux, &inv_inertia);
        gs.set_aux_state(omega_aux, &omf);

        let (e_eff, beta, g_eff, mu) = gpu_scalars(mt);
        let cfg = GranularConfig { e_eff, beta, g_eff, mu, dt };
        gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega_aux, &radius, cfg)));

        self.state = Some(gs);
        self.omega_aux = omega_aux;
        self.grid = grid;
        self.n = n;
        self.inv_mass = inv_mass;
    }
}

/// Extract the (single-material) plain Hertz-Mindlin scalars the GPU kernel uses.
fn gpu_scalars(mt: &MaterialTable) -> (f32, f32, f32, f32) {
    (
        mt.e_eff_ij[0][0] as f32,
        mt.beta_ij[0][0] as f32,
        mt.g_eff_ij[0][0] as f32,
        mt.friction_ij[0][0] as f32,
    )
}

/// Force-phase system: upload current state, evaluate contact force + torque on
/// the GPU (no integration), write force/torque back to the host so the CPU
/// integration and rotation phases consume them.
fn gpu_granular_force(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    mut g: ResMut<GpuGranular>,
) {
    let n = atoms.nlocal as usize;
    if n == 0 || g.ctx.is_none() {
        return;
    }
    if g.state.is_none() || g.n != n {
        g.build(&atoms, &registry, &material_table);
    }
    let grid = g.grid;
    let omega_aux = g.omega_aux;

    let posf: Vec<[f32; 3]> =
        (0..n).map(|i| [atoms.pos[i][0] as f32, atoms.pos[i][1] as f32, atoms.pos[i][2] as f32]).collect();
    let velf: Vec<[f32; 3]> =
        (0..n).map(|i| [atoms.vel[i][0] as f32, atoms.vel[i][1] as f32, atoms.vel[i][2] as f32]).collect();

    let mut dem = registry.expect_mut::<DemAtom>("gpu_granular_force");
    let omf: Vec<[f32; 3]> =
        (0..n).map(|i| [dem.omega[i][0] as f32, dem.omega[i][1] as f32, dem.omega[i][2] as f32]).collect();

    let (force, torque) = {
        let gs = g.state.as_ref().unwrap();
        gs.set_state(&posf, &velf, &g.inv_mass, grid); // re-upload moved positions/velocities
        gs.set_aux_state(omega_aux, &omf);
        gs.eval_force_once();
        (gs.download_force(), gs.download_aux_rate(omega_aux))
    };

    // Accumulate (like the CPU hertz path) so the GPU contact force composes with
    // gravity, walls, and any other force contributor.
    for i in 0..n {
        atoms.force[i][0] += force[i][0] as soil_core::Accum;
        atoms.force[i][1] += force[i][1] as soil_core::Accum;
        atoms.force[i][2] += force[i][2] as soil_core::Accum;
        dem.torque[i][0] += torque[i][0] as f64;
        dem.torque[i][1] += torque[i][1] as f64;
        dem.torque[i][2] += torque[i][2] as f64;
    }
}

/// GPU drop-in for [`HertzMindlinContactPlugin`](crate::HertzMindlinContactPlugin):
/// computes the plain Hertz-Mindlin particle-particle contact force on the GPU at
/// the `Force` phase, accumulating into the host force/torque so it composes with
/// gravity, walls, and other force fixes. Falls back to the CPU
/// `hertz_mindlin_contact_force` if no GPU adapter is available, so it is safe to
/// add on any machine. See the module docs for scope.
#[derive(Default)]
pub struct GpuGranularForcePlugin;

impl Plugin for GpuGranularForcePlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        match GpuContext::new() {
            Some(ctx) => {
                println!("GpuGranularForce: contact force on GPU adapter: {}", ctx.adapter_info);
                app.add_resource(GpuGranular::new(Some(ctx)));
                app.add_update_system(
                    gpu_granular_force.label("hertz_mindlin_contact"),
                    ParticleSimScheduleSet::Force,
                );
            }
            None => {
                eprintln!("GpuGranularForce: no GPU adapter — falling back to CPU hertz_mindlin");
                app.add_update_system(
                    crate::contact::hertz_mindlin_contact_force.label("hertz_mindlin_contact"),
                    ParticleSimScheduleSet::Force,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contact::hertz_mindlin_contact_force;
    use crate::rotational::{final_rotation, initial_rotation};
    use crate::tangential::ContactHistoryStore;
    use dirt_test_utils::make_material_table;
    use soil_core::Neighbor;
    use soil_verlet::{final_integration, initial_integration};

    const DENSITY: f64 = 2500.0;
    const R: f64 = 0.05;

    fn mass() -> f64 {
        DENSITY * 4.0 / 3.0 * std::f64::consts::PI * R * R * R
    }

    /// Two grains overlapping head-on along x, with tangential (y) approach and
    /// spin — exercises normal + tangential + torque. Gravity off, no walls, so
    /// the only force is particle-particle contact (the GPU plugin's domain).
    fn make_state() -> (Atom, DemAtom) {
        let m = mass();
        let inv_inertia = 1.0 / (0.4 * m * R * R);
        let mut atom = Atom::new();
        atom.dt = 1.0e-5;
        atom.push_test_atom(0, [-0.049, 0.0, 0.0], R, m);
        atom.push_test_atom(1, [0.049, 0.0, 0.0], R, m);
        atom.vel[0] = [0.05, 0.06, 0.0];
        atom.vel[1] = [-0.05, -0.06, 0.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut dem = DemAtom::new();
        for (i, om) in [[0.0, 0.0, 0.8], [0.0, 0.0, -0.5]].into_iter().enumerate() {
            dem.radius.push(R);
            dem.density.push(DENSITY);
            dem.inv_inertia.push(inv_inertia);
            dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
            dem.omega.push(om);
            dem.ang_mom.push([0.0; 3]);
            dem.torque.push([0.0; 3]);
            dem.body_id.push(0.0);
            let _ = i;
        }
        (atom, dem)
    }

    fn full_neighbor_list(n: usize) -> Neighbor {
        let mut nb = Neighbor::new();
        nb.newton = false;
        let mut offsets = vec![0u32];
        let mut indices = Vec::new();
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    indices.push(j as u32);
                }
            }
            offsets.push(indices.len() as u32);
        }
        nb.neighbor_offsets = offsets;
        nb.neighbor_indices = indices;
        nb
    }

    // Zero force + torque each step (mirrors soil_core's zero_all_forces at
    // PostInitialIntegration). Both the CPU hertz and the GPU contact accumulate
    // onto the zeroed arrays.
    fn zero_force_torque(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>) {
        let n = atoms.nlocal as usize;
        atoms.force[..n].fill([0.0; 3]);
        let mut dem = registry.expect_mut::<DemAtom>("zero_force_torque");
        for i in 0..n {
            dem.torque[i] = [0.0; 3];
        }
    }

    // A second force contributor (a constant body force), added to BOTH apps to
    // prove the GPU contact force composes additively with other force fixes.
    fn body_force(mut atoms: ResMut<Atom>) {
        let n = atoms.nlocal as usize;
        for i in 0..n {
            let m = atoms.mass[i] as f64;
            atoms.force[i][2] += (m * -9.81) as soil_core::Accum;
        }
    }

    #[test]
    fn gpu_plugin_matches_cpu_hertz_contact_only() {
        let Some(ctx) = GpuContext::new() else {
            eprintln!("no GPU adapter; skipping GPU-vs-CPU plugin comparison");
            return;
        };

        // ── CPU app: full schedule with CPU hertz at the Force phase ──────────
        let mut cpu = App::new();
        {
            let (atom, dem) = make_state();
            let n = atom.nlocal as usize;
            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(ContactHistoryStore::new());
            cpu.add_resource(atom);
            cpu.add_resource(registry);
            cpu.add_resource(make_material_table());
            cpu.add_resource(full_neighbor_list(n));
        }
        cpu.add_update_system(initial_integration, ParticleSimScheduleSet::InitialIntegration);
        cpu.add_update_system(initial_rotation, ParticleSimScheduleSet::InitialIntegration);
        cpu.add_update_system(zero_force_torque, ParticleSimScheduleSet::PostInitialIntegration);
        cpu.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        cpu.add_update_system(body_force, ParticleSimScheduleSet::Force);
        cpu.add_update_system(final_integration, ParticleSimScheduleSet::FinalIntegration);
        cpu.add_update_system(final_rotation, ParticleSimScheduleSet::FinalIntegration);
        cpu.organize_systems();

        // ── GPU app: same schedule, GPU contact force at the Force phase ──────
        let mut gpu = App::new();
        {
            let (atom, dem) = make_state();
            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            gpu.add_resource(atom);
            gpu.add_resource(registry);
            gpu.add_resource(make_material_table());
            gpu.add_resource(GpuGranular::new(Some(ctx)));
        }
        gpu.add_update_system(initial_integration, ParticleSimScheduleSet::InitialIntegration);
        gpu.add_update_system(initial_rotation, ParticleSimScheduleSet::InitialIntegration);
        gpu.add_update_system(zero_force_torque, ParticleSimScheduleSet::PostInitialIntegration);
        gpu.add_update_system(gpu_granular_force, ParticleSimScheduleSet::Force);
        gpu.add_update_system(body_force, ParticleSimScheduleSet::Force);
        gpu.add_update_system(final_integration, ParticleSimScheduleSet::FinalIntegration);
        gpu.add_update_system(final_rotation, ParticleSimScheduleSet::FinalIntegration);
        gpu.organize_systems();

        let steps = 150;
        for _ in 0..steps {
            cpu.run();
            gpu.run();
        }

        let ca = cpu.get_resource_ref::<Atom>().unwrap();
        let ga = gpu.get_resource_ref::<Atom>().unwrap();
        let mut max_diff = 0.0f64;
        for i in 0..2 {
            for d in 0..3 {
                let c = ca.pos[i][d] as f64;
                let g = ga.pos[i][d] as f64;
                max_diff = max_diff.max((c - g).abs());
                assert!(
                    (c - g).abs() < 1e-3,
                    "atom {i} pos[{d}]: cpu={c} gpu={g} (diff {})",
                    (c - g).abs()
                );
            }
        }
        eprintln!("GpuGranularForcePlugin vs CPU hertz: max pos diff = {max_diff:.2e} over {steps} steps");
    }
}
