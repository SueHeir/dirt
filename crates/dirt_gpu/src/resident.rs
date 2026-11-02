//! [`GpuResidentPlugin`]: the whole step on the GPU, for a window of steps at a time.
//!
//! Integration, neighbor binning, contact, planar walls and gravity all run on the device
//! for `window` steps per schedule tick, with nothing crossing to the CPU in between. So
//! everything else in the schedule runs once per window:
//!
//! - A `[[run]]` stage's `steps` count ticks, and the simulated steps are `steps * window`.
//! - Thermo, dumps and CPU fixes see the state every `window` steps. Reading it downloads
//!   positions, velocities and spins, through `soil_gpu`'s coherence bridge.
//! - A CPU system that changes those is uploaded before the next window, which then starts
//!   by evaluating the force afresh.
//! - When nothing changed, windows join without re-evaluating the force, so the window
//!   length changes the trajectory only by rounding (the cell grid can shift between
//!   windows, which reorders the sums).
//! - If the CPU re-sorts atoms or adds some, the springs follow their atoms by tag.
//!
//! It runs on one rank, for plain Hertz-Mindlin on one material, with gravity, periodic
//! axes, and fixed infinite planar walls. Anything else panics at the first step. Use
//! [`GpuContactPlugin`](crate::GpuContactPlugin) for those: under MPI, ghosts would go
//! stale inside a window, and other physics has no kernel here.
//!
//! Three things do not carry through. Springs are not written to restart files. Periodic
//! image counters are not updated, so unwrapped positions are wrong. Quaternions are not
//! integrated, so a run needing particle orientation should use the CPU rotation path.

use std::collections::HashMap;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_schedule::CONTACT_FORCE;
use dirt_wall::Walls;
use grass_app::prelude::*;
use grass_gpu::GpuDevice;
use grass_mpi::CommResource;
use grass_scheduler::prelude::*;
use grass_scheduler::CoherenceRegistry;
use soil_core::{Atom, AtomDataRegistry, Domain, ParticleSimScheduleSet};
use soil_fixes::GravityConfig;
use soil_gpu::neighbor_slots::{MAX_CONTACTS, PAYLOAD, SENTINEL};
use soil_gpu::resident::to_gpu;
use soil_gpu::{
    AtomColumnSync, Boundary, GpuReal, GpuState, Grid, Plane, ResidentAtoms, ResidentAtomsPlugin,
};

use crate::contact::require_supported;
use crate::granular_force::{GranularColumns, GranularConfig, GranularForce};
use crate::wall_force::{WallColumns, WallForce, MAX_WALLS};

const CELLS_PER_ATOM: usize = 8;

/// The rest of the resident state, beside `soil_gpu`'s [`ResidentAtoms`].
struct ResidentDem {
    gpu: GpuDevice,
    window: usize,
    omega: usize,
    granular: Option<GranularColumns>,
    walls: Option<WallColumns>,
    row_tags: Vec<u32>,
    primed: bool,
}

fn omega_to_host(registry: &AtomDataRegistry, omega: &[GpuReal]) {
    let mut dem = registry.expect_mut::<DemAtom>("omega_to_host");
    for (host, gpu) in dem.omega.iter_mut().zip(omega.chunks_exact(3)) {
        *host = [gpu[0] as f64, gpu[1] as f64, gpu[2] as f64];
    }
}

fn omega_from_host(registry: &AtomDataRegistry, n: usize) -> Vec<GpuReal> {
    let dem = registry.expect::<DemAtom>("omega_from_host");
    dem.omega[..n]
        .iter()
        .flat_map(|w| [w[0] as GpuReal, w[1] as GpuReal, w[2] as GpuReal])
        .collect()
}

/// The planes the kernels can run, or a panic naming what they cannot.
fn planes(walls: Option<&Walls>) -> Boundary {
    let mut boundary = Boundary::new();
    let Some(walls) = walls else { return boundary };
    let refuse = |what: &str| -> ! {
        panic!("GpuResidentPlugin runs fixed infinite planar walls, and this run has {what}. Use GpuContactPlugin.")
    };
    if walls.cylinder_active.iter().any(|&a| a) {
        refuse("a cylinder wall");
    }
    if walls.sphere_active.iter().any(|&a| a) {
        refuse("a sphere wall");
    }
    for (plane, _) in walls.planes.iter().zip(&walls.active).filter(|(_, &on)| on) {
        let bounds = [
            plane.bound_x_low,
            plane.bound_x_high,
            plane.bound_y_low,
            plane.bound_y_high,
            plane.bound_z_low,
            plane.bound_z_high,
        ];
        if bounds.iter().any(|b| b.is_finite()) {
            refuse("a bounded plane");
        }
        if plane.velocity != [0.0; 3] {
            refuse("a moving plane");
        }
        boundary.push(Plane::new(
            [
                plane.point_x as GpuReal,
                plane.point_y as GpuReal,
                plane.point_z as GpuReal,
            ],
            [
                plane.normal_x as GpuReal,
                plane.normal_y as GpuReal,
                plane.normal_z as GpuReal,
            ],
        ));
    }
    assert!(
        boundary.len() <= MAX_WALLS,
        "GpuResidentPlugin holds at most {MAX_WALLS} walls"
    );
    boundary
}

/// Remaps per-row state from `old` row order to `new`, by tag. Rows with no old
/// counterpart get `empty`.
fn remap_rows<T: Copy>(values: &[T], stride: usize, old: &[u32], new: &[u32], empty: T) -> Vec<T> {
    let from: HashMap<u32, usize> = old.iter().enumerate().map(|(i, &t)| (t, i)).collect();
    let mut out = vec![empty; new.len() * stride];
    for (row, tag) in new.iter().enumerate() {
        if let Some(&source) = from.get(tag) {
            out[row * stride..(row + 1) * stride]
                .copy_from_slice(&values[source * stride..(source + 1) * stride]);
        }
    }
    out
}

fn grid_for(domain: &Domain, pos: &[GpuReal], cutoff: GpuReal, max_cells: usize) -> Grid {
    Grid::for_box(
        domain.boundaries_low.map(|x| x as GpuReal),
        domain.boundaries_high.map(|x| x as GpuReal),
        [
            domain.is_periodic(0),
            domain.is_periodic(1),
            domain.is_periodic(2),
        ],
        pos,
        cutoff,
        max_cells,
    )
}

fn build(
    gpu: &GpuDevice,
    atoms: &Atom,
    mt: &MaterialTable,
    domain: &Domain,
    boundary: &Boundary,
    gravity: [GpuReal; 3],
    grid: &Grid,
    n: usize,
) -> (GpuState, usize, GranularColumns, Option<WallColumns>) {
    let dt = atoms.dt as GpuReal;
    let mut state = GpuState::new(gpu.clone(), n, *grid);
    state.set_params(dt, gravity);
    let periodic = [0, 1, 2].map(|d| domain.is_periodic(d));
    let lengths = [0, 1, 2].map(|d| {
        if periodic[d] {
            domain.size[d] as GpuReal
        } else {
            0.0
        }
    });
    if periodic.contains(&true) {
        state.set_box(
            lengths,
            domain.boundaries_low.map(|x| x as GpuReal),
            0.0,
            0.0,
        );
    }
    let omega = state.add_aux_dof();
    let (e, beta, g, mu) = (
        mt.e_eff_ij[0][0] as GpuReal,
        mt.beta_ij[0][0] as GpuReal,
        mt.g_eff_ij[0][0] as GpuReal,
        mt.friction_ij[0][0] as GpuReal,
    );
    let mut config = GranularConfig::new(e, beta, g, mu, dt);
    (config.lx, config.ly, config.lz) = (lengths[0], lengths[1], lengths[2]);
    let zeros = vec![0.0; n];
    let contact = GranularForce::new(&state, grid, omega, &zeros, config);
    let granular = contact.columns();
    state.add_force_hook(Box::new(contact));
    // Walls come after contact, since contact is what starts each torque at zero.
    let walls = (!boundary.is_empty()).then(|| {
        let hook = WallForce::new(&state, omega, &zeros, boundary, e, beta, g, mu, dt);
        let columns = hook.columns();
        state.add_force_hook(Box::new(hook));
        columns
    });
    (state, omega, granular, walls)
}

#[allow(clippy::too_many_arguments)]
fn gpu_resident_step(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    mt: Res<MaterialTable>,
    domain: Res<Domain>,
    walls: Option<Res<Walls>>,
    gravity: Option<Res<GravityConfig>>,
    mut resident: ResMut<ResidentAtoms>,
    mut dem_state: ResMut<ResidentDem>,
    mut coherence: ResMut<CoherenceRegistry>,
) {
    let n = atoms.nlocal() as usize;
    if n == 0 {
        return;
    }
    let host_dirty = coherence.take_host_dirty(std::any::TypeId::of::<Atom>());
    let tags = &atoms.tag[..n];
    let rows_changed = dem_state.row_tags != tags;
    let dem = &mut *dem_state;
    let max_cells = n * CELLS_PER_ATOM;

    if rows_changed {
        let pos = to_gpu(&atoms.pos, n);
        let cpu_dem = registry.expect::<DemAtom>("gpu_resident_step");
        require_supported(&mt, &cpu_dem, n);
        assert!(
            domain.tilt == [0.0; 3],
            "GpuResidentPlugin does not run a sheared or tilted box. Use GpuContactPlugin."
        );

        // Springs, read from the old rows before anything is replaced.
        let old = resident.state.as_ref().map(|_| {
            let rows = dem.row_tags.len();
            let (partner, payload) = dem
                .granular
                .as_ref()
                .expect("built with state")
                .read_slots(rows);
            let wall = dem.walls.as_ref().map(|w| w.read_springs(rows));
            (partner, payload, wall)
        });

        let radius: Vec<GpuReal> = cpu_dem.radius[..n].iter().map(|&r| r as GpuReal).collect();
        let cutoff = 2.0
            * radius
                .iter()
                .copied()
                .fold(GpuReal::MIN_POSITIVE, GpuReal::max);
        let grid = grid_for(&domain, &pos, cutoff, max_cells);
        if resident.state.is_none() {
            let g = gravity.as_ref().map_or([0.0; 3], |g| {
                [g.gx as GpuReal, g.gy as GpuReal, g.gz as GpuReal]
            });
            let boundary = planes(walls.as_deref());
            let (state, omega, granular, wall_columns) =
                build(&dem.gpu, &atoms, &mt, &domain, &boundary, g, &grid, n);
            resident.state = Some(state);
            resident.columns = vec![AtomColumnSync {
                aux: omega,
                to_host: omega_to_host,
                from_host: omega_from_host,
            }];
            dem.omega = omega;
            dem.granular = Some(granular);
            dem.walls = wall_columns;
        }

        let omega = dem.omega;
        let state = resident.state.as_mut().expect("just built");
        state.set_active(n, n, grid);
        state.write_pos(&pos);
        state.write_vel(&to_gpu(&atoms.vel, n));
        let inv_mass: Vec<GpuReal> = atoms.mass[..n]
            .iter()
            .map(|&m| 1.0 / m as GpuReal)
            .collect();
        state.write_inv_mass(&inv_mass);
        let inv_inertia: Vec<GpuReal> = cpu_dem.inv_inertia[..n]
            .iter()
            .map(|&x| x as GpuReal)
            .collect();
        state.set_aux_inv_coeff(omega, &inv_inertia);
        drop(cpu_dem);
        state.set_aux_state(omega, &omega_from_host(&registry, n));

        let granular = dem.granular.as_ref().expect("built with state");
        granular.write_radius(&radius);
        granular.write_tag(tags);
        if let Some(wall) = &dem.walls {
            wall.write_radius(&radius);
        }
        let (partner, payload, wall) = match old {
            Some((partner, payload, wall)) => (
                remap_rows(&partner, MAX_CONTACTS, &dem.row_tags, tags, SENTINEL),
                remap_rows(&payload, MAX_CONTACTS * PAYLOAD, &dem.row_tags, tags, 0.0),
                wall.map(|w| remap_rows(&w, MAX_WALLS * 3, &dem.row_tags, tags, 0.0)),
            ),
            None => (
                vec![SENTINEL; n * MAX_CONTACTS],
                vec![0.0; n * MAX_CONTACTS * PAYLOAD],
                dem.walls.as_ref().map(|_| vec![0.0; n * MAX_WALLS * 3]),
            ),
        };
        granular.write_slots(n, &partner, &payload);
        if let (Some(columns), Some(springs)) = (&dem.walls, wall) {
            columns.write_springs(&springs);
        }
        dem.row_tags.clear();
        dem.row_tags.extend_from_slice(tags);
        dem.primed = false;
    } else if host_dirty {
        // The CPU touched the atoms: upload only if it really changed them, so that plain
        // bookkeeping does not disturb the trajectory.
        let omega = dem.omega;
        let cpu_omega = omega_from_host(&registry, n);
        let cpu_pos = to_gpu(&atoms.pos, n);
        let cpu_vel = to_gpu(&atoms.vel, n);
        let state = resident.state.as_mut().expect("rows unchanged, so built");
        let changed = state.pos() != cpu_pos.as_slice()
            || state.vel() != cpu_vel.as_slice()
            || state.aux_state(omega) != cpu_omega.as_slice();
        if changed {
            state.write_pos(&cpu_pos);
            state.write_vel(&cpu_vel);
            state.set_aux_state(omega, &cpu_omega);
            dem.primed = false;
        }
        // The positions in hand are current, so a grid over them is too.
        if !(0..3).all(|d| domain.is_periodic(d)) {
            let cpu_dem = registry.expect::<DemAtom>("gpu_resident_step");
            let r_max = cpu_dem.radius[..n]
                .iter()
                .copied()
                .fold(f64::MIN_POSITIVE, f64::max) as GpuReal;
            state.set_grid(grid_for(&domain, &cpu_pos, 2.0 * r_max, max_cells));
        }
    }

    let window = dem.window;
    let state = resident.state.as_mut().expect("built above");
    if dem.primed {
        state.run_steps_continue(window);
    } else {
        state.run_steps(window);
    }
    dem.primed = true;
    ResidentAtoms::end_kernels(&mut coherence);
}

/// The whole DEM step on the GPU. Add it instead of the Verlet, rotation and contact
/// plugins; [`GranularResidentPlugins`] does that. See the [module docs](self).
pub struct GpuResidentPlugin {
    /// Steps on the device per schedule tick. A larger window means CPU systems and their
    /// transfers run less often.
    pub window: usize,
}

impl Plugin for GpuResidentPlugin {
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
        app.add_plugins(ResidentAtomsPlugin);
        let ranks = app
            .get_resource_ref::<CommResource>()
            .map_or(1, |comm| comm.0.size());
        assert_eq!(
            ranks, 1,
            "GpuResidentPlugin runs on one rank. Use GpuContactPlugin under MPI."
        );
        let gpu = app
            .get_resource_ref::<GpuDevice>()
            .expect("installed above")
            .clone();
        app.add_resource(ResidentDem {
            gpu,
            window: self.window.max(1),
            omega: 0,
            granular: None,
            walls: None,
            row_tags: Vec::new(),
            primed: false,
        });
        app.add_update_system(
            gpu_resident_step.label(CONTACT_FORCE),
            ParticleSimScheduleSet::Force,
        );
    }
}

/// `GranularDefaultPlugins` with the whole step on the GPU.
pub struct GranularResidentPlugins {
    /// See [`GpuResidentPlugin::window`].
    pub window: usize,
}

impl PluginGroup for GranularResidentPlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(dirt_atom::DemAtomPlugin)
            .add(dirt_atom::DemAtomInsertPlugin)
            .add(GpuResidentPlugin {
                window: self.window,
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_follow_their_tags() {
        let values = [10, 11, 20, 21, 30, 31];
        let out = remap_rows(&values, 2, &[1, 2, 3], &[3, 1, 7], -1);
        assert_eq!(out, [30, 31, 10, 11, -1, -1]);
    }
}
