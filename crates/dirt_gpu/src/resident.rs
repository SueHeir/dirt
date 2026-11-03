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
//! It runs on one rank, for what the kernels cover, with gravity, periodic axes, and fixed
//! infinite planar walls. Anything else panics at the first step. Use
//! [`GpuContactPlugin`](crate::GpuContactPlugin) for those: under MPI, ghosts would go
//! stale inside a window, and other physics has no kernel here.
//!
//! Three things do not carry through. Springs are not written to restart files. Periodic
//! image counters are not updated, so unwrapped positions are wrong. Quaternions are not
//! integrated, so a run needing particle orientation should use the CPU rotation path.
//!
//! Bonds run on the device too, from `DemBondPlugin`'s setup: elastic beams with
//! `combined_stress` breakage or none. Their history lives on the device for the run, and a broken
//! bond stops excluding its pair from contact at once, as on the CPU.

use std::collections::HashMap;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_bond::breakage::BreakageConfig;
use dirt_bond::{BondConfig, BondHistoryStore};
use dirt_schedule::CONTACT_FORCE;
use dirt_wall::Walls;
use grass_app::prelude::*;
use grass_gpu::GpuDevice;
use grass_mpi::CommResource;
use grass_scheduler::prelude::*;
use grass_scheduler::CoherenceRegistry;
use soil_core::{Atom, AtomDataRegistry, BondStore, Domain, ParticleSimScheduleSet};
use soil_fixes::GravityConfig;
use soil_gpu::neighbor_slots::{MAX_CONTACTS, PAYLOAD, SENTINEL};
use soil_gpu::resident::to_gpu;
use soil_gpu::{
    AtomColumnSync, Boundary, GpuReal, GpuState, Grid, Plane, ResidentAtoms, ResidentAtomsPlugin,
};

use crate::bond_force::{BondColumns, BondForce, BondTopology, BOND_STATE};
use crate::contact::{material_pairs, pair_row, require_supported};
use crate::granular_force::{GranularColumns, GranularConfig, GranularForce};
use crate::wall_force::{WallColumns, WallForce, MAX_WALLS};

const CELLS_PER_ATOM: usize = 8;

/// The rest of the resident state, beside `soil_gpu`'s [`ResidentAtoms`].
struct ResidentDem {
    gpu: GpuDevice,
    bonds: Option<BondColumns>,
    /// `(tag, partner tag)` at each bond end, in column order.
    bond_ends: Vec<(u32, u32)>,
    breakable: bool,
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

/// The contact table for each atom type against each active plane, in `planes` order.
fn wall_table(walls: Option<&Walls>, mt: &MaterialTable) -> Vec<GpuReal> {
    let Some(walls) = walls else {
        return Vec::new();
    };
    let mut rows = Vec::new();
    for i in 0..mt.e_eff_ij.len() {
        for (plane, _) in walls.planes.iter().zip(&walls.active).filter(|(_, &on)| on) {
            rows.extend(pair_row(mt, i, plane.material_index));
        }
    }
    rows
}

/// Why the kernels cannot run these bonds, if so.
fn require_supported_bonds(config: &BondConfig) {
    assert!(
        config.plasticity.is_none(),
        "GpuResidentPlugin runs elastic bonds. This run has plasticity. Use GpuContactPlugin, whose bonds stay on the CPU."
    );
    assert!(
        matches!(
            config.breakage,
            None | Some(BreakageConfig::Unbreakable) | Some(BreakageConfig::CombinedStress { .. })
        ),
        "GpuResidentPlugin runs combined_stress breakage or none. Use GpuContactPlugin, whose bonds stay on the CPU."
    );
    assert!(
        config.bond_radius_ratio > 0.0,
        "GpuResidentPlugin needs bond_radius_ratio > 0."
    );
}

/// The bonds as flat columns, carrying history from `old` where a bond's ends match, or else
/// from the CPU store. Fixed values follow the CPU's `bond_force`.
fn bond_topology(
    atoms: &Atom,
    registry: &AtomDataRegistry,
    config: &BondConfig,
    n: usize,
    old: Option<&HashMap<(u32, u32), [GpuReal; BOND_STATE]>>,
) -> BondTopology {
    let bonds = registry.expect::<BondStore>("bond_topology");
    let history = registry.expect::<BondHistoryStore>("bond_topology");
    let dem = registry.expect::<DemAtom>("bond_topology");
    let row: HashMap<u32, usize> = atoms.tag[..n]
        .iter()
        .enumerate()
        .map(|(i, &t)| (t, i))
        .collect();
    // History lives with the bond's lower tag, in that end's frame.
    let mut host: HashMap<(u32, u32), ([f64; 6], [f64; 2])> = HashMap::new();
    for (i, list) in history.history.iter().enumerate().take(n) {
        for entry in list {
            let d = entry.delta_t;
            let t = entry.delta_theta;
            host.insert(
                (atoms.tag[i], entry.partner_tag),
                (
                    [d[0], d[1], d[2], t[0], t[1], t[2]],
                    [entry.thresholds[0], entry.thresholds[1]],
                ),
            );
        }
    }
    let unbreakable = matches!(config.breakage, None | Some(BreakageConfig::Unbreakable));

    let mut topology = BondTopology {
        start: vec![0],
        ..Default::default()
    };
    for i in 0..n {
        for bond in bonds.bonds.get(i).map_or(&[][..], Vec::as_slice) {
            let tag_i = atoms.tag[i];
            let Some(&j) = row.get(&bond.partner_tag) else {
                continue;
            };
            let (lower, upper, sign) = if tag_i < bond.partner_tag {
                (tag_i, bond.partner_tag, 1.0)
            } else {
                (bond.partner_tag, tag_i, -1.0)
            };
            let seeded = host.get(&(lower, upper));
            let geometry = bond_values(atoms, &dem, config, i, j, bond.r0);
            let limits = match (unbreakable, seeded) {
                (false, Some((_, thresholds))) => *thresholds,
                _ => [f64::INFINITY; 2],
            };
            topology.partner.push(j as u32);
            topology.table.extend(geometry.map(|x| x as GpuReal));
            topology.table.extend(limits.map(|x| x as GpuReal));
            let state = match old.and_then(|old| old.get(&(tag_i, bond.partner_tag))) {
                Some(state) => *state,
                None => {
                    let start = seeded.map_or([0.0; 6], |(state, _)| *state);
                    start.map(|x| (sign * x) as GpuReal)
                }
            };
            topology.state.extend(state);
        }
        topology.start.push(topology.partner.len() as u32);
    }
    topology
}

/// `[r0, k_n, k_t, k_tor, k_bend, gamma_n, gamma_t, gamma_tor, gamma_bend, area, r_b, iben, jpol]`
fn bond_values(
    atoms: &Atom,
    dem: &DemAtom,
    config: &BondConfig,
    i: usize,
    j: usize,
    r0: f64,
) -> [f64; 13] {
    let r_b = config.bond_radius_ratio * dem.radius[i].min(dem.radius[j]);
    let area = std::f64::consts::PI * r_b * r_b;
    let jpol = 0.5 * std::f64::consts::PI * r_b.powi(4);
    let iben = 0.5 * jpol;
    let k_n = config
        .youngs_modulus
        .map_or(config.normal_stiffness, |e| e * area / r0);
    let k_t = config
        .shear_modulus
        .map_or(config.shear_stiffness, |g| g * area / r0);
    let k_tor = config
        .shear_modulus
        .map_or(config.twist_stiffness, |g| g * jpol / r0);
    let k_bend = config
        .youngs_modulus
        .map_or(config.bending_stiffness, |e| e * iben / r0);
    let (m_i, m_j) = (atoms.mass[i] as f64, atoms.mass[j] as f64);
    let m_red = if m_i + m_j > 0.0 {
        m_i * m_j / (m_i + m_j)
    } else {
        0.0
    };
    let moi = |x: f64| if x > 0.0 { 1.0 / x } else { 0.0 };
    let (moi_i, moi_j) = (moi(dem.inv_inertia[i]), moi(dem.inv_inertia[j]));
    let moi_red = if moi_i + moi_j > 0.0 {
        moi_i * moi_j / (moi_i + moi_j)
    } else {
        0.0
    };
    let gamma = |raw: Option<f64>, beta: f64, inertia: f64, k: f64| {
        raw.unwrap_or_else(|| 2.0 * beta * (inertia * k.max(0.0)).sqrt())
    };
    [
        r0,
        k_n,
        k_t,
        k_tor,
        k_bend,
        gamma(config.normal_damping, config.beta_normal, m_red, k_n),
        gamma(config.shear_damping, config.beta_shear, m_red, k_t),
        gamma(config.twist_damping, config.beta_twist, moi_red, k_tor),
        gamma(config.bending_damping, config.beta_bending, moi_red, k_bend),
        area,
        r_b,
        iben,
        jpol,
    ]
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
    wall_rows: Vec<GpuReal>,
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
    let mut config = GranularConfig::new(dt);
    (config.lx, config.ly, config.lz) = (lengths[0], lengths[1], lengths[2]);
    let zeros = vec![0.0; n];
    let contact = GranularForce::new(&state, grid, omega, &zeros, &material_pairs(mt), config);
    let granular = contact.columns();
    state.add_force_hook(Box::new(contact));
    // Walls come after contact, since contact is what starts each torque at zero.
    let walls = (!boundary.is_empty()).then(|| {
        let hook = WallForce::new(&state, omega, &zeros, boundary, &wall_rows, dt);
        let columns = hook.columns();
        state.add_force_hook(Box::new(hook));
        columns
    });
    (state, omega, granular, walls)
}

/// `(tag, partner tag)` at each bond end.
fn bond_ends(topology: &BondTopology, tags: &[u32]) -> Vec<(u32, u32)> {
    let mut ends = Vec::with_capacity(topology.partner.len());
    for (i, range) in topology.start.windows(2).enumerate() {
        for &j in &topology.partner[range[0] as usize..range[1] as usize] {
            ends.push((tags[i], tags[j as usize]));
        }
    }
    ends
}

/// Takes bonds the kernels broke off the CPU's lists.
fn sync_broken_bonds(dem: &ResidentDem, atoms: &Atom, registry: &AtomDataRegistry, n: usize) {
    let Some(columns) = dem.bonds.as_ref() else {
        return;
    };
    let partners = columns.read_partners();
    let mut bonds = registry.expect_mut::<BondStore>("sync_broken_bonds");
    let mut history = registry.expect_mut::<BondHistoryStore>("sync_broken_bonds");
    let local = n.min(atoms.nlocal() as usize);
    let row: HashMap<u32, usize> = atoms.tag[..local]
        .iter()
        .enumerate()
        .map(|(i, &t)| (t, i))
        .collect();
    for (&(tag, partner), &now) in dem.bond_ends.iter().zip(&partners) {
        if now == SENTINEL {
            if let Some(&i) = row.get(&tag) {
                if bonds.bonds[i].iter().any(|b| b.partner_tag == partner) {
                    bonds.bonds[i].retain(|b| b.partner_tag != partner);
                    if let Some(list) = history.history.get_mut(i) {
                        list.retain(|h| h.partner_tag != partner);
                    }
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn gpu_resident_step(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    mt: Res<MaterialTable>,
    domain: Res<Domain>,
    walls: Option<Res<Walls>>,
    gravity: Option<Res<GravityConfig>>,
    bond_config: Option<Res<BondConfig>>,
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

        // Bonds broken in the last window must leave the CPU's list first, so the topology
        // built below does not bring them back.
        if dem.breakable && resident.state.is_some() {
            sync_broken_bonds(dem, &atoms, &registry, dem.row_tags.len());
        }
        let old_bonds = dem.bonds.as_ref().map(|columns| {
            let state = columns.read_state();
            let mut by_end = HashMap::new();
            for (end, values) in dem.bond_ends.iter().zip(state.chunks_exact(BOND_STATE)) {
                let mut copy = [0.0; BOND_STATE];
                copy.copy_from_slice(values);
                by_end.insert(*end, copy);
            }
            by_end
        });

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
            let wall_rows = wall_table(walls.as_deref(), &mt);
            let (state, omega, granular, wall_columns) = build(
                &dem.gpu, &atoms, &mt, &domain, &boundary, wall_rows, g, &grid, n,
            );
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
        granular.write_type(&atoms.atom_type[..n]);
        if let Some(wall) = &dem.walls {
            wall.write_radius(&radius);
            wall.write_type(&atoms.atom_type[..n]);
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
        if let Some(config) = bond_config.as_deref() {
            require_supported_bonds(config);
            let topology = bond_topology(&atoms, &registry, config, n, old_bonds.as_ref());
            dem.bond_ends = bond_ends(&topology, tags);
            dem.breakable = matches!(config.breakage, Some(BreakageConfig::CombinedStress { .. }));
            match &dem.bonds {
                Some(columns) => columns.write(&topology),
                None => {
                    let periodic = [0, 1, 2].map(|d| {
                        if domain.is_periodic(d) {
                            domain.size[d] as GpuReal
                        } else {
                            0.0
                        }
                    });
                    let params = [
                        atoms.dt as GpuReal,
                        periodic[0],
                        periodic[1],
                        periodic[2],
                        0.0,
                    ];
                    let state = resident.state.as_mut().expect("just built");
                    let hook = BondForce::new(state, dem.omega, &topology, params);
                    dem.bonds = Some(hook.columns());
                    state.add_force_hook(Box::new(hook));
                }
            }
            // Contact reads the same partner column the bond kernel marks when a bond breaks.
            let partner_tags: Vec<u32> = dem.bond_ends.iter().map(|end| end.1).collect();
            let alive = dem.bonds.as_ref().expect("set above").alive();
            granular.write_bonds(&topology.start, &partner_tags, alive);
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

    // Bonds the kernels broke come off the CPU's list, and those two atoms can touch again.
    if dem.breakable && !rows_changed {
        sync_broken_bonds(dem, &atoms, &registry, n);
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
        // The CPU bond force would otherwise run each tick on a stale view.
        if app.has_resource::<BondConfig>() {
            app.remove_update_system(dirt_bond::bond_force);
        }
        app.add_resource(ResidentDem {
            gpu,
            bonds: None,
            bond_ends: Vec::new(),
            breakable: false,
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
    /// Adds `DemBondPlugin`, with its force running on the GPU.
    pub bonds: bool,
}

impl PluginGroup for GranularResidentPlugins {
    fn build(self) -> PluginGroupBuilder {
        let group = PluginGroupBuilder::start::<Self>()
            .add(dirt_atom::DemAtomPlugin)
            .add(dirt_atom::DemAtomInsertPlugin);
        let group = if self.bonds {
            group.add(dirt_bond::DemBondPlugin)
        } else {
            group
        };
        group.add(GpuResidentPlugin {
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
