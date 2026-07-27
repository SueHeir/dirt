//! Single-body sphero-cylinder (capsule) rigid particles for non-spherical DEM.
//!
//! A **capsule** is *one* rigid particle: a cylinder of radius `r` and length `L`
//! capped by two hemispheres of the same radius, carried by a full rigid-body
//! state (center of mass, orientation quaternion, angular velocity) with an
//! anisotropic inertia tensor. This is distinct from a `dirt_clump` "rod" built
//! from many bonded/overlapping sub-spheres: a capsule is a single atom with a
//! single contact geometry.
//!
//! A rod or fiber is assembled by **bonding capsules center-of-mass to
//! center-of-mass** with the existing `dirt_bond` beam model — no new bond
//! machinery is needed, because bonds already act between two atom centers.
//!
//! # Design
//! - The capsule's body symmetry axis is the body-frame **z**-axis; the world
//!   axis is `q · ẑ`. Its inertia tensor is diagonal in the body frame, so the
//!   principal axes are the identity rotation.
//! - A capsule atom sets `DemAtom::inv_inertia = 0`, so the scalar per-atom
//!   rotation integrator in `dirt_granular` skips it; [`CapsulePlugin`] advances
//!   its orientation with the full-tensor angular-momentum integrator in
//!   [`rigid`]. Translation is handled by the ordinary velocity-Verlet integrator.
//! - The world axis and half-length are `#[forward]` so ghost copies (periodic
//!   images and MPI halos) carry a fresh contact geometry every step.
//!
//! # Modules
//! - [`geometry`] — closest-point / overlap geometry for capsule contact.
//! - [`inertia`] — analytic mass and principal moments of a homogeneous capsule.
//! - [`rigid`]    — anisotropic free-atom rotational integration.

#![deny(missing_docs)]

use std::f64::consts::TAU;

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Deserialize;
use soil_derive::AtomData;

use soil_core::{
    register_atom_data, Atom, AtomData, AtomDataRegistry, CommResource, Config, Domain,
    ParticleSimScheduleSet, ParticleStore, ParticlesWith, Read, Real, Region, ScheduleSetupSet,
    Write,
};

use dirt_atom::{DemAtom, MaterialTable};

pub mod contact;
pub mod geometry;
pub mod inertia;
pub mod rigid;

pub use contact::{CapsuleContactHistory, CapsuleContactPlugin};
pub use geometry::{capsule_overlap, closest_segment_point, closest_segment_segment};
pub use inertia::capsule_mass_and_moments;
pub use rigid::{angmom_to_omega, richardson, quat_rotate, RigidRotState};

// ── Per-atom capsule data ────────────────────────────────────────────────────

/// Per-atom capsule geometry and rigid-body inertia.
///
/// Non-capsule atoms carry `is_capsule = 0` and `half_length = 0` (they are left
/// untouched by the capsule systems). `half_length` and `axis` are `#[forward]`
/// so ghost copies always present a current contact geometry; `principal_moments`
/// and `principal_axes` are only used by the owning rank's integrator.
#[derive(AtomData)]
pub struct CapsuleAtom {
    /// `1.0` for a capsule particle, `0.0` otherwise.
    #[forward]
    pub is_capsule: Vec<f64>,
    /// Half the cylindrical length `L/2` (m); the hemispherical caps add `radius`
    /// beyond each endpoint. `0.0` for non-capsules.
    #[forward]
    pub half_length: Vec<f64>,
    /// Current world-frame unit symmetry axis (`q · ẑ`), refreshed each step.
    #[forward]
    pub axis: Vec<[f64; 3]>,
    /// Principal moments of inertia `[I_x, I_y, I_z]` (kg·m²) about the COM.
    pub principal_moments: Vec<[f64; 3]>,
    /// Body→principal rotation quaternion (identity for an axis-aligned capsule).
    pub principal_axes: Vec<[f64; 4]>,
}

impl Default for CapsuleAtom {
    fn default() -> Self {
        Self::new()
    }
}

impl CapsuleAtom {
    /// Create an empty `CapsuleAtom` column set.
    pub fn new() -> Self {
        CapsuleAtom {
            is_capsule: Vec::new(),
            half_length: Vec::new(),
            axis: Vec::new(),
            principal_moments: Vec::new(),
            principal_axes: Vec::new(),
        }
    }
}

/// Returns `true` if atom `i` is a capsule particle.
#[inline]
pub fn is_capsule(caps: &CapsuleAtom, i: usize) -> bool {
    i < caps.is_capsule.len() && caps.is_capsule[i] != 0.0
}

/// World-frame centerline endpoints of capsule `i`: `(COM − h·axis, COM + h·axis)`.
#[inline]
pub fn capsule_endpoints(pos: [f64; 3], axis: [f64; 3], half_length: f64) -> ([f64; 3], [f64; 3]) {
    (
        [
            pos[0] - half_length * axis[0],
            pos[1] - half_length * axis[1],
            pos[2] - half_length * axis[2],
        ],
        [
            pos[0] + half_length * axis[0],
            pos[1] + half_length * axis[1],
            pos[2] + half_length * axis[2],
        ],
    )
}

// ── Configuration ────────────────────────────────────────────────────────────

/// A named capsule shape (`[[capsule.definitions]]`).
#[derive(Deserialize, Clone, Debug)]
#[serde(deny_unknown_fields)]
pub struct CapsuleDef {
    /// Shape name, referenced by insertion.
    pub name: String,
    /// Sphere/cap radius `r` (m).
    pub radius: f64,
    /// Cylindrical length `L` (m); the tip-to-tip length is `L + 2r`.
    pub length: f64,
}

/// A capsule insertion block (`[[capsule.insert]]`).
#[derive(Deserialize, Clone, Debug)]
pub struct CapsuleInsertConfig {
    /// Name of the capsule definition to insert.
    pub definition: String,
    /// Number of capsules to insert.
    pub count: u32,
    /// Density (kg/m³).
    pub density: f64,
    /// Material name (must match a `[[dem.materials]]` entry).
    pub material: String,
    /// Random per-component velocity magnitude (m/s); each component uniform in
    /// `[-v, v]`.
    #[serde(default)]
    pub velocity: Option<f64>,
    /// Insertion region; defaults to the domain inset by the capsule half-extent.
    #[serde(default)]
    pub region: Option<Region>,
    /// Give each capsule an independent uniformly-random 3-D orientation.
    #[serde(default)]
    pub random_orientation: bool,
    /// Deterministic insertion seed (default 0).
    #[serde(default)]
    pub seed: Option<u64>,
}

/// TOML `[capsule]` top-level configuration.
#[derive(Deserialize, Clone, Default)]
pub struct CapsuleTopConfig {
    /// Capsule shape definitions.
    #[serde(default)]
    pub definitions: Option<Vec<CapsuleDef>>,
    /// Capsule insertion commands.
    #[serde(default)]
    pub insert: Option<Vec<CapsuleInsertConfig>>,
}

/// Runtime store of capsule definitions, looked up during insertion.
#[derive(Default)]
pub struct CapsuleRegistry {
    /// All registered capsule definitions.
    pub defs: Vec<CapsuleDef>,
}

impl CapsuleRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        CapsuleRegistry { defs: Vec::new() }
    }
    /// Find a definition by name.
    pub fn find(&self, name: &str) -> Option<&CapsuleDef> {
        self.defs.iter().find(|d| d.name == name)
    }
}

// ── Plugin ───────────────────────────────────────────────────────────────────

/// Plugin adding single-body capsule (sphero-cylinder) particles.
pub struct CapsulePlugin;

impl Plugin for CapsulePlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn build(&self, app: &mut App) {
        register_atom_data!(app, CapsuleAtom::new());

        let mut registry = CapsuleRegistry::new();
        let config = Config::load::<CapsuleTopConfig>(app, "capsule");
        if let Some(defs) = config.definitions.clone() {
            for def in defs {
                assert!(def.radius > 0.0, "capsule '{}' radius must be > 0", def.name);
                assert!(def.length >= 0.0, "capsule '{}' length must be >= 0", def.name);
                registry.defs.push(def);
            }
        }
        app.add_resource(registry);
        app.add_resource(config);

        app.add_setup_system(capsule_insert_atoms, ScheduleSetupSet::Setup);
        app.add_setup_system(extend_ghost_cutoff_for_capsules, ScheduleSetupSet::Setup);

        // Advance orientation (full-tensor) alongside the translational integrator.
        app.add_update_system(
            capsule_initial_rotation,
            ParticleSimScheduleSet::InitialIntegration,
        );
        app.add_update_system(
            capsule_final_rotation,
            ParticleSimScheduleSet::FinalIntegration,
        );
        // Refresh the forwarded world axis after the quaternion moves, before the
        // ghost forward-comm at the next PreNeighbor.
        app.add_update_system(
            refresh_capsule_axis,
            ParticleSimScheduleSet::PostInitialIntegration,
        );
    }
}

/// Extend `Domain::ghost_cutoff` so a capsule near a subdomain/periodic boundary
/// is visible along its whole length. Mirrors the clump ghost-cutoff extension.
fn extend_ghost_cutoff_for_capsules(
    registry: Res<CapsuleRegistry>,
    mut domain: ResMut<Domain>,
) {
    let mut max_reach = 0.0_f64;
    for def in &registry.defs {
        max_reach = max_reach.max(0.5 * def.length + def.radius);
    }
    if max_reach > 0.0 {
        domain.ghost_cutoff += 2.0 * max_reach;
    }
}

/// Initial half-step orientation advance for every local capsule atom.
fn capsule_initial_rotation(
    atoms: Res<Atom>,
    particles: ParticlesWith<'_, (Write<DemAtom>, Read<CapsuleAtom>)>,
) {
    let dt = atoms.dt;
    let nlocal = atoms.nlocal() as usize;
    particles.with(|(mut dem, caps)| {
        for i in 0..nlocal.min(caps.is_capsule.len()) {
            if caps.is_capsule[i] == 0.0 {
                continue;
            }
            let torque = dem.torque[i];
            let mut q = dem.quaternion[i];
            let mut omega = dem.omega[i];
            let mut angmom = dem.ang_mom[i];
            let mut st = RigidRotState {
                quaternion: &mut q,
                omega: &mut omega,
                angmom: &mut angmom,
                moments: caps.principal_moments[i],
                principal_axes: caps.principal_axes[i],
            };
            st.integrate_initial(torque, dt);
            dem.quaternion[i] = q;
            dem.omega[i] = omega;
            dem.ang_mom[i] = angmom;
        }
    });
}

/// Final half-step orientation advance (after new torques) for local capsules.
fn capsule_final_rotation(
    atoms: Res<Atom>,
    particles: ParticlesWith<'_, (Write<DemAtom>, Read<CapsuleAtom>)>,
) {
    let dt = atoms.dt;
    let nlocal = atoms.nlocal() as usize;
    particles.with(|(mut dem, caps)| {
        for i in 0..nlocal.min(caps.is_capsule.len()) {
            if caps.is_capsule[i] == 0.0 {
                continue;
            }
            let torque = dem.torque[i];
            let mut q = dem.quaternion[i];
            let mut omega = dem.omega[i];
            let mut angmom = dem.ang_mom[i];
            let mut st = RigidRotState {
                quaternion: &mut q,
                omega: &mut omega,
                angmom: &mut angmom,
                moments: caps.principal_moments[i],
                principal_axes: caps.principal_axes[i],
            };
            st.integrate_final(torque, dt);
            dem.quaternion[i] = q;
            dem.omega[i] = omega;
            dem.ang_mom[i] = angmom;
        }
    });
}

/// Recompute each local capsule's world-frame axis from its quaternion so the
/// forwarded `axis` is current before the ghost forward-comm.
fn refresh_capsule_axis(
    atoms: Res<Atom>,
    particles: ParticlesWith<'_, (Read<DemAtom>, Write<CapsuleAtom>)>,
) {
    let nlocal = atoms.nlocal() as usize;
    particles.with(|(dem, mut caps)| {
        for i in 0..nlocal.min(caps.is_capsule.len()) {
            if caps.is_capsule[i] == 0.0 {
                continue;
            }
            caps.axis[i] = quat_rotate(dem.quaternion[i], [0.0, 0.0, 1.0]);
        }
    });
}

// ── Insertion ────────────────────────────────────────────────────────────────

/// Setup system: insert capsules from `[[capsule.insert]]`.
fn capsule_insert_atoms(
    comm: Res<CommResource>,
    domain: Res<Domain>,
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    capsule_registry: Res<CapsuleRegistry>,
    config: Res<CapsuleTopConfig>,
    material_table: Res<MaterialTable>,
    scheduler_manager: Res<SchedulerManager>,
) {
    if scheduler_manager.index != 0 || comm.rank() != 0 {
        return;
    }
    let inserts = match config.insert {
        Some(ref v) => v,
        None => return,
    };

    for insert in inserts {
        let def = capsule_registry
            .find(&insert.definition)
            .unwrap_or_else(|| panic!("capsule definition '{}' not found", insert.definition));
        let mat_idx = material_table
            .names
            .iter()
            .position(|n| n == &insert.material)
            .unwrap_or_else(|| panic!("material '{}' not found", insert.material))
            as u32;

        let reach = 0.5 * def.length + def.radius; // COM-to-tip extent
        let region = insert.region.clone().unwrap_or_else(|| Region::Block {
            min: [
                domain.boundaries_low[0] + reach,
                domain.boundaries_low[1] + reach,
                domain.boundaries_low[2] + reach,
            ],
            max: [
                domain.boundaries_high[0] - reach,
                domain.boundaries_high[1] - reach,
                domain.boundaries_high[2] - reach,
            ],
        });

        let mut rng = StdRng::seed_from_u64(insert.seed.unwrap_or(0));
        let inserted = insert_capsules(
            &mut atoms, &registry, def, insert, mat_idx, reach, &region, &mut rng,
        );
        if inserted < insert.count {
            eprintln!(
                "WARNING: inserted only {}/{} '{}' capsules",
                inserted, insert.count, insert.definition
            );
        }
        println!(
            "CapsuleInsert: inserted {} '{}' capsules (r={:.4}mm L={:.4}mm rho={})",
            inserted,
            insert.definition,
            def.radius * 1e3,
            def.length * 1e3,
            insert.density
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn insert_capsules<R: Rng>(
    atoms: &mut Atom,
    registry: &AtomDataRegistry,
    def: &CapsuleDef,
    insert: &CapsuleInsertConfig,
    mat_idx: u32,
    reach: f64,
    region: &Region,
    rng: &mut R,
) -> u32 {
    let mut centers: Vec<[f64; 3]> = Vec::new();
    let mut inserted = 0u32;
    let mut attempts = 0u64;
    let max_attempts = insert.count as u64 * 1_000_000;
    let min_sep = 2.0 * reach * 1.05;

    while inserted < insert.count && attempts < max_attempts {
        attempts += 1;
        let pos = match region.random_point_inside(rng) {
            Ok(p) => p,
            Err(_) => continue,
        };
        // Bounding-sphere overlap rejection against prior capsules and any atoms.
        let mut overlaps = centers.iter().any(|c| {
            let d = [pos[0] - c[0], pos[1] - c[1], pos[2] - c[2]];
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < min_sep * min_sep
        });
        if !overlaps {
            for k in 0..atoms.len() {
                let dx = pos[0] - atoms.pos[k][0] as f64;
                let dy = pos[1] - atoms.pos[k][1] as f64;
                let dz = pos[2] - atoms.pos[k][2] as f64;
                let md = reach + atoms.cutoff_radius[k] as f64;
                if dx * dx + dy * dy + dz * dz < md * md {
                    overlaps = true;
                    break;
                }
            }
        }
        if overlaps {
            continue;
        }

        let quat = if insert.random_orientation {
            random_unit_quaternion(rng)
        } else {
            [1.0, 0.0, 0.0, 0.0]
        };
        let vel = match insert.velocity {
            Some(v) => [
                rng.random_range(-v..v),
                rng.random_range(-v..v),
                rng.random_range(-v..v),
            ],
            None => [0.0; 3],
        };
        insert_capsule(
            atoms, registry, def, pos, vel, quat, insert.density, mat_idx,
        );
        centers.push(pos);
        inserted += 1;
    }
    inserted
}

/// Shoemake uniform random unit quaternion.
fn random_unit_quaternion<R: Rng>(rng: &mut R) -> [f64; 4] {
    let u1 = rng.random_range(0.0..1.0f64);
    let u2 = rng.random_range(0.0..1.0f64);
    let u3 = rng.random_range(0.0..1.0f64);
    [
        (1.0 - u1).sqrt() * (TAU * u2).sin(),
        (1.0 - u1).sqrt() * (TAU * u2).cos(),
        u1.sqrt() * (TAU * u3).sin(),
        u1.sqrt() * (TAU * u3).cos(),
    ]
}

/// Insert a single capsule atom at `com_pos` with orientation `quat`.
#[allow(clippy::too_many_arguments)]
pub fn insert_capsule(
    atoms: &mut Atom,
    registry: &AtomDataRegistry,
    def: &CapsuleDef,
    com_pos: [f64; 3],
    com_vel: [f64; 3],
    quat: [f64; 4],
    density: f64,
    atom_type: u32,
) {
    let (mass, moments) = capsule_mass_and_moments(def.radius, def.length, density);

    let global = atoms.natoms() + 1;
    ParticleStore::new(atoms, registry)
        .push_default_local(global)
        .expect("capsule row must be insertable");
    let i = atoms.len() - 1;

    atoms.tag[i] = atoms.get_max_tag() + 1;
    atoms.atom_type[i] = atom_type;
    atoms.pos[i] = [com_pos[0] as Real, com_pos[1] as Real, com_pos[2] as Real];
    atoms.vel[i] = [com_vel[0] as Real, com_vel[1] as Real, com_vel[2] as Real];
    atoms.force[i] = [0.0; 3];
    atoms.mass[i] = mass as Real;
    atoms.inv_mass[i] = if mass > 0.0 { (1.0 / mass) as Real } else { 0.0 };
    // Neighbor reach must cover the capsule's whole length plus its radius.
    atoms.cutoff_radius[i] = (0.5 * def.length + def.radius) as Real;
    atoms.image[i] = [0, 0, 0];
    atoms.is_ghost[i] = false;

    let axis = quat_rotate(quat, [0.0, 0.0, 1.0]);

    {
        let mut dem = registry.expect_mut::<DemAtom>("insert_capsule");
        dem.radius[i] = def.radius;
        dem.density[i] = density;
        dem.inv_inertia[i] = 0.0; // full-tensor path: scalar integrator skips this atom
        dem.quaternion[i] = quat;
        dem.omega[i] = [0.0; 3];
        dem.ang_mom[i] = [0.0; 3];
        dem.torque[i] = [0.0; 3];
        dem.body_id[i] = 0.0;
    }
    {
        let mut caps = registry.expect_mut::<CapsuleAtom>("insert_capsule");
        caps.is_capsule[i] = 1.0;
        caps.half_length[i] = 0.5 * def.length;
        caps.axis[i] = axis;
        caps.principal_moments[i] = moments;
        caps.principal_axes[i] = [1.0, 0.0, 0.0, 0.0];
    }
}
