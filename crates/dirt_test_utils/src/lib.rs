//! Shared test utilities for all DIRT crates.
//!
//! This crate provides helper functions to quickly set up [`Atom`], [`GroupRegistry`],
//! [`CommResource`], [`DemAtom`](dirt_atom::DemAtom), and [`MaterialTable`](dirt_atom::MaterialTable)
//! instances for use in unit tests. Each helper produces a minimal, valid object so tests can
//! focus on the logic under test rather than boilerplate setup.
//!
//! # Quick Start
//!
//! ```
//! use dirt_test_utils::{make_atoms, make_group_registry, make_single_comm};
//!
//! let atom = make_atoms(3); // 3 atoms at (0,0,0), (1,0,0), (2,0,0)
//! let groups = make_group_registry("all", vec![true, true, true]);
//! let _comm = make_single_comm();
//! assert_eq!(atom.nlocal, 3);
//! assert_eq!(groups.groups[0].count, 3);
//! ```
//!
//! For DEM-specific tests, use [`push_dem_test_atom`] and [`make_material_table`]:
//!
//! ```
//! use dirt_test_utils::{push_dem_test_atom, make_material_table};
//!
//! let mut atom = soil_core::Atom::new();
//! let mut dem = dirt_atom::DemAtom::default();
//! push_dem_test_atom(&mut atom, &mut dem, 0, [0.0, 0.0, 0.0], 0.5);
//! // NOTE: push_dem_test_atom does NOT set nlocal/natoms — do it yourself:
//! atom.nlocal = 1;
//! atom.natoms = 1;
//!
//! let materials = make_material_table(); // single "glass" material
//! assert_eq!(materials.names.len(), 1);
//! ```
//!
//! # How to write a DIRT test
//!
//! These helpers cover atom/material *construction* only. A force-law test also
//! needs an `App` (from `grass_app`) with a neighbor list and the system under
//! test scheduled into it. That assembly lives in the *consumer* crate (which
//! depends on `grass_app`/`grass_scheduler`), not here — `dirt_test_utils`
//! deliberately depends only on `soil_core` + `dirt_atom`, so the snippet below
//! is illustrative (`text`, not a runnable doctest). It mirrors the real
//! pattern in `dirt_granular`'s `contact.rs` tests (e.g.
//! `fused_contact_repulsive_for_overlap`):
//!
//! ```text
//! use dirt_test_utils::{push_dem_test_atom, make_material_table};
//! use dirt_atom::DemAtom;
//! use soil_core::{Atom, AtomDataRegistry, Neighbor, ParticleSimScheduleSet};
//! use grass_app::prelude::*;          // App, add_update_system, organize_systems, run, get_resource_ref
//!
//! let mut app = App::new();
//! let radius = 0.001;
//!
//! // 1. Build Atom + DemAtom in parallel with the helper.
//! let mut atom = Atom::new();
//! let mut dem = DemAtom::new();
//! atom.dt = 1e-7;                     // DEM needs a tiny dt (~1e-7 s); the helper does NOT set it
//! push_dem_test_atom(&mut atom, &mut dem, 0, [0.0, 0.0, 0.0], radius);
//! push_dem_test_atom(&mut atom, &mut dem, 1, [0.0019, 0.0, 0.0], radius);
//!
//! // 2. Manually set the counts — push_dem_test_atom does NOT touch them.
//! atom.nlocal = 2;
//! atom.natoms = 2;
//!
//! // 3. Build a Neighbor list by hand (CSR offsets/indices for the pair 0–1).
//! let mut neighbor = Neighbor::new();
//! neighbor.neighbor_offsets = vec![0, 1, 1];
//! neighbor.neighbor_indices = vec![1];
//!
//! // 4. Register per-atom data (and any history store) into the registry.
//! let mut registry = AtomDataRegistry::new();
//! registry.register(dem);
//!
//! // 5. Add resources, schedule the system under test at the Force set, organize, run.
//! app.add_resource(atom);
//! app.add_resource(neighbor);
//! app.add_resource(registry);
//! app.add_resource(make_material_table());
//! app.add_update_system(my_force_system, ParticleSimScheduleSet::Force);
//! app.organize_systems();
//! app.run();
//!
//! // 6. Read results back off the App resources and assert.
//! let atom = app.get_resource_ref::<Atom>().unwrap();
//! assert!(atom.force[0][0] < 0.0);
//! ```
//!
//! ## Two footguns
//!
//! - **`nlocal` / `natoms` are NOT set by [`push_dem_test_atom`].** It only
//!   *appends* one atom to the `Atom` arrays and the parallel `DemAtom` arrays;
//!   it leaves `atom.nlocal` and `atom.natoms` untouched. Force systems iterate
//!   `0..nlocal`, so if you forget to set them (they default to 0) **your system
//!   silently does nothing** — no panic, just an empty loop and a passing-but-
//!   meaningless test. Set both after your last `push_*` call.
//! - **DEM needs a tiny timestep.** [`make_atoms`] sets `dt = 0.001`, which is
//!   far too large for stiff DEM contacts. For any contact/bond test, override
//!   `atom.dt` to roughly `1e-7` (the value the real `dirt_granular`/`dirt_bond`
//!   tests use) so a single `app.run()` step doesn't blow the contact up.
//!
//! # What this crate does NOT provide
//!
//! - **No `App` / `Neighbor` builders.** You construct and wire those in the
//!   consuming crate (see the assembly above). The neighbor list in particular
//!   is built by hand (raw CSR `neighbor_offsets` / `neighbor_indices`).
//! - **No custom assertions or test macros.** Use plain `assert!` /
//!   `assert_eq!`; these helpers only build inputs.
//! - **No scheduling / plugin helpers.** Schedule sets and systems are added
//!   directly via the `grass_app` / `grass_scheduler` API in the consumer.

use soil_core::group::{Group, GroupDef};
use soil_core::{Atom, CommResource, GroupRegistry, SingleProcessComm};

/// Create an [`Atom`] with `n` test atoms arranged along the x-axis.
///
/// Each atom `i` is placed at position `(i, 0, 0)` with radius `0.5`, mass `1.0`,
/// and a timestep of `0.001`. Atom tags are `0..n`.
///
/// # Examples
///
/// ```
/// use dirt_test_utils::make_atoms;
///
/// let atom = make_atoms(5);
/// assert_eq!(atom.nlocal, 5);
/// assert_eq!(atom.natoms, 5);
/// assert_eq!(atom.pos[2], [2.0, 0.0, 0.0]);
/// ```
pub fn make_atoms(n: usize) -> Atom {
    let mut atom = Atom::new();
    for i in 0..n {
        atom.push_test_atom(i as u32, [i as f64, 0.0, 0.0], 0.5, 1.0);
    }
    atom.nlocal = n as u32;
    atom.natoms = n as u64;
    atom.dt = 0.001;
    atom
}

/// Create a [`GroupRegistry`] containing a single named group with the given membership mask.
///
/// The `mask` vector should have one entry per atom: `true` means the atom belongs to the
/// group, `false` means it does not. The group count is computed automatically from the mask.
///
/// # Examples
///
/// ```
/// use dirt_test_utils::make_group_registry;
///
/// // Group "mobile" includes atoms 0 and 2, but not atom 1
/// let registry = make_group_registry("mobile", vec![true, false, true]);
/// assert_eq!(registry.groups[0].count, 2);
/// ```
pub fn make_group_registry(name: &str, mask: Vec<bool>) -> GroupRegistry {
    let count = mask.iter().filter(|&&m| m).count();
    let mut registry = GroupRegistry::new();
    registry.groups.push(Group {
        name: name.to_string(),
        def: GroupDef {
            name: name.to_string(),
            atom_types: None,
            region: None,
            dynamic: None,
        },
        mask,
        count,
    });
    registry
}

/// Create a single-process [`CommResource`] for testing.
///
/// Returns a communication resource backed by [`SingleProcessComm`], which requires no
/// MPI initialization. Suitable for all single-process unit tests.
///
/// # Examples
///
/// ```
/// use dirt_test_utils::make_single_comm;
///
/// let _comm = make_single_comm();
/// ```
pub fn make_single_comm() -> CommResource {
    CommResource(Box::new(SingleProcessComm::new()))
}

/// Push a DEM test atom with all [`DemAtom`](dirt_atom::DemAtom) fields populated.
///
/// Creates a solid sphere using a density of `2500 kg/m³`. Mass is computed from the
/// given `radius` via `m = ρ · (4/3)πr³`, and the moment of inertia assumes a uniform
/// solid sphere (`I = 0.4 · m · r²`). Rotational fields (quaternion, omega, angular
/// momentum, torque) are initialized to zero/identity defaults.
///
/// Both `atom` and `dem` are extended in parallel — callers must ensure they stay in sync.
///
/// # Examples
///
/// ```
/// use dirt_test_utils::push_dem_test_atom;
///
/// let mut atom = soil_core::Atom::new();
/// let mut dem = dirt_atom::DemAtom::default();
/// push_dem_test_atom(&mut atom, &mut dem, 0, [1.0, 2.0, 3.0], 0.5);
/// assert_eq!(dem.radius[0], 0.5);
/// assert_eq!(dem.density[0], 2500.0);
/// // push_dem_test_atom does NOT set nlocal/natoms — the caller must:
/// atom.nlocal = 1;
/// atom.natoms = 1;
/// ```
pub fn push_dem_test_atom(
    atom: &mut Atom,
    dem: &mut dirt_atom::DemAtom,
    tag: u32,
    pos: [f64; 3],
    radius: f64,
) {
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    atom.push_test_atom(tag, pos, radius, mass);
    dem.radius.push(radius);
    dem.density.push(density);
    dem.inv_inertia.push(1.0 / (0.4 * mass * radius * radius));
    dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
    dem.omega.push([0.0; 3]);
    dem.ang_mom.push([0.0; 3]);
    dem.torque.push([0.0; 3]);
    dem.body_id.push(0.0);
}

/// Create a single-material "glass" [`MaterialTable`](dirt_atom::MaterialTable) for testing.
///
/// Returns a material table with one material named `"glass"` that has the following
/// properties:
///
/// | Property          | Value   |
/// |-------------------|---------|
/// | Young's modulus   | 8.7 GPa |
/// | Poisson ratio     | 0.3     |
/// | Restitution       | 0.95    |
/// | Friction          | 0.4     |
/// | Rolling friction  | 0.0     |
/// | Cohesion energy   | 0.0     |
///
/// Pair tables are pre-built, so the table is ready for contact computations immediately.
///
/// # Examples
///
/// ```
/// use dirt_test_utils::make_material_table;
///
/// let mt = make_material_table();
/// assert_eq!(mt.names.len(), 1);
/// assert_eq!(mt.find_material("glass"), Some(0));
/// ```
pub fn make_material_table() -> dirt_atom::MaterialTable {
    let mut mt = dirt_atom::MaterialTable::new();
    mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0);
    mt.build_pair_tables();
    mt
}
