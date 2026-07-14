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
//! For DEM-specific tests, prefer [`ParticleFixture`]:
//!
//! ```
//! use dirt_test_utils::ParticleFixture;
//!
//! let fixture = ParticleFixture::new().build();
//! assert_eq!(fixture.atom.nlocal, 2);
//! assert_eq!(fixture.neighbor.neighbor_indices, vec![1]);
//! ```
//!
//! # How to write a DIRT test
//!
//! [`ParticleFixture::build`] creates synchronized `Atom`/`DemAtom` rows, non-zero
//! counts, a stable timestep, material pair tables, and CSR neighbours. Add specialized
//! extension data with [`SimulationFixture::register_atom_data`], then either use
//! [`SimulationFixture::into_app`] or [`SimulationFixture::into_scheduled_app`].
//! `push_dem_test_atom` remains available for tests that deliberately exercise malformed
//! or partially assembled data.

use grass_app::App;
use grass_scheduler::{IntoScheduledSystem, ScheduleSet};
use soil_core::group::{Group, GroupDef};
use soil_core::{
    Atom, AtomData, AtomDataRegistry, CommResource, GroupRegistry, Neighbor, SingleProcessComm,
};

/// A conservative timestep used by [`ParticleFixture`] unless a smaller positive value is chosen.
pub const DEFAULT_DEM_TIMESTEP: f64 = 1e-7;

/// Largest timestep accepted by [`ParticleFixture::with_timestep`].
///
/// This is deliberately conservative for the stiff, small contact fixtures used by unit tests;
/// a test that needs a different integration regime should build its resources explicitly.
pub const MAX_TEST_DEM_TIMESTEP: f64 = 1e-6;

/// One particle declaration for a [`ParticleFixture`].
#[derive(Clone, Debug)]
pub struct ParticleSpec {
    /// Stable particle tag.
    pub tag: u32,
    /// Initial Cartesian position.
    pub position: [f64; 3],
    /// Sphere radius in metres.
    pub radius: f64,
}

impl ParticleSpec {
    /// Construct a particle declaration.
    pub fn new(tag: u32, position: [f64; 3], radius: f64) -> Self {
        assert!(
            radius.is_finite() && radius > 0.0,
            "fixture particle radius must be finite and positive"
        );
        Self {
            tag,
            position,
            radius,
        }
    }
}

/// Declarative builder for the common particle portion of a DIRT test.
///
/// [`Default`] is a two-particle overlapping pair, rather than an empty fixture. This makes a
/// forgotten atom insertion visible to force tests: the fixture always has non-zero counts and a
/// declared neighbour pair. Use [`single`](Self::single) for wall/fix tests, or
/// [`pair`](Self::pair) for a precisely positioned contact pair.
#[derive(Clone, Debug)]
pub struct ParticleFixture {
    particles: Vec<ParticleSpec>,
    pairs: Vec<(usize, usize)>,
    timestep: f64,
}

impl Default for ParticleFixture {
    fn default() -> Self {
        Self::pair(
            ParticleSpec::new(0, [0.0, 0.0, 0.0], 0.001),
            ParticleSpec::new(1, [0.0019, 0.0, 0.0], 0.001),
        )
    }
}

impl ParticleFixture {
    /// Start from the non-empty default contact fixture.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a valid one-particle fixture, useful for wall and fix tests.
    pub fn single(particle: ParticleSpec) -> Self {
        Self {
            particles: vec![particle],
            pairs: Vec::new(),
            timestep: DEFAULT_DEM_TIMESTEP,
        }
    }

    /// Build a valid two-particle fixture with one directed half-neighbour pair `0 -> 1`.
    pub fn pair(first: ParticleSpec, second: ParticleSpec) -> Self {
        Self {
            particles: vec![first, second],
            pairs: vec![(0, 1)],
            timestep: DEFAULT_DEM_TIMESTEP,
        }
    }

    /// Set a conservative positive DEM timestep.
    ///
    /// Panics for zero, negative, non-finite, or overly large values, preventing the usual
    /// one-step contact test from silently using an unstable timestep.
    pub fn with_timestep(mut self, timestep: f64) -> Self {
        assert!(
            timestep.is_finite() && timestep > 0.0 && timestep <= MAX_TEST_DEM_TIMESTEP,
            "fixture timestep must be finite, positive, and <= {MAX_TEST_DEM_TIMESTEP:e}"
        );
        self.timestep = timestep;
        self
    }

    /// Append a particle and return its local index.
    pub fn push_particle(&mut self, particle: ParticleSpec) -> usize {
        self.particles.push(particle);
        self.particles.len() - 1
    }

    /// Declare a directed CSR neighbour pair. Both indices must name fixture particles.
    pub fn add_pair(&mut self, pair: (usize, usize)) {
        assert!(
            pair.0 < self.particles.len() && pair.1 < self.particles.len(),
            "fixture pair index is outside the particle list"
        );
        assert_ne!(
            pair.0, pair.1,
            "fixture neighbour pair cannot be self-referential"
        );
        self.pairs.push(pair);
    }

    /// Materialize synchronized atom, DEM extension, material, neighbour, and registry resources.
    pub fn build(self) -> SimulationFixture {
        let mut atom = Atom::new();
        let mut dem = dirt_atom::DemAtom::new();
        for particle in &self.particles {
            push_dem_test_atom(
                &mut atom,
                &mut dem,
                particle.tag,
                particle.position,
                particle.radius,
            );
        }
        // The builder is never empty; keep this assertion next to the count assignment so future
        // changes cannot reintroduce the silent 0..nlocal force-loop footgun.
        assert!(
            !self.particles.is_empty(),
            "SimulationFixture requires at least one particle"
        );
        atom.nlocal = self.particles.len() as u32;
        atom.natoms = self.particles.len() as u64;
        atom.dt = self.timestep;

        let mut neighbor = Neighbor::new();
        let mut rows = vec![Vec::new(); self.particles.len()];
        for (i, j) in self.pairs {
            rows[i].push(j as u32);
        }
        neighbor.neighbor_offsets.push(0);
        for row in rows {
            neighbor.neighbor_indices.extend(row);
            neighbor
                .neighbor_offsets
                .push(neighbor.neighbor_indices.len() as u32);
        }

        let mut registry = AtomDataRegistry::new();
        registry
            .try_register(dem, atom.len())
            .expect("fixture DEM extension must match Atom length");
        SimulationFixture {
            atom,
            neighbor,
            registry,
            materials: make_material_table(),
        }
    }
}

/// Ready-to-use base resources for a single-process DIRT unit test.
///
/// The fixture owns `Atom`, its registered [`dirt_atom::DemAtom`] extension, a material table with
/// pair tables, and CSR neighbours. Domain-specific extensions such as bond histories or clump
/// data remain the responsibility of their owning crates and can be added with
/// [`register_atom_data`](Self::register_atom_data).
pub struct SimulationFixture {
    /// Synchronized particle arrays with non-zero `nlocal` and `natoms`.
    pub atom: Atom,
    /// CSR neighbours for the pairs declared on [`ParticleFixture`].
    pub neighbor: Neighbor,
    /// Atom-data registry containing a correctly sized `DemAtom`.
    pub registry: AtomDataRegistry,
    /// One valid glass material with built contact pair tables.
    pub materials: dirt_atom::MaterialTable,
}

impl SimulationFixture {
    /// Register another atom extension, checking it has exactly one row per fixture particle.
    pub fn register_atom_data<T: AtomData + 'static>(&mut self, data: T) {
        self.registry
            .try_register(data, self.atom.len())
            .expect("fixture atom extension must match Atom length");
    }

    /// Return the resources for a consumer that needs to add specialized resources before making an app.
    pub fn into_parts(self) -> (Atom, Neighbor, AtomDataRegistry, dirt_atom::MaterialTable) {
        (self.atom, self.neighbor, self.registry, self.materials)
    }

    /// Move the base resources into a new [`App`].
    pub fn into_app(self) -> App {
        let mut app = App::new();
        app.add_resource(self.atom);
        app.add_resource(self.neighbor);
        app.add_resource(self.registry);
        app.add_resource(self.materials);
        app
    }

    /// Move resources into an [`App`] and schedule one test system.
    ///
    /// The caller may add further resources or systems before calling `organize_systems`.
    pub fn into_scheduled_app<M>(
        self,
        system: impl IntoScheduledSystem<M>,
        set: impl ScheduleSet,
    ) -> App {
        let mut app = self.into_app();
        app.add_update_system(system, set);
        app
    }
}

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
        // `Atom` has no registered extensions in this core-only fixture.
        unsafe { atom.push_test_atom(i as u32, [i as f64, 0.0, 0.0], 0.5, 1.0) };
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
    // The extension rows below are appended in the same transaction.
    unsafe { atom.push_test_atom(tag, pos, radius, mass) };
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
    mt.add(
        dirt_atom::Material::new("glass", dirt_atom::Elastic::new(8.7e9, 0.3, 0.95)).with_friction(
            dirt_atom::Friction {
                sliding: 0.4,
                ..Default::default()
            },
        ),
    )
    .expect("test glass material is valid");
    mt.build_pair_tables();
    mt
}

#[cfg(test)]
mod fixture_tests {
    use super::*;
    use soil_core::ParticleSimScheduleSet;

    #[test]
    fn default_fixture_is_nonempty_and_synchronized() {
        let fixture = ParticleFixture::default().build();
        assert_eq!(fixture.atom.nlocal, 2);
        assert_eq!(fixture.atom.natoms, 2);
        assert_eq!(fixture.atom.len(), 2);
        assert_eq!(
            fixture
                .registry
                .expect::<dirt_atom::DemAtom>("fixture test")
                .len(),
            2
        );
        assert_eq!(fixture.neighbor.neighbor_offsets, vec![0, 1, 1]);
        assert_eq!(fixture.neighbor.neighbor_indices, vec![1]);
        assert_eq!(fixture.atom.dt, DEFAULT_DEM_TIMESTEP);
        assert_eq!(fixture.materials.e_eff_ij.len(), 1);
        assert_eq!(fixture.materials.e_eff_ij[0].len(), 1);
    }

    #[test]
    fn declared_pairs_become_csr_rows() {
        let mut builder = ParticleFixture::single(ParticleSpec::new(4, [0.0; 3], 0.002));
        let second = builder.push_particle(ParticleSpec::new(9, [1.0, 0.0, 0.0], 0.002));
        builder.add_pair((0, second));
        let fixture = builder.build();
        assert_eq!(fixture.neighbor.neighbor_offsets, vec![0, 1, 1]);
        assert_eq!(fixture.neighbor.neighbor_indices, vec![1]);
        assert_eq!(fixture.atom.tag, vec![4, 9]);
    }

    #[test]
    #[should_panic(expected = "fixture timestep")]
    fn rejects_unstable_timestep() {
        let _ = ParticleFixture::new().with_timestep(0.01);
    }

    fn mark_first_atom(mut atom: grass_scheduler::ResMut<Atom>) {
        atom.force[0][0] = 42.0;
    }

    #[test]
    fn can_schedule_a_system_without_manual_resource_wiring() {
        let mut app = ParticleFixture::new()
            .build()
            .into_scheduled_app(mark_first_atom, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();
        assert_eq!(app.get_resource_ref::<Atom>().unwrap().force[0][0], 42.0);
    }
}
