//! # DIRT — Discrete Element Method Simulation Framework
//!
//! DIRT is a plugin-based simulation framework for **Discrete Element Method
//! (DEM)** simulations of granular materials.
//!
//! This crate is the main entry point. It re-exports every plugin crate and
//! provides a convenience plugin group ([`CorePlugins`]) so that a complete
//! simulation can be set up in just a few lines.
//!
//! # Quick Start
//!
//! ## DEM granular simulation
//!
//! ```no_run
//! use dirt_core::prelude::*;
//!
//! fn main() {
//!     let mut app = App::new();
//!     app.add_plugins(CorePlugins)
//!        .add_plugins(GranularDefaultPlugins);
//!     app.start();
//! }
//! ```
//!
//! The simulation reads its parameters from a TOML config file passed as a
//! command-line argument (e.g. `cargo run -- config.toml`).
//!
//! # Architecture
//!
//! Every feature is implemented as a [`Plugin`](grass_app::Plugin) that
//! registers systems into a schedule. Plugins are grouped into **plugin groups**
//! for convenience:
//!
//! | Plugin Group | Purpose |
//! |---|---|
//! | [`CorePlugins`] | Infrastructure: input, communication, domain, neighbor lists, run loop, output |
//! | [`GranularDefaultPlugins`](dirt_granular::GranularDefaultPlugins) | DEM physics: atom data, insertion, Verlet integration, Hertz–Mindlin contact, rotational dynamics |
//!
//! ## Plugin-group membership
//!
//! The two groups are complementary — a typical run adds **both**. Membership:
//!
//! | Plugin | `CorePlugins` | `GranularDefaultPlugins` |
//! |---|:---:|:---:|
//! | `InputPlugin` (CLI + TOML config) | ✅ | |
//! | `CommunicationPlugin` (MPI / single-process) | ✅ | |
//! | `DomainPlugin` (decomposition, PBC, shrink-wrap) | ✅ | |
//! | `NeighborPlugin` (bin-based neighbor lists) | ✅ | |
//! | `GroupPlugin` (atom groups) | ✅ | |
//! | `RunPlugin` (run/cycle loop) | ✅ | |
//! | `PrintPlugin` (thermo, dump, restart) | ✅ | |
//! | `DemAtomPlugin` (radius, density, `MaterialTable`) | | ✅ |
//! | `DemAtomInsertPlugin` (`[[particles.insert]]`) | | ✅ |
//! | `VelocityVerletPlugin` (**translational integration**) | | ✅ |
//! | `HertzMindlinContactPlugin` (normal + tangential contact) | | ✅ |
//! | `RotationalDynamicsPlugin` (quaternion angular integration) | | ✅ |
//!
//! > ⚠️ **Neither group alone integrates motion.** Velocity Verlet lives in
//! > `GranularDefaultPlugins`, not `CorePlugins`. `CorePlugins` builds the
//! > simulation infrastructure but leaves atoms static; adding only
//! > `CorePlugins` gives a run that reads config, builds neighbor lists, and
//! > prints output but never moves a particle. For a non-granular method you
//! > would pair `CorePlugins` with `soil_verlet::VelocityVerletPlugin` (and
//! > your own force plugins) instead of `GranularDefaultPlugins`.
//! >
//! > `GranularTempPlugin` is **opt-in** — it is *not* in
//! > `GranularDefaultPlugins`; add it explicitly when you want
//! > granular-temperature output.
//!
//! For finer control, add individual plugins instead of a plugin group.
//!
//! # Crate Organization
//!
//! ## Framework (`grass_*`) and substrate (`soil_*`) crates
//!
//! | Crate | Description |
//! |---|---|
//! | [`grass_app`] | Application framework: [`App`](grass_app::App), [`Plugin`](grass_app::Plugin) trait, ECS-style resources, [`ScheduleSetupSet`](grass_app::ScheduleSetupSet) |
//! | [`grass_scheduler`] | System scheduler with [`ScheduleSet`](grass_scheduler::ScheduleSet) ordering, hierarchical schedules |
//! | [`grass_io`] | TOML config loading, multi-stage [`RunPlugin`](grass_io::RunPlugin), `SimClock` / `TermOut` / `Dump` plugins |
//! | [`grass_mpi`] | MPI abstraction (`CommBackend`, `SingleProcessComm`, `MpiCommBackend`) |
//! | [`grass_derive`] | Derive macros: `#[derive(ScheduleSet)]`, `#[derive(StageEnum)]`, `#[derive(Namespace)]` |
//! | [`soil_core`] | Domain decomposition, atom data, regions, groups; re-exports `Config` / `InputPlugin` / `RunPlugin` from `grass_io` |
//! | `neighbor` (in `soil_core`) | Bin-based neighbor list construction |
//! | [`soil_verlet`] | Velocity Verlet time integration |
//! | [`soil_print`] | Thermo output, dump files (CSV/binary/VTP), restart files |
//! | [`soil_derive`] | `#[derive(AtomData)]` proc macro |
//! | [`soil_deform`] | Box deformation: engineering strain rate, velocity, target size |
//!
//! ## DEM crates (`dirt_*`)
//!
//! | Crate | Description |
//! |---|---|
//! | [`dirt_atom`] | DEM per-atom data (`DemAtom`), material table, particle insertion, size distributions |
//! | [`dirt_fixes`] | General-purpose fixes: gravity, addforce, setforce, freeze (full immobilization), movelinear, viscous |
//! | [`soil_fixes`] | Method-agnostic translational position constraint: pin |
//! | [`dirt_granular`] | Hertz normal, Mindlin tangential, rotational dynamics, granular temperature |
//! | [`dirt_bond`] | Inter-particle bonds: normal/tangential/bending, auto-bonding, breakage |
//! | [`dirt_wall`] | Wall definitions: plane, cylinder, sphere, cone, region surfaces; wall motion |
//! | [`dirt_contact_analysis`] | Contact statistics: coordination number, fabric tensor, rattlers, per-contact CSV |
//! | [`dirt_measure_plane`] | Measurement planes for flux and profile sampling |
//!
//! # Feature Flags
//!
//! | Feature | Default | Description |
//! |---|---|---|
//! | `mpi_backend` | **yes** | Enable MPI-based parallel communication. Disable with `--no-default-features` for single-process runs. |
//!
//! # The Prelude
//!
//! The [`prelude`] module re-exports the most commonly used types so you can
//! write `use dirt_core::prelude::*;` and have everything you need for a typical
//! simulation. See the [`prelude`] module docs for the full list.

// ---------------------------------------------------------------------------
// Sub-crate re-exports
// ---------------------------------------------------------------------------

// --- DEM crates ---

/// DEM per-atom data, material properties, particle insertion, and size distributions.
pub use dirt_atom;

/// Inter-particle bonds for DEM: normal, tangential, bending modes with auto-bonding and breakage.
pub use dirt_bond;
/// Rigid clump (multisphere) composites for non-spherical DEM particles.
pub use dirt_clump;

/// Contact analysis: coordination number, fabric tensor, rattler detection, per-contact CSV output.
pub use dirt_contact_analysis;

/// Measurement planes for sampling particle flux and spatial profiles.
pub use dirt_measure_plane;

/// DEM granular contact models: Hertz normal, Mindlin tangential, rotational dynamics.
pub use dirt_granular;

/// Wall definitions (plane, cylinder, sphere, cone, region surfaces) and wall motion.
pub use dirt_wall;

// --- Infrastructure crates ---

/// Core simulation types: atoms, config, domain decomposition, communication, regions, groups.
pub use soil_core;

/// Box deformation: engineering strain rate, velocity, and target-size modes.
pub use soil_deform;

/// General-purpose fixes: gravity, addforce, setforce, freeze, movelinear, viscous.
pub use dirt_fixes;

/// Method-agnostic translational position constraint: pin.
pub use soil_fixes;

/// Thermo output, dump files (CSV/binary/VTP), and restart file I/O.
pub use soil_print;

/// Velocity Verlet time integration for translational degrees of freedom.
pub use soil_verlet;


use grass_app::prelude::*;

/// Core simulation infrastructure plugin group.
///
/// Includes, in registration order:
/// - [`InputPlugin`](soil_core::InputPlugin) — CLI parsing, banner printing, TOML config loading
///   (skipped if `Config` is already present)
/// - [`CommunicationPlugin`](soil_core::CommunicationPlugin) —
///   Unified MPI or single-process communication backend (selected by `mpi_backend` feature)
/// - [`DomainPlugin`](soil_core::DomainPlugin) — domain decomposition, PBC, and shrink-wrap
/// - [`NeighborPlugin`](soil_core::NeighborPlugin) — bin-based neighbor lists
/// - [`GroupPlugin`](soil_core::GroupPlugin) — atom group definitions and filtering
/// - [`RunPlugin`](soil_core::RunPlugin) — run/cycle management
/// - [`PrintPlugin`](soil_print::PrintPlugin) — thermo output, dump files, restart files
///
/// **Note:** Velocity Verlet integration is **not** included here. Use
/// [`VelocityVerletPlugin`](soil_verlet::VelocityVerletPlugin) directly.
/// The DEM plugin group [`GranularDefaultPlugins`](dirt_granular::GranularDefaultPlugins)
/// includes Velocity Verlet automatically.
///
/// MPI finalization is registered as a cleanup callback and runs automatically
/// at the end of [`App::start()`](grass_app::App::start).
///
/// # Usage
/// ```no_run
/// use dirt_core::prelude::*;
///
/// let mut app = App::new();
/// app.add_plugins(CorePlugins).add_plugins(GranularDefaultPlugins);
/// app.start();
/// ```
pub struct CorePlugins;

impl PluginGroup for CorePlugins {
    fn build(self) -> PluginGroupBuilder {
        let builder = PluginGroupBuilder::start::<Self>()
            .add(|app: &mut App| {
                app.set_warning_fn(soil_core::verlet_schedule_warnings);
            })
            .add(soil_core::InputPlugin);

        builder
            .add(soil_core::CommunicationPlugin)
            .add(soil_core::DomainPlugin)
            .add(soil_core::NeighborPlugin)
            .add(soil_core::GroupPlugin)
            .add(soil_core::RunPlugin)
            .add(soil_print::PrintPlugin)
    }
}

/// The DIRT prelude — import everything you need for a typical simulation.
///
/// ```no_run
/// use dirt_core::prelude::*;
/// # let _ = || { let mut app = App::new(); app.add_plugins(CorePlugins); };
/// ```
///
/// This re-exports:
///
/// ## Plugin groups
/// - [`CorePlugins`] — infrastructure (input, comm, domain, neighbor, run, print)
/// - [`GranularDefaultPlugins`] — Hertz–Mindlin + rotational dynamics + particle insertion
///
/// ## DEM types
/// - [`DemAtomPlugin`], [`DemConfig`], [`MaterialTable`] — DEM atom data and materials
/// - [`DemAtomInsertPlugin`], [`ParticlesConfig`] — particle insertion
/// - [`DemBondPlugin`] — inter-particle bonds
/// - [`WallPlugin`], [`Walls`], [`WallDef`], [`WallPlane`] — wall definitions
/// - [`ContactAnalysisPlugin`], [`ContactAnalysisConfig`] — contact statistics
/// - [`MeasurePlanePlugin`], [`MeasurePlanes`], [`MeasurePlaneDef`] — measurement planes
///
/// ## Shared types
/// - [`DeformPlugin`], [`DeformConfig`], [`DeformState`] — box deformation
/// - [`FixesPlugin`], [`GravityPlugin`] — general-purpose fixes
/// - [`StageEnum`] — derive macro for multi-stage runs
///
/// ## Core framework (via glob re-exports)
/// - [`App`](grass_app::App), [`Plugin`](grass_app::Plugin),
///   [`PluginGroup`](grass_app::PluginGroup) — application framework
/// - [`Atom`](soil_core::Atom), [`Config`](soil_core::Config),
///   [`RunState`](soil_core::RunState) — core simulation types
/// - [`Res`](grass_scheduler::Res), [`ResMut`](grass_scheduler::ResMut) — resource accessors
/// - [`ScheduleSet`](grass_scheduler::ScheduleSet) — system ordering
pub mod prelude {
    // Plugin groups defined in this crate
    pub use crate::CorePlugins;

    // DEM plugins and config types
    pub use dirt_atom::{DemAtomPlugin, DemConfig, MaterialTable};
    pub use dirt_atom::{DemAtomInsertPlugin, ParticlesConfig};
    pub use dirt_bond::DemBondPlugin;
    pub use dirt_clump::{ClumpPlugin, ClumpRegistry, ClumpAtom, ClumpDef, MultisphereBody, MultisphereBodyStore};
    pub use dirt_granular::{GranularDefaultPlugins, GranularGpuPlugins, GpuGranularResidentMpiPlugin, Boundary, Plane, GpuGranularForcePlugin, HertzMindlinContactPlugin, RotationalDynamicsPlugin, GranularTempPlugin};
    pub use dirt_contact_analysis::{ContactAnalysisConfig, ContactAnalysisPlugin};
    pub use dirt_measure_plane::{MeasurePlaneDef, MeasurePlanePlugin, MeasurePlanes};
    pub use dirt_wall::{WallDef, WallMotion, WallPlane, WallPlugin, Walls};

    // Shared infrastructure plugins
    pub use soil_deform::{DeformConfig, DeformPlugin, DeformState};
    pub use dirt_fixes::{AddForceDef, FixesPlugin, FixesRegistry, FreezeDef, GravityConfig, GravityPlugin, MoveLinearDef, SetForceDef, ViscousDef};
    pub use soil_fixes::{PinDef, PinRegistry, PinState, SoilFixesPlugin};

    // Derive macros
    pub use grass_derive::{ScheduleSet, StageEnum};
    pub use soil_derive::AtomData;

    // Core framework re-exports (glob)
    pub use grass_app::prelude::*;
    pub use soil_core::*;
    pub use soil_print::*;
    // Re-export the ParticleSimScheduleSet enum explicitly so downstream users
    // can access it without ambiguity with the ScheduleSet trait from grass_scheduler.
    pub use soil_core::ParticleSimScheduleSet;
    pub use grass_scheduler::prelude::*;
    pub use soil_verlet::*;
}
