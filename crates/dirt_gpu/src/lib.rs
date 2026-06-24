//! GPU DEM kernels for DIRT via `wgpu` (Metal on macOS, Vulkan/DX12 elsewhere).
//!
//! DIRT's discrete-element physics is expressed as **Force hooks** over soil's
//! resident GPU loop ([`soil_gpu::GpuState`]):
//!   - [`GranularForce`] — Hertz-normal + Mindlin-tangential contacts with
//!     persistent per-contact spring history and rotational torque,
//!   - [`WallForce`] — sphere-plane Hertz/Mindlin response against planar walls.
//!
//! All generic GPU substrate (device/context, cell-list neighbor build, the
//! resident velocity-Verlet loop + auxiliary-DOF rotation integration, the
//! ping-pong neighbor-state slots, planar boundary geometry) lives in `soil_gpu`
//! and is reused here — dirt holds only DEM-specific constitutive physics.
//!
//! GPU kernels are intrinsically f32 (Apple GPUs have no f64); host data is cast
//! to f32 on upload regardless of the build's [`soil_core::Real`] precision.

mod granular_force;
pub use granular_force::{GranularConfig, GranularForce};

mod wall_force;
pub use wall_force::{WallForce, MAX_WALLS};

mod bond_force;
pub use bond_force::{BondConfig, BondForce, BondTopology};

mod beam_bond_force;
pub use beam_bond_force::{BeamBondConfig, BeamBondForce};

// Generic GPU substrate owned by soil_gpu, re-exported for dirt_gpu consumers so
// they can build the whole hook stack from a single crate.
pub use soil_gpu::{Boundary, GpuContext, GpuState, Grid, Plane};
