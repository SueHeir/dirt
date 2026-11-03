//! GPU contact kernels for DIRT, on `grass_gpu` and `soil_gpu`.
//!
//! Hooks a resident [`GpuState`] runs:
//! - [`GranularForce`]: Hertz normal and Mindlin tangential contact with persistent
//!   spring history and torque.
//! - [`WallForce`]: sphere-plane contact with the same law.
//!
//! And two ways to use them in an ordinary DIRT App:
//! - [`GpuResidentPlugin`]: the whole step on the GPU for many steps at a time. Fastest,
//!   on one rank, for runs whose physics the kernels fully cover.
//! - [`GpuContactPlugin`]: only the contact force on the GPU. Everything else stays on the
//!   CPU, so it composes with MPI, CPU walls, fixes, sorting, migration and restart.
//!
//! Kernels are CubeCL `#[cube]` functions over `GpuReal`: `f32` on the Mac and `f64` with the
//! `cuda` feature. Compare `f32` runs against a CPU run within the single-precision band.
//!
//! The kernels cover Hertz-Mindlin with constant-torque rolling and twisting friction, and any
//! number of materials. Anything else is refused at startup, not silently ignored.

mod contact;
pub use contact::{GpuContact, GpuContactPlugin, GranularGpuPlugins};

mod resident;
pub use resident::{GpuResidentPlugin, GranularResidentPlugins};

mod granular_force;
pub use granular_force::{
    GranularColumns, GranularConfig, GranularForce, MaterialPairs, PAIR_FIELDS,
};

mod wall_force;
pub use wall_force::{WallColumns, WallForce, MAX_WALLS};

pub use soil_gpu::{Boundary, GpuDevice, GpuReal, GpuState, Grid, Plane};
