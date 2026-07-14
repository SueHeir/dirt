//! Wall contact forces for DEM simulations using Hertz or Hooke normal contact
//! with viscous damping, tangential/rolling/twisting friction, and optional
//! adhesion (JKR, DMT, SJKR cohesion).
//!
//! # Contact mechanics
//!
//! Wall contacts reuse the same per-pair mixing tables as particle–particle
//! contacts: the particle's material and the wall's `material` index a row of
//! [`MaterialTable`] (`e_eff_ij`, `g_eff_ij`, `beta_ij`, `friction_ij`,
//! `rolling_friction_ij`, `twisting_friction_ij`, …). Because a wall is treated
//! as infinitely massive and infinitely flat, the **effective radius is the
//! particle radius** (`R* = particle_radius`) and the reduced mass is the
//! particle mass.
//!
//! Beyond the normal force + damping, walls apply:
//!
//! - **Tangential (Mindlin sliding) friction** — incremental spring-history
//!   model with a per-contact tangential spring, Coulomb-capped at `μ |F_n|`.
//!   Supported by **all** wall types (plane, cylinder, sphere, region).
//! - **Rolling resistance** — `constant` (default) or `sds` (spring–dashpot–
//!   slider) model, mirroring the particle–particle rolling model with the wall
//!   as a zero-spin second body. Supported by **all** wall types.
//! - **Twisting friction** — constant-torque model opposing spin about the
//!   local contact normal. Supported by **all** wall types.
//!
//! Frictionless walls (`friction = 0`) are byte-for-byte unchanged from a
//! pure-normal contact. Tangential and rolling spring histories are stored on
//! the [`Walls`] resource (`tangential_springs`, `rolling_springs`).
//!
//! ## Adhesion-model asymmetry
//!
//! Adhesion support differs by wall geometry:
//!
//! - **Plane walls** support JKR and DMT (`surface_energy`) *and* SJKR cohesion
//!   (`cohesion_energy`), including the JKR extended-range pull-off regime.
//! - **Cylinder, sphere, and region walls** support **SJKR cohesion only**
//!   (`cohesion_energy`); their `surface_energy` is not consulted, so JKR/DMT
//!   pull-off is unavailable on curved/region walls. Use a plane wall if you
//!   need JKR/DMT against a wall.
//!
//! ## Wall temperature
//!
//! Every wall config accepts an optional `temperature` (K). This crate
//! **stores** it on the wall but never reads it — it is a hook for an external
//! heat-transfer system (e.g. a thermal-conduction plugin) to consult a wall's
//! temperature. It has no effect on the contact force computed here.
//!
//! # Wall Types
//!
//! | Type | Description | Config key |
//! |------|-------------|------------|
//! | **Plane** | Infinite flat plane defined by a point and unit normal | `type = "plane"` |
//! | **Cylinder** | Infinite cylinder along X/Y/Z axis with finite axial bounds | `type = "cylinder"` |
//! | **Sphere** | Sphere defined by center and radius | `type = "sphere"` |
//! | **Region** | Any [`Region`] shape used as a wall surface | `type = "region"` |
//!
//! All wall types treat the wall as having infinite mass and infinite radius
//! for contact mechanics, so the effective radius equals the particle radius
//! and the reduced mass equals the particle mass.
//!
//! # Motion Types
//!
//! | Motion | Description |
//! |--------|-------------|
//! | **Static** | Wall does not move (default) |
//! | **Constant velocity** | Wall translates at a fixed velocity each timestep |
//! | **Oscillating** | Sinusoidal displacement along the wall normal |
//! | **Servo** | Proportional controller adjusting velocity to reach a target force |
//!
//! Motion is currently supported only for plane walls.
//!
//! # TOML Configuration
//!
//! Walls are defined as `[[wall]]` array-of-tables entries. Each entry requires
//! a `material` field matching a name in `[[dem.materials]]`.
//!
//! ```toml
//! # Plane wall (floor at z=0, normal pointing up)
//! [[wall]]
//! type = "plane"
//! point_z = 0.0
//! normal_z = 1.0
//! material = "glass"
//! name = "floor"                  # optional, for runtime enable/disable
//!
//! # Cylinder wall (particles confined inside a z-aligned cylinder)
//! [[wall]]
//! type = "cylinder"
//! axis = "z"
//! center = [0.005, 0.005]         # center in the XY plane
//! radius = 0.004
//! lo = 0.0                        # axial lo bound (default: -inf)
//! hi = 0.01                       # axial hi bound (default: +inf)
//! inside = true                   # particles live inside the cylinder
//! material = "glass"
//!
//! # Sphere wall (particles confined inside a sphere)
//! [[wall]]
//! type = "sphere"
//! center = [0.005, 0.005, 0.005]
//! radius = 0.004
//! inside = true
//! material = "glass"
//!
//! # Region wall (any Region shape as a wall surface)
//! [[wall]]
//! type = "region"
//! inside = true
//! material = "glass"
//! region = { type = "cone", center = [0.005, 0.005], axis = "z",
//!            rad_lo = 0.004, rad_hi = 0.002, lo = 0.0, hi = 0.01 }
//!
//! # Moving wall with constant velocity
//! [[wall]]
//! type = "plane"
//! normal_z = 1.0
//! material = "glass"
//! velocity = [0.0, 0.0, -0.01]    # [vx, vy, vz]
//!
//! # Oscillating wall (sinusoidal along normal)
//! [[wall]]
//! type = "plane"
//! point_z = 0.1
//! normal_z = 1.0
//! material = "glass"
//! oscillate = { amplitude = 0.001, frequency = 50.0 }
//!
//! # Servo-controlled wall (adjusts velocity to reach target force)
//! [[wall]]
//! type = "plane"
//! point_z = 0.1
//! normal_z = -1.0
//! material = "glass"
//! servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }
//! ```
//!
//! # Plugin Registration
//!
//! Add [`WallPlugin`] to your app. It depends on `DemAtomPlugin` (for
//! [`MaterialTable`] and [`DemAtom`] data).
//!
//! Wall configuration is preflighted through [`Plugin::try_build`]. A malformed
//! `[[wall]]` entry is returned as an [`AppError`] to the outer runner, which
//! owns formatting and MPI-consistent termination.
//!
//! # Systems
//!
//! | System | Schedule | Purpose |
//! |--------|----------|---------|
//! | [`wall_move`] | `PreInitialIntegration` | Updates wall positions from motion modes |
//! | [`wall_zero_force_accumulators`] | `PreForce` | Zeros per-wall force accumulators |
//! | [`wall_contact_force`] | `Force` | Computes normal contact + damping + adhesion |

// Public API documentation-completeness gate: every public item in this crate
// must carry a doc comment. Enforced on both `cargo build` (rustc) and
// `cargo doc` (rustdoc; e.g. `RUSTDOCFLAGS="-D missing_docs"`). Document real
// API intent here — do not add empty doc comments just to satisfy the lint.
#![deny(missing_docs)]

mod config;
mod contact;
mod geometry;
mod motion;
mod plugin;

pub use config::{OscillateDef, ServoDef, WallDef};
pub use contact::wall_contact_force;
pub use geometry::{WallCylinder, WallPlane, WallRegion, WallSphere, Walls};
pub use motion::{wall_move, wall_zero_force_accumulators, WallMotion};
pub use plugin::WallPlugin;

#[cfg(test)]
mod tests;
