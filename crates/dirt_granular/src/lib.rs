//! DEM granular contact models for spherical particle simulations.
//!
//! This crate provides the core physics for Discrete Element Method (DEM) simulations
//! of granular materials. It implements contact force models, rotational dynamics,
//! and granular temperature output.
//!
//! # Contact models
//!
//! ## Normal contact
//! - **Hertz** (default) — nonlinear elastic contact: `F_n = 4/3 E* √(R* δ) · δ`
//!   with viscoelastic damping proportional to `β √(S_n m_r)`
//! - **Hooke** — linear spring contact: `F_n = k_n δ` with linear damping `γ_n v_n`
//!
//! ## Tangential contact
//! - **Mindlin** — incremental spring-history model with Coulomb friction cap `μ |F_n|`.
//!   Spring displacement is stored per-contact and rotated to stay in the tangent plane
//!   each step. Damping: `γ_t = 2 β √(5/6) √(k_t m_r)`.
//!
//! ## Rolling resistance
//! - **Constant torque** (default) — `τ_r = μ_r |F_n| R*` opposing relative rolling
//! - **SDS** (spring-dashpot-slider) — incremental rolling displacement with spring
//!   stiffness, viscous damping, and Coulomb-style slider cap
//!
//! ## Twisting friction
//! - **Constant torque** (default) — `τ_tw = μ_tw |F_n| R*` opposing relative twisting
//! - **SDS** (spring-dashpot-slider) — incremental twist angle with spring, damping, cap
//!
//! ## Adhesion / cohesion
//! - **JKR** — Johnson-Kendall-Roberts adhesion with extended interaction range beyond
//!   geometric contact; pull-off force `F = 3/2 π γ R*`
//! - **DMT** — Derjaguin-Muller-Toporov adhesion with constant attractive force
//!   `F = 2π γ R*` during contact only
//! - **SJKR** — simplified cohesion proportional to contact area: `F = k_coh π δ R*`
//!
//! > ⚠️ **Hooke/Hertz adhesion asymmetry.** JKR and DMT adhesion (driven by
//! > `surface_energy`) are implemented **only on the Hertz contact path**. Under
//! > `contact_model = "hooke"` the `surface_energy` term is *silently ignored* —
//! > the linear-spring path applies SJKR cohesion (`cohesion_energy`) only. If
//! > you need JKR/DMT pull-off, use the default Hertz model.
//!
//! # TOML configuration
//!
//! Contact model parameters are set per-material in the `[[dem.materials]]` array:
//!
//! ```toml
//! [[dem.materials]]
//! name = "glass"
//! youngs_modulus = 8.7e9      # Pa — Young's modulus E
//! poisson_ratio = 0.3         # dimensionless — Poisson's ratio ν
//! restitution = 0.95          # dimensionless — coefficient of restitution (0–1)
//! friction = 0.4              # dimensionless — sliding friction coefficient μ
//! rolling_friction = 0.1      # dimensionless — rolling friction coefficient μ_r
//! cohesion_energy = 0.0       # J/m² — SJKR cohesion energy density (0 = disabled)
//! surface_energy = 0.0        # J/m² — JKR/DMT surface energy γ (0 = disabled)
//! ```
//!
//! Global model selection:
//!
//! ```toml
//! [dem]
//! contact_model = "hertz"     # "hertz" (default) or "hooke"
//! adhesion_model = "jkr"      # "jkr" (default) or "dmt" (only when surface_energy > 0)
//! rolling_model = "constant"  # "constant" (default) or "sds"
//! twisting_model = "constant" # "constant" (default) or "sds"
//! ```
//!
//! # Material-parameter reference
//!
//! Every parameter above is stored per-material in [`dirt_atom::MaterialTable`]
//! and mixed into per-pair tables by `MaterialTable::build_pair_tables()`. Which
//! `MaterialTable` fields each model branch reads:
//!
//! | Model branch | `MaterialTable` inputs (per-pair table) |
//! |---|---|
//! | Hertz normal | `e_eff_ij` (E*), `beta_ij` (from `restitution`) |
//! | Hooke normal | `kn_ij` (harmonic mean of per-material `kn`), `beta_ij` |
//! | Mindlin tangential | `g_eff_ij` (G*), `friction_ij` (μ), `beta_ij` |
//! | Hooke tangential | `kt_ij`, `friction_ij` |
//! | Rolling (constant / SDS) | `rolling_friction_ij`, `rolling_stiffness_ij`, `rolling_damping_ij` |
//! | Twisting (constant / SDS) | `twisting_friction_ij`, `twisting_stiffness_ij`, `twisting_damping_ij` |
//! | JKR / DMT adhesion (Hertz only) | `surface_energy_ij` (γ), `adhesion_model` |
//! | SJKR cohesion | `cohesion_energy_ij` |
//!
//! `restitution` is the **target COR**: `beta_ij` is derived by inverting the
//! exact head-on Hertz collision (see [`dirt_atom::hertz_beta_for_cor`]).
//!
//! # Tangential / rolling / twisting history (canonical frame)
//!
//! The Mindlin tangential force and the SDS rolling/twisting variants are
//! **incremental, history-dependent** springs: a displacement is integrated
//! across timesteps, rotated to stay in the current tangent plane, and capped at
//! a Coulomb limit. That history lives in
//! [`tangential::ContactHistoryStore`], which stores **7 `f64` per contact** —
//! `[0..3]` tangential spring vector, `[3..6]` rolling spring vector, `[6]`
//! twisting scalar (rolling/twisting slots are zero under the constant-torque
//! models).
//!
//! Each entry is kept in **canonical form** (from the lower-tag particle's
//! perspective) so the spring is frame-consistent no matter which particle is
//! `i` vs `j` in the neighbor list; a `sign` factor of `±1` flips the canonical
//! spring into the local `(i, j)` frame each step. The Coulomb friction limit is
//! applied in **two stages**: the stored spring is first capped so that
//! `|k_t s| ≤ μ|F_n|` (truncating the history that survives to the next step),
//! then the assembled force `F_t = k_t s − γ_t v_t` is capped again at `μ|F_n|`.
//!
//! # Modules
//!
//! - [`contact`] — Fused Hertz-Mindlin + Hooke contact force (primary code path)
//! - [`tangential`] — Per-contact tangential spring-history store (`ContactHistoryStore`)
//! - [`rotational`] — Quaternion-based velocity Verlet for angular degrees of freedom
//! - [`granular_temp`] — Granular temperature (velocity fluctuation) output;
//!   the plugin ([`GranularTempPlugin`]) is **opt-in** and is *not* part of
//!   [`GranularDefaultPlugins`] — add it explicitly when you want
//!   `data/GranularTemp.txt` written.

pub mod granular_temp;
pub mod rotational;
pub mod tangential;

pub use granular_temp::GranularTempPlugin;
pub use rotational::RotationalDynamicsPlugin;

pub mod contact;

use grass_app::prelude::*;

use dirt_atom::DemAtomPlugin;
use dirt_atom::DemAtomInsertPlugin;
use soil_verlet::VelocityVerletPlugin;

pub use contact::HertzMindlinContactPlugin;

/// Re-export from [`dirt_atom`] for convenience.
pub use dirt_atom::SQRT_5_6;
/// Small epsilon to avoid division by zero when normalizing tangential,
/// rolling, or twisting spring displacements.
pub const TANGENTIAL_EPSILON: f64 = 1e-30;

/// Warn when `distance / (r1 + r2)` falls below this threshold.
///
/// A ratio near 0.0 means nearly full overlap, which indicates an unstable
/// simulation (timestep too large or bad initial packing). Contacts with
/// overlap exceeding this threshold trigger a warning but still compute
/// forces (capped at half the smaller radius) to prevent runaway penetration.
pub const LARGE_OVERLAP_WARN_THRESHOLD: f64 = 0.90;

/// Maximum overlap warnings per timestep before the simulation panics.
///
/// If more than this many pairs exceed [`LARGE_OVERLAP_WARN_THRESHOLD`],
/// the simulation aborts with an actionable error message suggesting the
/// user reduce the timestep or fix the initial configuration.
pub const MAX_OVERLAP_WARNINGS: usize = 500;

/// Default DEM granular physics plugin group.
///
/// Includes, in registration order:
/// - [`DemAtomPlugin`] — per-atom material properties (radius, density) and
///   `MaterialTable` for per-material Young's modulus, Poisson ratio, restitution,
///   friction with geometric-mean mixing
/// - [`DemAtomInsertPlugin`] — random particle insertion from `[[particles.insert]]` config
/// - [`RotationalDynamicsPlugin`] — quaternion Velocity Verlet for angular degrees of freedom
///   (I = 2/5 m r² for solid spheres)
///
/// Granular temperature output ([`GranularTempPlugin`], which writes
/// `data/GranularTemp.txt`) is **not** bundled — add it explicitly only in the
/// examples that consume that file, so ordinary runs don't emit it.
///
/// Does **not** include infrastructure plugins (input, comm, domain, neighbor,
/// run, print). Use [`CorePlugins`] to get all infrastructure.
///
/// # Usage
/// ```rust,ignore
/// use dirt::prelude::*;
///
/// let mut app = App::new();
/// app.add_plugins(CorePlugins).add_plugins(GranularDefaultPlugins);
/// app.start();
/// ```
pub struct GranularDefaultPlugins;

impl PluginGroup for GranularDefaultPlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(DemAtomPlugin)
            .add(DemAtomInsertPlugin)
            .add(VelocityVerletPlugin::new())
            .add(HertzMindlinContactPlugin)
            .add(RotationalDynamicsPlugin)
    }
}
