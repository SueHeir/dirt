//! # dirt_atom — Per-atom DEM data and material property tables
//!
//! This crate provides the core per-atom data structures and material property
//! system for Discrete Element Method (DEM) simulations in DIRT:
//!
//! - **[`DemAtom`]** — Per-atom extension data (radius, density, angular velocity, torque, etc.)
//!   registered via the `AtomData` derive macro with `#[forward]`, `#[reverse]`, and `#[zero]`
//!   attributes that control MPI pack/unpack and per-timestep zeroing behavior.
//! - **[`MaterialTable`]** — Named materials with mechanical properties (Young's modulus,
//!   Poisson's ratio, restitution, friction) and precomputed per-pair mixing tables
//!   (geometric-mean or harmonic-mean) for efficient contact force evaluation.
//! - **[`DemAtomInsertPlugin`]** — Particle insertion from `[[particles.insert]]` TOML config:
//!   random placement, rate-based trickle insertion, or file-based loading (CSV, LAMMPS dump/data).
//! - **[`RadiusSpec`]** — Particle radius specification: fixed value or statistical distribution
//!   (uniform, gaussian, lognormal, discrete).
//!
//! ## TOML Configuration
//!
//! Materials are defined under `[[dem.materials]]`:
//!
//! ```toml
//! [dem]
//! contact_model = "hertz"  # "hertz" (default), "hooke", or "mdr"
//!
//! [[dem.materials]]
//! name = "glass"
//! youngs_mod = 8.7e9       # Young's modulus (Pa)
//! poisson_ratio = 0.3      # Poisson's ratio (dimensionless, 0–0.5)
//! restitution = 0.95       # Coefficient of restitution (0–1)
//! friction = 0.4           # Coulomb friction coefficient (default 0.4)
//! ```
//!
//! ## AtomData Derive Attributes
//!
//! The `#[forward]`, `#[reverse]`, and `#[zero]` attributes on [`DemAtom`] fields
//! control how data is handled during MPI communication and timestep integration:
//!
//! - **`#[forward]`** — Field is packed and sent forward from owning processor to ghost
//!   atoms on neighboring processors (e.g., angular velocity `omega`).
//! - **`#[reverse]`** — Field is accumulated in reverse: ghost contributions are sent back
//!   and summed into the owning atom's value (e.g., `torque`).
//! - **`#[zero]`** — Field is zeroed at the start of each force computation step (e.g., `torque`),
//!   before new contact forces accumulate.
//!
//! ## Where this fits in an app
//!
//! `dirt_atom` is registered for you by [`DemAtomPlugin`] / [`DemAtomInsertPlugin`],
//! which are bundled in `dirt_granular::GranularDefaultPlugins`. You rarely
//! construct a `MaterialTable` directly in application code — it is built from
//! the `[[dem.materials]]` TOML section — but the
//! [`MaterialTable::new`] docs show the by-hand two-phase build for tests and
//! tools. For the smallest complete runnable application (drop particles in a
//! box under gravity), see the **`hello_bed`** example in the DIRT repo's
//! `examples/` directory, which assembles `CorePlugins + GranularDefaultPlugins`
//! and lets this crate's insertion + material code run from a config file.

// Public API documentation-completeness gate: every public item in this crate
// must carry a doc comment. Enforced on both `cargo build` (rustc) and
// `cargo doc` (rustdoc; e.g. `RUSTDOCFLAGS="-D missing_docs"`). Document real
// API intent here — do not add empty doc comments just to satisfy the lint.
#![deny(missing_docs)]

pub mod insert;
pub mod radius;

pub use insert::*;
pub use radius::*;

use std::f64::consts::PI;

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use serde::Deserialize;
use soil_derive::AtomData;

use soil_core::{register_atom_data, Atom, AtomData, AtomPlugin, Config, ScheduleSetupSet};

// ── Shared physics constants ────────────────────────────────────────────────

/// Precomputed `sqrt(5/6)`, the coefficient relating tangential to normal
/// damping in the Hertz–Mindlin contact model.
pub const SQRT_5_6: f64 = 0.9128709291752768;

// ── Tsuji (1992) restitution → damping mapping (Hertz) ───────────────────────

/// Tsuji–Tanaka–Ishida (1992) damping coefficient α(e) for a nonlinear Hertz
/// contact. This is the **same polynomial LAMMPS uses** for `damping tsuji`
/// (`src/GRANULAR/gran_sub_mod_damping.cpp`, `GranSubModDampingTsuji::init`):
///
/// `α(e) = 1.2728 − 4.2783 e + 11.087 e² − 22.348 e³ + 27.467 e⁴ − 18.022 e⁵ + 4.8218 e⁶`
///
/// where `e` is the restitution input (the third `hertz/material` argument in
/// LAMMPS). α(e) is the damping prefactor in `F_damp = α(e)·√(mₑ Fₙ/δ)·vₙ`.
fn tsuji_alpha(e: f64) -> f64 {
    1.2728 - 4.2783 * e + 11.087 * e.powi(2) - 22.348 * e.powi(3) + 27.467 * e.powi(4)
        - 18.022 * e.powi(5)
        + 4.8218 * e.powi(6)
}

/// COR of a head-on Hertz collision with the damping DIRT applies
/// (`f_diss = 2β√(5/6)√(Sₙ mᵣ) vₙ`, `Sₙ = 2E*√(R*δ)`), as a function of β.
///
/// Computed by integrating the dimensionless 1D collision (E*=R*=mᵣ=1, v₀=1):
/// `δ̈ = −(4/3)δ^{3/2} − 2β√(5/6)√2 · δ^{1/4} · δ̇`. The Tsuji scaling makes this
/// velocity-independent, so one integration fixes the β↔COR map. RK4, fixed dt.
/// Used by the tests to verify DIRT's realized COR reproduces LAMMPS `damping tsuji`.
#[cfg(test)]
fn hertz_cor_of_beta(beta: f64) -> f64 {
    if beta <= 0.0 {
        return 1.0;
    }
    let c = 2.0 * beta * SQRT_5_6 * std::f64::consts::SQRT_2; // damping prefactor
                                                              // acceleration of the overlap coordinate (only while in contact, δ>0)
    let acc = |d: f64, v: f64| -> f64 {
        if d <= 0.0 {
            0.0
        } else {
            -(4.0 / 3.0) * d.powf(1.5) - c * d.powf(0.25) * v
        }
    };
    let dt = 1.0e-4;
    // δ = overlap (grows on approach): start at contact with δ̇ = +1 (approaching).
    let (mut d, mut v) = (0.0_f64, 1.0_f64);
    for _ in 0..2_000_000 {
        // RK4 on (d, v) with d' = v, v' = acc(d, v)
        let (k1d, k1v) = (v, acc(d, v));
        let (k2d, k2v) = (
            v + 0.5 * dt * k1v,
            acc(d + 0.5 * dt * k1d, v + 0.5 * dt * k1v),
        );
        let (k3d, k3v) = (
            v + 0.5 * dt * k2v,
            acc(d + 0.5 * dt * k2d, v + 0.5 * dt * k2v),
        );
        let (k4d, k4v) = (v + dt * k3v, acc(d + dt * k3d, v + dt * k3v));
        d += dt / 6.0 * (k1d + 2.0 * k2d + 2.0 * k3d + k4d);
        v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        if d <= 0.0 && v < 0.0 {
            return v.abs(); // separated (overlap back to 0, receding); COR = |v_out| (v_in = 1)
        }
    }
    v.abs()
}

/// Damping ratio β for DIRT's Hertz viscoelastic damping that **reproduces
/// LAMMPS `damping tsuji`** for a restitution input `e_target`.
///
/// DIRT applies `F_damp = 2β√(5/6)·√(Sₙ mᵣ)·vₙ` with `Sₙ = 2E*√(R*δ)`, while
/// LAMMPS `damping tsuji` applies `F_damp = α(e)·√(mₑ Fₙ/δ)·vₙ` with the same
/// Hertz stiffness (`Fₙ/δ = (4/3)E*√R*·δ^{1/2}`). Equating the two forces
/// (`mₑ = mᵣ` for a binary collision) gives, after the δ^{1/4} terms cancel,
///
/// `2β√(5/6)·√2 = α(e)·√(4/3)`  ⇒  `β = α(e)/√5`.
///
/// So feeding both codes the **same nominal restitution** makes them integrate
/// an *identical* damping force and realize the same (velocity-independent) COR
/// to integrator precision — the cross-code check in `bench_hertz_rebound`.
/// (`e` is the Tsuji restitution parameter, LAMMPS's convention; because the
/// Tsuji polynomial is a fit, the realized COR sits slightly above `e` below
/// e≈0.9 in *both* codes — a shared property, not a code error.)
pub fn hertz_beta_for_cor(e_target: f64) -> f64 {
    if e_target >= 0.9999 {
        return 0.0;
    }
    let e = e_target.clamp(1.0e-3, 0.9999);
    tsuji_alpha(e) / 5.0_f64.sqrt()
}

// ── Config structs ──────────────────────────────────────────────────────────

fn default_friction() -> f64 {
    0.4
}

fn default_contact_model() -> String {
    "hertz".to_string()
}

#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
/// A single material definition from `[[dem.materials]]`.
pub struct MaterialConfig {
    /// Material name, referenced by particle insert blocks.
    pub name: String,
    /// Young's modulus (Pa).
    pub youngs_mod: f64,
    /// Poisson's ratio (dimensionless, 0–0.5).
    pub poisson_ratio: f64,
    /// Coefficient of restitution (0–1).
    pub restitution: f64,
    /// Coulomb friction coefficient.
    #[serde(default = "default_friction")]
    pub friction: f64,
    /// Rolling friction coefficient (0 = disabled).
    #[serde(default)]
    pub rolling_friction: f64,
    /// Cohesion energy density (J/m², 0 = disabled). SJKR model.
    #[serde(default)]
    pub cohesion_energy: f64,
    /// Surface energy (J/m², 0 = disabled). Activates JKR adhesion model.
    /// Cannot be used together with cohesion_energy on the same material.
    #[serde(default)]
    pub surface_energy: f64,
    /// Twisting friction coefficient (0 = disabled).
    #[serde(default)]
    pub twisting_friction: f64,
    /// Linear normal stiffness for Hooke model (N/m, 0 = use Hertz).
    #[serde(default)]
    pub kn: f64,
    /// Linear tangential stiffness for Hooke model (N/m, 0 = use Mindlin).
    #[serde(default)]
    pub kt: f64,
    /// Rolling spring stiffness for SDS rolling model (N·m/rad, 0 = use constant model).
    #[serde(default)]
    pub rolling_stiffness: f64,
    /// Rolling viscous damping coefficient for SDS rolling model.
    #[serde(default)]
    pub rolling_damping: f64,
    /// Twisting spring stiffness for SDS twisting model (N·m/rad, 0 = use constant model).
    #[serde(default)]
    pub twisting_stiffness: f64,
    /// Twisting viscous damping coefficient for SDS twisting model.
    #[serde(default)]
    pub twisting_damping: f64,
    /// Yield stress for the MDR elastic-plastic normal model (Pa).
    #[serde(default)]
    pub mdr_yield_stress: f64,
    /// Critical free-surface ratio for the full LAMMPS MDR bulk branch.
    ///
    /// DIRT stores and documents this value for LAMMPS input parity. The current
    /// pairwise MDR implementation does not update apparent particle radii or
    /// activate the multi-contact bulk response, so this parameter is not used
    /// in force evaluation yet.
    #[serde(default)]
    pub mdr_psi_b: f64,
    /// Dimensionless MDR normal damping prefactor.
    #[serde(default)]
    pub mdr_damping: f64,
    /// Liquid bridge volume per contact (m^3, 0 = disabled).
    #[serde(default)]
    pub liquid_bridge_volume: f64,
    /// Liquid-vapor surface tension for pendular bridges (N/m).
    #[serde(default)]
    pub liquid_surface_tension: f64,
    /// Solid-liquid contact angle for pendular bridges (radians).
    #[serde(default)]
    pub liquid_contact_angle: f64,
    /// Optional rupture distance for liquid bridges (m). If 0, the Willett/Lian
    /// volume-scaled estimate is used.
    #[serde(default)]
    pub liquid_rupture_distance: f64,
}

fn default_adhesion_model() -> String {
    "jkr".to_string()
}

fn default_rolling_model() -> String {
    "constant".to_string()
}

fn default_tangential_model() -> String {
    "history".to_string()
}

fn default_twisting_model() -> String {
    "constant".to_string()
}

fn default_limit_damping() -> bool {
    true
}

fn default_liquid_bridge_model() -> String {
    "off".to_string()
}

/// TOML `[dem]` — top-level DEM configuration containing material definitions.
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct DemConfig {
    /// List of material definitions, each corresponding to a `[[dem.materials]]` block.
    pub materials: Option<Vec<MaterialConfig>>,
    /// Contact model: "hertz" (default), "hooke", or "mdr".
    #[serde(default = "default_contact_model")]
    pub contact_model: String,
    /// Adhesion model when `surface_energy > 0`: "jkr" (default) or "dmt".
    ///
    /// - **JKR** (Johnson-Kendall-Roberts): Modified contact area with pull-off force = 1.5*pi*gamma*R*.
    ///   Suitable for soft materials with high surface energy (Tabor parameter > 5).
    /// - **DMT** (Derjaguin-Muller-Toporov): Pure Hertz contact area with constant attractive
    ///   force = 2*pi*gamma*R*. Suitable for stiff materials with low surface energy (Tabor parameter < 0.1).
    #[serde(default = "default_adhesion_model")]
    pub adhesion_model: String,
    /// Rolling resistance model: "constant" (default) or "sds".
    #[serde(default = "default_rolling_model")]
    pub rolling_model: String,
    /// Twisting friction model: "constant" (default) or "sds".
    #[serde(default = "default_twisting_model")]
    pub twisting_model: String,
    /// Tangential friction model: "history" (default), "linear_nohistory",
    /// "mindlin_rescale", or "mindlin_rescale/force".
    ///
    /// - **history** (default): the incremental Mindlin tangential spring — the
    ///   force accumulates the sliding displacement `ξ` over the contact lifetime
    ///   (`F_t = k_t ξ + γ_t v_t`, Coulomb-capped). This is DIRT's historical path
    ///   and mirrors LAMMPS `pair_granular … tangential linear_history/mindlin`.
    /// - **linear_nohistory**: a history-free velocity-Coulomb law mirroring LAMMPS
    ///   `pair_granular … tangential linear_nohistory` (and the classic
    ///   `pair gran/hooke`). The tangential force is `F_t = -min(μ |F_n|, γ_t |v_t|) t̂`
    ///   with **no** accumulated spring displacement — the force depends only on the
    ///   instantaneous relative tangential velocity. Rolling/twisting are unaffected.
    /// - **mindlin_rescale**: the Mindlin displacement-history spring, with the
    ///   LAMMPS unloading rule `ξ <- ξ * a/a_prev` when contact radius decreases.
    /// - **mindlin_rescale/force** (alias: `mindlin_rescale_force`): stores the
    ///   elastic tangential force as history and applies the same unloading gate to
    ///   that force, mirroring LAMMPS `mindlin_rescale/force`.
    #[serde(default = "default_tangential_model")]
    pub tangential_model: String,
    /// Track per-sphere orientation (quaternion). Default `false`.
    ///
    /// A sphere is rotationally symmetric, so its orientation never enters any
    /// contact force law (those depend on angular velocity `ω`, not absolute
    /// orientation). The per-sphere quaternion is therefore causally inert for
    /// pure-sphere runs, and integrating it each step (sqrt + division + sin +
    /// cos + Hamilton product per atom) is pure overhead. Leave this `false`
    /// unless something downstream actually reads the orientation (e.g. a future
    /// surface-marker visualization). Non-spherical bodies track orientation in
    /// their own `BodyData`, not here, so this flag does not affect them.
    #[serde(default)]
    pub track_orientation: bool,
    /// Clamp the total normal force to be repulsive-only (≥ 0), i.e. forbid the
    /// viscoelastic damping term from producing a net *attractive* force near the
    /// end of a collision. Mirrors LAMMPS `pair granular … limit_damping`.
    ///
    /// Default `true` (DIRT's historical behavior). Set `false` to reproduce
    /// LAMMPS's **default** viscoelastic contact (no tensile cutoff), which is
    /// required for exact cross-code COR agreement at low restitution — below
    /// e≈0.9 the damping impulse near separation is large enough that clamping it
    /// away raises the realized COR by up to ~0.03. See `bench_hertz_rebound`.
    #[serde(default = "default_limit_damping")]
    pub limit_damping: bool,
    /// Pendular capillary bridge model: `"off"` (default) or `"willett2000"`.
    #[serde(default = "default_liquid_bridge_model")]
    pub liquid_bridge_model: String,
}

impl Default for DemConfig {
    fn default() -> Self {
        DemConfig {
            materials: None,
            contact_model: default_contact_model(),
            adhesion_model: default_adhesion_model(),
            rolling_model: default_rolling_model(),
            twisting_model: default_twisting_model(),
            tangential_model: default_tangential_model(),
            track_orientation: false,
            limit_damping: default_limit_damping(),
            liquid_bridge_model: default_liquid_bridge_model(),
        }
    }
}

/// Diagnostic for the Hooke/Hertz adhesion asymmetry footgun.
///
/// JKR/DMT adhesion (driven by per-material `surface_energy`) is implemented
/// **only** on the Hertz contact path; the `contact_model = "hooke"`
/// linear-spring path ignores `surface_energy` entirely (see the module docs of
/// `dirt_granular` and `dirt_granular::contact`). Under `contact_model = "hooke"`
/// a nonzero `surface_energy` is therefore *silently dropped*, which reads as a
/// physics change that never happens.
///
/// Returns `Some(message)` when the loaded config is in that silent-drop state —
/// `contact_model == "hooke"` and at least one material sets `surface_energy > 0`
/// — naming every offending material and pointing at the Hertz path. Returns
/// `None` otherwise (Hertz with any `surface_energy`, or Hooke with all
/// `surface_energy == 0`), so a clean config produces no diagnostic. This
/// function is pure (no I/O); the plugin build wires it to `eprintln!`.
///
/// The `"hooke"` comparison is exact and mirrors the runtime dispatch in
/// `HertzMindlinContactPlugin` (any non-`"hooke"` string selects the Hertz path).
pub fn hooke_surface_energy_warning(config: &DemConfig) -> Option<String> {
    if config.contact_model != "hooke" {
        return None;
    }
    let offenders: Vec<&str> = match config.materials {
        Some(ref mats) => mats
            .iter()
            .filter(|m| m.surface_energy > 0.0)
            .map(|m| m.name.as_str())
            .collect(),
        None => Vec::new(),
    };
    if offenders.is_empty() {
        return None;
    }
    Some(format!(
        "WARNING: contact_model = \"hooke\" ignores `surface_energy` \
         (JKR/DMT adhesion is only implemented on the Hertz contact path). \
         surface_energy > 0 on material(s) [{}] will be SILENTLY DROPPED — no \
         adhesion/pull-off force will be applied. To get JKR/DMT adhesion set \
         contact_model = \"hertz\" (the default); otherwise set surface_energy = 0 \
         on these material(s) to silence this warning. For linear-spring cohesion \
         under Hooke, use `cohesion_energy` (SJKR) instead.",
        offenders.join(", ")
    ))
}

fn liquid_bridge_model_error(config: &DemConfig) -> Option<String> {
    match config.liquid_bridge_model.as_str() {
        "off" | "willett2000" => None,
        other => Some(format!(
            "ERROR: invalid [dem].liquid_bridge_model = {:?}. Supported values are \
             \"off\" and \"willett2000\".",
            other
        )),
    }
}

// ── MaterialTable — per-material and per-pair precomputed properties ────────

/// Per-material properties and precomputed per-pair mixing tables for contact force evaluation.
///
/// # Two-phase build contract
///
/// A `MaterialTable` is filled in **two phases**, and the contact force code
/// reads only the second-phase output:
///
/// 1. **Register materials.** Each `add_material*` call appends one row to the
///    per-material vectors (`youngs_mod`, `restitution`, `friction`, …) and
///    returns its integer index. During this phase **every `*_ij` pair table is
///    empty** (`Vec::new()`).
/// 2. **Build pair tables.** [`build_pair_tables`](Self::build_pair_tables)
///    allocates the `N×N` `*_ij` tables and fills them from the registered
///    per-material values using the mixing rules below. It **must** be called
///    once, after the last material is added and before any contact force is
///    evaluated. Indexing a `*_ij` table before this is an out-of-bounds panic.
///
/// Mixing rules used in phase 2:
/// - **Geometric mean** for friction, restitution, cohesion/surface energy, rolling/twisting friction
/// - **Harmonic mean** (2·ki·kj/(ki+kj)) for Hooke stiffnesses and SDS spring stiffnesses
/// - **Effective modulus** formulas for Hertz (`E*`) and Mindlin (`G*`) contact models
///
/// Per-material rows are indexed by the index returned from `add_material*` (or
/// looked up via [`find_material`](Self::find_material)); pair tables are indexed
/// `table_ij[i][j]`.
///
/// # The `add_material*` ladder
///
/// Four constructors form a wrapping ladder from fewest to most arguments; each
/// delegates to the next with sensible zero defaults:
///
/// - [`add_material`](Self::add_material) — basics (E, ν, restitution, friction,
///   rolling friction, cohesion); sets `surface_energy = 0` (no adhesion).
/// - [`add_material_full`](Self::add_material_full) — adds `surface_energy`.
/// - [`add_material_extended`](Self::add_material_extended) — adds twisting
///   friction and Hooke linear stiffnesses `kn`/`kt`.
/// - [`add_material_with_sds`](Self::add_material_with_sds) — adds the SDS
///   rolling/twisting spring–dashpot parameters; the full constructor.
///
/// Use the shortest one that covers the parameters you need.
///
/// # Restitution → damping
///
/// `restitution` is stored as the **restitution input `e`** (the same parameter
/// LAMMPS's `hertz/material` takes), not a damping ratio. In phase 2,
/// `beta_ij[i][j]` is derived via [`hertz_beta_for_cor`] as `β = α(e)/√5` with the
/// Tsuji (1992) polynomial α(e), so DIRT's Hertz normal-damping force is IDENTICAL
/// to LAMMPS `damping tsuji`. Both codes then realize the same velocity-independent
/// COR for the same nominal `e` (validated by `bench_hertz_rebound`), keeping
/// DIRT's shear/cooling results on the same cross-code line as LAMMPS/LIGGGHTS.
/// Because the Tsuji polynomial is a fit, the realized COR sits slightly above the
/// nominal `e` below e≈0.9 — a property shared with LAMMPS, not a code error.
///
/// # Config-error convention
///
/// `add_material*` validates physically inconsistent input (e.g. both
/// `cohesion_energy` and `surface_energy` set) with [`MaterialError`].
/// [`DemAtomPlugin`] maps that error to the fallible app boundary, allowing a
/// runner to report the malformed table without terminating the process.
///
/// # Hooke vs. Hertz adhesion asymmetry
///
/// JKR/DMT adhesion (`surface_energy`) is only honored under the Hertz contact
/// model. Under `contact_model = "hooke"` the surface-energy term is silently
/// ignored — only SJKR-style cohesion (`cohesion_energy`) is applied. See
/// `dirt_granular` for the per-branch parameter reference.
pub struct MaterialTable {
    /// Material names, indexed by material ID (the index into every per-material vector).
    pub names: Vec<String>,
    /// Per-material Young's modulus E (Pa).
    pub youngs_mod: Vec<f64>,
    /// Per-material Poisson ratio ν (dimensionless).
    pub poisson_ratio: Vec<f64>,
    /// Per-material sliding friction coefficient μ (dimensionless).
    pub friction: Vec<f64>,
    /// Per-material coefficient of restitution e (dimensionless).
    pub restitution: Vec<f64>,
    /// Per-material rolling friction coefficient (dimensionless).
    pub rolling_friction: Vec<f64>,
    /// Per-material twisting friction coefficient (dimensionless).
    pub twisting_friction: Vec<f64>,
    /// Per-material SJKR cohesion energy density (J/m³).
    pub cohesion_energy: Vec<f64>,
    /// Per-material JKR/DMT surface energy γ (J/m²).
    pub surface_energy: Vec<f64>,
    /// Per-pair Tsuji damping coefficient β, derived from the pair restitution.
    pub beta_ij: Vec<Vec<f64>>,
    /// Per-pair sliding friction coefficient (overrides the geometric mean when set).
    pub friction_ij: Vec<Vec<f64>>,
    /// Per-pair rolling friction coefficient.
    pub rolling_friction_ij: Vec<Vec<f64>>,
    /// Per-pair SJKR cohesion energy density (J/m³).
    pub cohesion_energy_ij: Vec<Vec<f64>>,
    /// Per-pair surface energy for JKR adhesion (geometric mean mixing).
    pub surface_energy_ij: Vec<Vec<f64>>,
    /// Precomputed effective Young's modulus for each material pair (Hertz contact).
    pub e_eff_ij: Vec<Vec<f64>>,
    /// Precomputed effective shear modulus for each material pair (Mindlin contact).
    pub g_eff_ij: Vec<Vec<f64>>,
    /// Per-pair twisting friction (geometric mean mixing).
    pub twisting_friction_ij: Vec<Vec<f64>>,
    /// Per-material linear normal stiffness for Hooke model.
    pub kn: Vec<f64>,
    /// Per-material linear tangential stiffness for Hooke model.
    pub kt: Vec<f64>,
    /// Per-pair Hooke normal stiffness (harmonic mean: 2*ki*kj/(ki+kj)).
    pub kn_ij: Vec<Vec<f64>>,
    /// Per-pair Hooke tangential stiffness (harmonic mean).
    pub kt_ij: Vec<Vec<f64>>,
    /// Contact model: "hertz" or "hooke".
    pub contact_model: String,
    /// Adhesion model: "jkr" (default) or "dmt".
    pub adhesion_model: String,
    /// Rolling resistance model: "constant" or "sds".
    pub rolling_model: String,
    /// Twisting friction model: "constant" or "sds".
    pub twisting_model: String,
    /// Tangential friction model. See [`DemConfig::tangential_model`].
    pub tangential_model: String,
    /// Track per-sphere orientation (quaternion). Default `false`; see [`DemConfig::track_orientation`].
    pub track_orientation: bool,
    /// Clamp total normal force to repulsive-only (≥ 0). Default `true`; see
    /// [`DemConfig::limit_damping`]. `false` reproduces LAMMPS's default
    /// (no tensile cutoff) for exact cross-code COR at low restitution.
    pub limit_damping: bool,
    /// Per-material rolling spring stiffness (SDS model).
    pub rolling_stiffness: Vec<f64>,
    /// Per-material rolling damping coefficient (SDS model).
    pub rolling_damping: Vec<f64>,
    /// Per-material twisting spring stiffness (SDS model).
    pub twisting_stiffness: Vec<f64>,
    /// Per-material twisting damping coefficient (SDS model).
    pub twisting_damping: Vec<f64>,
    /// Per-material MDR yield stress (Pa).
    pub mdr_yield_stress: Vec<f64>,
    /// Per-material MDR critical free-surface ratio.
    pub mdr_psi_b: Vec<f64>,
    /// Per-material MDR damping prefactor.
    pub mdr_damping: Vec<f64>,
    /// Per-material liquid bridge volume (m^3).
    pub liquid_bridge_volume: Vec<f64>,
    /// Per-material liquid-vapor surface tension (N/m).
    pub liquid_surface_tension: Vec<f64>,
    /// Per-material liquid bridge contact angle (radians).
    pub liquid_contact_angle: Vec<f64>,
    /// Per-material liquid bridge rupture distance override (m).
    pub liquid_rupture_distance: Vec<f64>,
    /// Per-pair rolling stiffness (harmonic mean).
    pub rolling_stiffness_ij: Vec<Vec<f64>>,
    /// Per-pair rolling damping (geometric mean).
    pub rolling_damping_ij: Vec<Vec<f64>>,
    /// Per-pair twisting stiffness (harmonic mean).
    pub twisting_stiffness_ij: Vec<Vec<f64>>,
    /// Per-pair twisting damping (geometric mean).
    pub twisting_damping_ij: Vec<Vec<f64>>,
    /// Per-pair MDR yield stress (geometric mean).
    pub mdr_yield_stress_ij: Vec<Vec<f64>>,
    /// Per-pair MDR critical free-surface ratio (arithmetic mean).
    pub mdr_psi_b_ij: Vec<Vec<f64>>,
    /// Per-pair MDR damping prefactor (geometric mean).
    pub mdr_damping_ij: Vec<Vec<f64>>,
    /// Per-pair liquid bridge volume (geometric mean, m^3).
    pub liquid_bridge_volume_ij: Vec<Vec<f64>>,
    /// Per-pair liquid-vapor surface tension (geometric mean, N/m).
    pub liquid_surface_tension_ij: Vec<Vec<f64>>,
    /// Per-pair liquid bridge contact angle (arithmetic mean, radians).
    pub liquid_contact_angle_ij: Vec<Vec<f64>>,
    /// Per-pair liquid bridge rupture distance override (geometric mean, m).
    pub liquid_rupture_distance_ij: Vec<Vec<f64>>,
    /// Pendular capillary bridge model: "off" or "willett2000".
    pub liquid_bridge_model: String,
}

/// Error returned when a material definition cannot be represented by a
/// [`MaterialTable`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MaterialError {
    /// A material enabled mutually exclusive cohesion models.
    ConflictingCohesion {
        /// Name of the rejected material.
        name: String,
    },
    /// A named material property was non-finite or outside the model's domain.
    InvalidProperty {
        /// Name of the rejected material.
        name: String,
        /// Fully qualified input property name.
        property: &'static str,
        /// Physical domain required by the contact models.
        requirement: &'static str,
    },
}

/// Stable index of a material registered in a [`MaterialTable`].
///
/// The raw integer remains available through [`MaterialId::index`] for atom
/// type fields, but named APIs use this type to avoid mixing material IDs with
/// unrelated counts or indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MaterialId(u32);

impl MaterialId {
    /// Returns the zero-based table index used by hot-path pair tables.
    pub const fn index(self) -> usize {
        self.0 as usize
    }

    /// Returns the raw atom material type value.
    pub const fn raw(self) -> u32 {
        self.0
    }
}

/// Elastic material constants shared by Hertz, Mindlin, and Hooke contacts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elastic {
    /// Young's modulus (Pa).
    pub youngs_mod: f64,
    /// Poisson ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Coefficient of restitution (dimensionless).
    pub restitution: f64,
    /// Linear Hooke normal stiffness (N/m); zero selects Hertz.
    pub normal_stiffness: f64,
    /// Linear Hooke tangential stiffness (N/m); zero selects Mindlin.
    pub tangential_stiffness: f64,
}

impl Elastic {
    /// Creates Hertz/Mindlin elastic constants with no Hooke stiffness override.
    pub const fn new(youngs_mod: f64, poisson_ratio: f64, restitution: f64) -> Self {
        Self {
            youngs_mod,
            poisson_ratio,
            restitution,
            normal_stiffness: 0.0,
            tangential_stiffness: 0.0,
        }
    }

    /// Adds linear Hooke normal and tangential stiffnesses.
    pub const fn with_hooke_stiffness(mut self, normal: f64, tangential: f64) -> Self {
        self.normal_stiffness = normal;
        self.tangential_stiffness = tangential;
        self
    }
}

/// Sliding and torsional Coulomb friction coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Friction {
    /// Sliding Coulomb friction coefficient.
    pub sliding: f64,
    /// Rolling Coulomb friction coefficient.
    pub rolling: f64,
    /// Twisting Coulomb friction coefficient.
    pub twisting: f64,
}

impl Default for Friction {
    fn default() -> Self {
        Self {
            sliding: default_friction(),
            rolling: 0.0,
            twisting: 0.0,
        }
    }
}

/// Adhesion choice for one material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Adhesion {
    /// No adhesive contribution.
    None,
    /// SJKR cohesion energy density (J/m³).
    Sjkr {
        /// Cohesion energy density.
        energy: f64,
    },
    /// JKR/DMT surface energy (J/m²), selected by the table's adhesion model.
    SurfaceEnergy {
        /// Surface energy.
        energy: f64,
    },
}

impl Default for Adhesion {
    fn default() -> Self {
        Self::None
    }
}

/// Optional spring-dashpot rolling model constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Rolling {
    /// The constant rolling-friction model.
    Constant,
    /// SDS rolling spring stiffness (N m/rad) and damping coefficient.
    Sds {
        /// Rolling spring stiffness.
        stiffness: f64,
        /// Rolling damping coefficient.
        damping: f64,
    },
}

impl Default for Rolling {
    fn default() -> Self {
        Self::Constant
    }
}

/// Optional spring-dashpot twisting model constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Twisting {
    /// The constant twisting-friction model.
    Constant,
    /// SDS twisting spring stiffness (N m/rad) and damping coefficient.
    Sds {
        /// Twisting spring stiffness.
        stiffness: f64,
        /// Twisting damping coefficient.
        damping: f64,
    },
}

impl Default for Twisting {
    fn default() -> Self {
        Self::Constant
    }
}

/// MDR elastic-plastic parameters.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Mdr {
    /// Yield stress (Pa).
    pub yield_stress: f64,
    /// Critical free-surface ratio.
    pub psi_b: f64,
    /// Dimensionless normal damping prefactor.
    pub damping: f64,
}

/// Pendular liquid-bridge parameters.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LiquidBridge {
    /// Bridge volume per contact (m³).
    pub volume: f64,
    /// Liquid-vapor surface tension (N/m).
    pub surface_tension: f64,
    /// Solid-liquid contact angle (radians).
    pub contact_angle: f64,
    /// Optional rupture-distance override (m); zero uses the model estimate.
    pub rupture_distance: f64,
}

/// Complete, named input for one material row.
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    /// Human-readable material name.
    pub name: String,
    /// Elastic constants.
    pub elastic: Elastic,
    /// Friction constants.
    pub friction: Friction,
    /// Adhesion selection.
    pub adhesion: Adhesion,
    /// Rolling model constants.
    pub rolling: Rolling,
    /// Twisting model constants.
    pub twisting: Twisting,
    /// MDR constants.
    pub mdr: Mdr,
    /// Liquid-bridge constants.
    pub liquid_bridge: LiquidBridge,
}

impl Material {
    /// Starts a material definition with required elastic constants and named defaults.
    pub fn new(name: impl Into<String>, elastic: Elastic) -> Self {
        Self {
            name: name.into(),
            elastic,
            friction: Friction::default(),
            adhesion: Adhesion::None,
            rolling: Rolling::Constant,
            twisting: Twisting::Constant,
            mdr: Mdr::default(),
            liquid_bridge: LiquidBridge::default(),
        }
    }

    /// Sets sliding, rolling, and twisting friction coefficients.
    pub const fn with_friction(mut self, friction: Friction) -> Self {
        self.friction = friction;
        self
    }
    /// Sets the material adhesion model and its parameter.
    pub const fn with_adhesion(mut self, adhesion: Adhesion) -> Self {
        self.adhesion = adhesion;
        self
    }
    /// Sets rolling model constants.
    pub const fn with_rolling(mut self, rolling: Rolling) -> Self {
        self.rolling = rolling;
        self
    }
    /// Sets twisting model constants.
    pub const fn with_twisting(mut self, twisting: Twisting) -> Self {
        self.twisting = twisting;
        self
    }
    /// Sets MDR elastic-plastic constants.
    pub const fn with_mdr(mut self, mdr: Mdr) -> Self {
        self.mdr = mdr;
        self
    }
    /// Sets liquid-bridge constants.
    pub const fn with_liquid_bridge(mut self, liquid_bridge: LiquidBridge) -> Self {
        self.liquid_bridge = liquid_bridge;
        self
    }
}

impl Material {
    fn from_config(config: &MaterialConfig) -> Result<Self, MaterialError> {
        if config.cohesion_energy > 0.0 && config.surface_energy > 0.0 {
            return Err(MaterialError::ConflictingCohesion {
                name: config.name.clone(),
            });
        }
        Ok(Material::new(
            &config.name,
            Elastic::new(config.youngs_mod, config.poisson_ratio, config.restitution)
                .with_hooke_stiffness(config.kn, config.kt),
        )
        .with_friction(Friction {
            sliding: config.friction,
            rolling: config.rolling_friction,
            twisting: config.twisting_friction,
        })
        .with_adhesion(if config.cohesion_energy > 0.0 {
            Adhesion::Sjkr {
                energy: config.cohesion_energy,
            }
        } else if config.surface_energy > 0.0 {
            Adhesion::SurfaceEnergy {
                energy: config.surface_energy,
            }
        } else {
            Adhesion::None
        })
        .with_rolling(Rolling::Sds {
            stiffness: config.rolling_stiffness,
            damping: config.rolling_damping,
        })
        .with_twisting(Twisting::Sds {
            stiffness: config.twisting_stiffness,
            damping: config.twisting_damping,
        })
        .with_mdr(Mdr {
            yield_stress: config.mdr_yield_stress,
            psi_b: config.mdr_psi_b,
            damping: config.mdr_damping,
        })
        .with_liquid_bridge(LiquidBridge {
            volume: config.liquid_bridge_volume,
            surface_tension: config.liquid_surface_tension,
            contact_angle: config.liquid_contact_angle,
            rupture_distance: config.liquid_rupture_distance,
        }))
    }
}

/// Reject input that would make a constitutive or mixing calculation undefined.
///
/// This is deliberately performed before `MaterialTable::add` touches any
/// column, so a failed registration cannot leave the parallel columns out of
/// sync.  The bounds are the isotropic, non-negative parameter domains used by
/// DIRT's DEM models (rather than an attempt to represent every exotic material).
fn validate_material(material: &Material) -> Result<(), MaterialError> {
    let invalid = |property, requirement| MaterialError::InvalidProperty {
        name: material.name.clone(),
        property,
        requirement,
    };
    let finite = |value: f64| value.is_finite();
    let nonnegative = |value: f64| finite(value) && value >= 0.0;

    if !finite(material.elastic.youngs_mod) || material.elastic.youngs_mod <= 0.0 {
        return Err(invalid("elastic.youngs_mod", "finite and > 0"));
    }
    if !finite(material.elastic.poisson_ratio)
        || !(0.0..0.5).contains(&material.elastic.poisson_ratio)
    {
        return Err(invalid("elastic.poisson_ratio", "finite and in [0, 0.5)"));
    }
    if !finite(material.elastic.restitution) || !(0.0..=1.0).contains(&material.elastic.restitution)
    {
        return Err(invalid("elastic.restitution", "finite and in [0, 1]"));
    }
    for (property, value) in [
        (
            "elastic.normal_stiffness",
            material.elastic.normal_stiffness,
        ),
        (
            "elastic.tangential_stiffness",
            material.elastic.tangential_stiffness,
        ),
        ("friction.sliding", material.friction.sliding),
        ("friction.rolling", material.friction.rolling),
        ("friction.twisting", material.friction.twisting),
        ("mdr.yield_stress", material.mdr.yield_stress),
        ("mdr.damping", material.mdr.damping),
        ("liquid_bridge.volume", material.liquid_bridge.volume),
        (
            "liquid_bridge.surface_tension",
            material.liquid_bridge.surface_tension,
        ),
        (
            "liquid_bridge.rupture_distance",
            material.liquid_bridge.rupture_distance,
        ),
    ] {
        if !nonnegative(value) {
            return Err(invalid(property, "finite and >= 0"));
        }
    }
    if !finite(material.mdr.psi_b) || !(0.0..=1.0).contains(&material.mdr.psi_b) {
        return Err(invalid("mdr.psi_b", "finite and in [0, 1]"));
    }
    if !finite(material.liquid_bridge.contact_angle)
        || !(0.0..=std::f64::consts::PI).contains(&material.liquid_bridge.contact_angle)
    {
        return Err(invalid(
            "liquid_bridge.contact_angle",
            "finite and in [0, pi]",
        ));
    }
    match material.adhesion {
        Adhesion::None => {}
        Adhesion::Sjkr { energy } if !nonnegative(energy) => {
            return Err(invalid("adhesion.sjkr.energy", "finite and >= 0"));
        }
        Adhesion::SurfaceEnergy { energy } if !nonnegative(energy) => {
            return Err(invalid("adhesion.surface_energy", "finite and >= 0"));
        }
        _ => {}
    }
    for (property, value) in match material.rolling {
        Rolling::Constant => Vec::new(),
        Rolling::Sds { stiffness, damping } => vec![
            ("rolling.sds.stiffness", stiffness),
            ("rolling.sds.damping", damping),
        ],
    }
    .into_iter()
    .chain(match material.twisting {
        Twisting::Constant => Vec::new(),
        Twisting::Sds { stiffness, damping } => vec![
            ("twisting.sds.stiffness", stiffness),
            ("twisting.sds.damping", damping),
        ],
    }) {
        if !nonnegative(value) {
            return Err(invalid(property, "finite and >= 0"));
        }
    }
    Ok(())
}

impl std::fmt::Display for MaterialError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ConflictingCohesion { name } => write!(
                f,
                "material '{}' has both cohesion_energy and surface_energy > 0; use only one",
                name
            ),
            Self::InvalidProperty {
                name,
                property,
                requirement,
            } => write!(
                f,
                "material '{}' property '{}' must be {}",
                name, property, requirement
            ),
        }
    }
}

impl std::error::Error for MaterialError {}

impl Default for MaterialTable {
    fn default() -> Self {
        Self::new()
    }
}

#[allow(deprecated)] // compatibility wrappers deliberately delegate down the old ladder
impl MaterialTable {
    /// Creates an empty `MaterialTable` with default contact/adhesion/rolling/twisting models.
    ///
    /// # Building a table by hand
    ///
    /// The full two-phase pattern — register materials, then build the pair
    /// tables once before any contact force is evaluated:
    ///
    /// ```
    /// use dirt_atom::{Elastic, Friction, Material, MaterialTable};
    ///
    /// let mut mat = MaterialTable::new();
    ///
    /// // Phase 1 — register materials. Configuration input is fallible.
    /// let glass = mat.add(
    ///     Material::new("glass", Elastic::new(8.7e9, 0.3, 0.95))
    ///         .with_friction(Friction { sliding: 0.5, ..Friction::default() }),
    /// ).expect("valid glass material");
    /// assert_eq!(glass.raw(), 0);
    /// assert!(mat.beta_ij.is_empty()); // pair tables still empty in phase 1
    ///
    /// // Phase 2 — build the per-pair mixing tables. Required before contact eval.
    /// mat.build_pair_tables();
    ///
    /// // restitution 0.95 inverts to a small (but non-zero) Hertz damping ratio.
    /// let beta = mat.beta_ij[glass.index()][glass.index()];
    /// assert!(beta > 0.0 && beta < 0.1);
    /// ```
    pub fn new() -> Self {
        MaterialTable {
            names: Vec::new(),
            youngs_mod: Vec::new(),
            poisson_ratio: Vec::new(),
            friction: Vec::new(),
            restitution: Vec::new(),
            rolling_friction: Vec::new(),
            twisting_friction: Vec::new(),
            cohesion_energy: Vec::new(),
            surface_energy: Vec::new(),
            beta_ij: Vec::new(),
            friction_ij: Vec::new(),
            rolling_friction_ij: Vec::new(),
            cohesion_energy_ij: Vec::new(),
            surface_energy_ij: Vec::new(),
            e_eff_ij: Vec::new(),
            g_eff_ij: Vec::new(),
            twisting_friction_ij: Vec::new(),
            kn: Vec::new(),
            kt: Vec::new(),
            kn_ij: Vec::new(),
            kt_ij: Vec::new(),
            contact_model: "hertz".to_string(),
            adhesion_model: "jkr".to_string(),
            rolling_model: "constant".to_string(),
            twisting_model: "constant".to_string(),
            tangential_model: "history".to_string(),
            track_orientation: false,
            limit_damping: true,
            rolling_stiffness: Vec::new(),
            rolling_damping: Vec::new(),
            twisting_stiffness: Vec::new(),
            twisting_damping: Vec::new(),
            mdr_yield_stress: Vec::new(),
            mdr_psi_b: Vec::new(),
            mdr_damping: Vec::new(),
            liquid_bridge_volume: Vec::new(),
            liquid_surface_tension: Vec::new(),
            liquid_contact_angle: Vec::new(),
            liquid_rupture_distance: Vec::new(),
            rolling_stiffness_ij: Vec::new(),
            rolling_damping_ij: Vec::new(),
            twisting_stiffness_ij: Vec::new(),
            twisting_damping_ij: Vec::new(),
            mdr_yield_stress_ij: Vec::new(),
            mdr_psi_b_ij: Vec::new(),
            mdr_damping_ij: Vec::new(),
            liquid_bridge_volume_ij: Vec::new(),
            liquid_surface_tension_ij: Vec::new(),
            liquid_contact_angle_ij: Vec::new(),
            liquid_rupture_distance_ij: Vec::new(),
            liquid_bridge_model: "off".to_string(),
        }
    }

    /// Adds a material with basic properties. Returns the material index or a
    /// [`MaterialError`].
    ///
    /// This is a convenience wrapper around [`add_material_full`](Self::add_material_full)
    /// with `surface_energy = 0.0` (no adhesion).
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn add_material(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
    ) -> Result<u32, MaterialError> {
        self.add_material_full(
            name,
            youngs_mod,
            poisson_ratio,
            restitution,
            friction,
            rolling_friction,
            cohesion_energy,
            0.0,
        )
    }

    /// Adds a material with basic properties plus surface energy. Returns the material index or a
    /// [`MaterialError`].
    ///
    /// Wraps [`add_material_extended`](Self::add_material_extended) with twisting/Hooke stiffness = 0.
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn add_material_full(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
        surface_energy: f64,
    ) -> Result<u32, MaterialError> {
        self.add_material_extended(
            name,
            youngs_mod,
            poisson_ratio,
            restitution,
            friction,
            rolling_friction,
            cohesion_energy,
            surface_energy,
            0.0,
            0.0,
            0.0,
        )
    }

    /// Adds a material with twisting friction and Hooke stiffnesses. Returns the material index or a
    /// [`MaterialError`].
    ///
    /// Wraps [`add_material_with_sds`](Self::add_material_with_sds) with SDS parameters = 0.
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn add_material_extended(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
        surface_energy: f64,
        twisting_friction: f64,
        kn: f64,
        kt: f64,
    ) -> Result<u32, MaterialError> {
        self.add_material_with_sds(
            name,
            youngs_mod,
            poisson_ratio,
            restitution,
            friction,
            rolling_friction,
            cohesion_energy,
            surface_energy,
            twisting_friction,
            kn,
            kt,
            0.0,
            0.0,
            0.0,
            0.0,
        )
    }

    /// Add a material with all fields including SDS rolling/twisting parameters.
    /// Returns [`MaterialError`] for incompatible material options.
    #[allow(clippy::too_many_arguments)]
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn add_material_with_sds(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
        surface_energy: f64,
        twisting_friction: f64,
        kn: f64,
        kt: f64,
        rolling_stiffness: f64,
        rolling_damping: f64,
        twisting_stiffness: f64,
        twisting_damping: f64,
    ) -> Result<u32, MaterialError> {
        self.add_material_with_mdr(
            name,
            youngs_mod,
            poisson_ratio,
            restitution,
            friction,
            rolling_friction,
            cohesion_energy,
            surface_energy,
            twisting_friction,
            kn,
            kt,
            rolling_stiffness,
            rolling_damping,
            twisting_stiffness,
            twisting_damping,
            0.0,
            0.0,
            0.0,
        )
    }

    /// Add a material with all fields including SDS and MDR parameters.
    /// Returns [`MaterialError`] for incompatible material options.
    #[allow(clippy::too_many_arguments)]
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn add_material_with_mdr(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
        surface_energy: f64,
        twisting_friction: f64,
        kn: f64,
        kt: f64,
        rolling_stiffness: f64,
        rolling_damping: f64,
        twisting_stiffness: f64,
        twisting_damping: f64,
        mdr_yield_stress: f64,
        mdr_psi_b: f64,
        mdr_damping: f64,
    ) -> Result<u32, MaterialError> {
        self.add_material_with_liquid_bridge(
            name,
            youngs_mod,
            poisson_ratio,
            restitution,
            friction,
            rolling_friction,
            cohesion_energy,
            surface_energy,
            twisting_friction,
            kn,
            kt,
            rolling_stiffness,
            rolling_damping,
            twisting_stiffness,
            twisting_damping,
            mdr_yield_stress,
            mdr_psi_b,
            mdr_damping,
            0.0,
            0.0,
            0.0,
            0.0,
        )
    }

    /// Add a material with all fields including SDS, MDR, and liquid bridge parameters.
    /// Returns [`MaterialError`] for incompatible material options.
    #[allow(clippy::too_many_arguments)]
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn add_material_with_liquid_bridge(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
        surface_energy: f64,
        twisting_friction: f64,
        kn: f64,
        kt: f64,
        rolling_stiffness: f64,
        rolling_damping: f64,
        twisting_stiffness: f64,
        twisting_damping: f64,
        mdr_yield_stress: f64,
        mdr_psi_b: f64,
        mdr_damping: f64,
        liquid_bridge_volume: f64,
        liquid_surface_tension: f64,
        liquid_contact_angle: f64,
        liquid_rupture_distance: f64,
    ) -> Result<u32, MaterialError> {
        self.try_add_material_with_liquid_bridge(
            name,
            youngs_mod,
            poisson_ratio,
            restitution,
            friction,
            rolling_friction,
            cohesion_energy,
            surface_energy,
            twisting_friction,
            kn,
            kt,
            rolling_stiffness,
            rolling_damping,
            twisting_stiffness,
            twisting_damping,
            mdr_yield_stress,
            mdr_psi_b,
            mdr_damping,
            liquid_bridge_volume,
            liquid_surface_tension,
            liquid_contact_angle,
            liquid_rupture_distance,
        )
    }

    /// Backward-compatible alias for the fallible full material constructor.
    #[allow(clippy::too_many_arguments)]
    #[deprecated(note = "use MaterialTable::add(Material::new(...)) for named material properties")]
    pub fn try_add_material_with_liquid_bridge(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
        surface_energy: f64,
        twisting_friction: f64,
        kn: f64,
        kt: f64,
        rolling_stiffness: f64,
        rolling_damping: f64,
        twisting_stiffness: f64,
        twisting_damping: f64,
        mdr_yield_stress: f64,
        mdr_psi_b: f64,
        mdr_damping: f64,
        liquid_bridge_volume: f64,
        liquid_surface_tension: f64,
        liquid_contact_angle: f64,
        liquid_rupture_distance: f64,
    ) -> Result<u32, MaterialError> {
        if cohesion_energy > 0.0 && surface_energy > 0.0 {
            return Err(MaterialError::ConflictingCohesion {
                name: name.to_string(),
            });
        }
        self.add(
            Material::new(
                name,
                Elastic::new(youngs_mod, poisson_ratio, restitution).with_hooke_stiffness(kn, kt),
            )
            .with_friction(Friction {
                sliding: friction,
                rolling: rolling_friction,
                twisting: twisting_friction,
            })
            .with_adhesion(if cohesion_energy > 0.0 {
                Adhesion::Sjkr {
                    energy: cohesion_energy,
                }
            } else if surface_energy > 0.0 {
                Adhesion::SurfaceEnergy {
                    energy: surface_energy,
                }
            } else {
                Adhesion::None
            })
            .with_rolling(Rolling::Sds {
                stiffness: rolling_stiffness,
                damping: rolling_damping,
            })
            .with_twisting(Twisting::Sds {
                stiffness: twisting_stiffness,
                damping: twisting_damping,
            })
            .with_mdr(Mdr {
                yield_stress: mdr_yield_stress,
                psi_b: mdr_psi_b,
                damping: mdr_damping,
            })
            .with_liquid_bridge(LiquidBridge {
                volume: liquid_bridge_volume,
                surface_tension: liquid_surface_tension,
                contact_angle: liquid_contact_angle,
                rupture_distance: liquid_rupture_distance,
            }),
        )
        .map(MaterialId::raw)
    }

    /// Validates and appends a named material input, returning its typed ID.
    ///
    /// This is the single registration path; it updates every per-material
    /// column together, so a newly added property cannot drift out of alignment
    /// with the pair-table index used by the contact hot path.
    pub fn add(&mut self, material: Material) -> Result<MaterialId, MaterialError> {
        validate_material(&material)?;
        let (cohesion_energy, surface_energy) = match material.adhesion {
            Adhesion::None => (0.0, 0.0),
            Adhesion::Sjkr { energy } => (energy, 0.0),
            Adhesion::SurfaceEnergy { energy } => (0.0, energy),
        };
        let (rolling_stiffness, rolling_damping) = match material.rolling {
            Rolling::Constant => (0.0, 0.0),
            Rolling::Sds { stiffness, damping } => (stiffness, damping),
        };
        let (twisting_stiffness, twisting_damping) = match material.twisting {
            Twisting::Constant => (0.0, 0.0),
            Twisting::Sds { stiffness, damping } => (stiffness, damping),
        };
        let id = MaterialId(self.names.len() as u32);
        self.names.push(material.name);
        self.youngs_mod.push(material.elastic.youngs_mod);
        self.poisson_ratio.push(material.elastic.poisson_ratio);
        self.restitution.push(material.elastic.restitution);
        self.friction.push(material.friction.sliding);
        self.rolling_friction.push(material.friction.rolling);
        self.twisting_friction.push(material.friction.twisting);
        self.cohesion_energy.push(cohesion_energy);
        self.surface_energy.push(surface_energy);
        self.kn.push(material.elastic.normal_stiffness);
        self.kt.push(material.elastic.tangential_stiffness);
        self.rolling_stiffness.push(rolling_stiffness);
        self.rolling_damping.push(rolling_damping);
        self.twisting_stiffness.push(twisting_stiffness);
        self.twisting_damping.push(twisting_damping);
        self.mdr_yield_stress.push(material.mdr.yield_stress);
        self.mdr_psi_b.push(material.mdr.psi_b);
        self.mdr_damping.push(material.mdr.damping);
        self.liquid_bridge_volume
            .push(material.liquid_bridge.volume);
        self.liquid_surface_tension
            .push(material.liquid_bridge.surface_tension);
        self.liquid_contact_angle
            .push(material.liquid_bridge.contact_angle);
        self.liquid_rupture_distance
            .push(material.liquid_bridge.rupture_distance);
        Ok(id)
    }

    /// Looks up a material by name, returning its index if found.
    pub fn find_material(&self, name: &str) -> Option<u32> {
        self.names.iter().position(|n| n == name).map(|i| i as u32)
    }

    /// Extra per-atom neighbor cutoff needed for opt-in liquid bridges.
    ///
    /// SOIL owns only a generic per-atom `cutoff_radius`; DIRT must expand that
    /// radius when a DEM force law acts across a finite gap. Returning the full
    /// maximum range for this material is conservative: for any pair, at least
    /// one side contributes enough padding for the configured bridge rupture
    /// distance, while dry/off materials keep the historical radius-only cutoff.
    pub fn liquid_bridge_cutoff_padding(&self, material_idx: u32) -> f64 {
        if self.liquid_bridge_model != "willett2000" {
            return 0.0;
        }
        let i = material_idx as usize;
        if i >= self.names.len()
            || i >= self.liquid_bridge_volume_ij.len()
            || i >= self.liquid_surface_tension_ij.len()
            || i >= self.liquid_contact_angle_ij.len()
            || i >= self.liquid_rupture_distance_ij.len()
        {
            return 0.0;
        }
        let mut padding = 0.0_f64;
        for j in 0..self.names.len() {
            let volume = self.liquid_bridge_volume_ij[i][j];
            let gamma = self.liquid_surface_tension_ij[i][j];
            if volume <= 0.0 || gamma <= 0.0 {
                continue;
            }
            let theta = self.liquid_contact_angle_ij[i][j];
            let rupture = if self.liquid_rupture_distance_ij[i][j] > 0.0 {
                self.liquid_rupture_distance_ij[i][j]
            } else {
                (1.0 + 0.5 * theta) * volume.cbrt()
            };
            padding = padding.max(rupture.max(0.0));
        }
        padding
    }

    /// Computes all per-pair mixing tables from the registered per-material properties.
    ///
    /// Must be called after all materials have been added. Populates `*_ij` fields using
    /// geometric-mean or harmonic-mean mixing rules as appropriate for each property.
    pub fn build_pair_tables(&mut self) {
        let n = self.names.len();
        self.beta_ij = vec![vec![0.0; n]; n];
        self.friction_ij = vec![vec![0.0; n]; n];
        self.rolling_friction_ij = vec![vec![0.0; n]; n];
        self.cohesion_energy_ij = vec![vec![0.0; n]; n];
        self.surface_energy_ij = vec![vec![0.0; n]; n];
        self.e_eff_ij = vec![vec![0.0; n]; n];
        self.g_eff_ij = vec![vec![0.0; n]; n];
        self.twisting_friction_ij = vec![vec![0.0; n]; n];
        self.kn_ij = vec![vec![0.0; n]; n];
        self.kt_ij = vec![vec![0.0; n]; n];
        self.rolling_stiffness_ij = vec![vec![0.0; n]; n];
        self.rolling_damping_ij = vec![vec![0.0; n]; n];
        self.twisting_stiffness_ij = vec![vec![0.0; n]; n];
        self.twisting_damping_ij = vec![vec![0.0; n]; n];
        self.mdr_yield_stress_ij = vec![vec![0.0; n]; n];
        self.mdr_psi_b_ij = vec![vec![0.0; n]; n];
        self.mdr_damping_ij = vec![vec![0.0; n]; n];
        self.liquid_bridge_volume_ij = vec![vec![0.0; n]; n];
        self.liquid_surface_tension_ij = vec![vec![0.0; n]; n];
        self.liquid_contact_angle_ij = vec![vec![0.0; n]; n];
        self.liquid_rupture_distance_ij = vec![vec![0.0; n]; n];
        // Pad optional fields if old API was used
        while self.surface_energy.len() < n {
            self.surface_energy.push(0.0);
        }
        while self.twisting_friction.len() < n {
            self.twisting_friction.push(0.0);
        }
        while self.kn.len() < n {
            self.kn.push(0.0);
        }
        while self.kt.len() < n {
            self.kt.push(0.0);
        }
        while self.rolling_stiffness.len() < n {
            self.rolling_stiffness.push(0.0);
        }
        while self.rolling_damping.len() < n {
            self.rolling_damping.push(0.0);
        }
        while self.twisting_stiffness.len() < n {
            self.twisting_stiffness.push(0.0);
        }
        while self.twisting_damping.len() < n {
            self.twisting_damping.push(0.0);
        }
        while self.mdr_yield_stress.len() < n {
            self.mdr_yield_stress.push(0.0);
        }
        while self.mdr_psi_b.len() < n {
            self.mdr_psi_b.push(0.0);
        }
        while self.mdr_damping.len() < n {
            self.mdr_damping.push(0.0);
        }
        while self.liquid_bridge_volume.len() < n {
            self.liquid_bridge_volume.push(0.0);
        }
        while self.liquid_surface_tension.len() < n {
            self.liquid_surface_tension.push(0.0);
        }
        while self.liquid_contact_angle.len() < n {
            self.liquid_contact_angle.push(0.0);
        }
        while self.liquid_rupture_distance.len() < n {
            self.liquid_rupture_distance.push(0.0);
        }
        for i in 0..n {
            for j in 0..n {
                // Geometric mean mixing for restitution
                let e_ij = (self.restitution[i] * self.restitution[j]).sqrt();
                let log_e = e_ij.ln();
                // Damping coefficient β derived from the restitution so DIRT
                // integrates the SAME normal-damping force as the reference DEM
                // code (LAMMPS) for the same nominal restitution input:
                //   - Hooke (linear): β = -ln(e)/√(π²+ln²e), exact for a
                //     constant-stiffness spring-dashpot (matches LAMMPS
                //     `damping coeff_restitution` for the linear contact).
                //   - Hertz (nonlinear): β = α(e)/√5 with the Tsuji (1992)
                //     polynomial α(e) (`hertz_beta_for_cor`), which makes DIRT's
                //     `2β√(5/6)√(Sₙmᵣ)vₙ` damping IDENTICAL to LAMMPS
                //     `damping tsuji` `α(e)√(mₑFₙ/δ)vₙ`. Both codes then realize
                //     the same velocity-independent COR (bench_hertz_rebound
                //     cross-check). A previous exact-COR-inversion mapping made
                //     realized COR = nominal e, but that DISAGREED with LAMMPS's
                //     Tsuji model by up to 0.067 — the bug this reconciles.
                self.beta_ij[i][j] = if self.contact_model == "hooke" {
                    -log_e / (PI * PI + log_e * log_e).sqrt()
                } else {
                    hertz_beta_for_cor(e_ij)
                };

                // Geometric mean mixing for friction
                self.friction_ij[i][j] = (self.friction[i] * self.friction[j]).sqrt();

                // Geometric mean mixing for rolling friction
                self.rolling_friction_ij[i][j] =
                    (self.rolling_friction[i] * self.rolling_friction[j]).sqrt();

                // Geometric mean mixing for cohesion energy
                self.cohesion_energy_ij[i][j] =
                    (self.cohesion_energy[i] * self.cohesion_energy[j]).sqrt();

                // Geometric mean mixing for surface energy (JKR)
                self.surface_energy_ij[i][j] =
                    (self.surface_energy[i] * self.surface_energy[j]).sqrt();

                // Geometric mean mixing for twisting friction
                self.twisting_friction_ij[i][j] = (self.twisting_friction[i].max(0.0)
                    * self.twisting_friction[j].max(0.0))
                .sqrt();

                // Effective Young's modulus (Hertz)
                let nu_i = self.poisson_ratio[i];
                let nu_j = self.poisson_ratio[j];
                self.e_eff_ij[i][j] = 1.0
                    / ((1.0 - nu_i * nu_i) / self.youngs_mod[i]
                        + (1.0 - nu_j * nu_j) / self.youngs_mod[j]);

                // Effective shear modulus (Mindlin)
                self.g_eff_ij[i][j] = 1.0
                    / (2.0 * (2.0 - nu_i) * (1.0 + nu_i) / self.youngs_mod[i]
                        + 2.0 * (2.0 - nu_j) * (1.0 + nu_j) / self.youngs_mod[j]);

                // Harmonic mean mixing for Hooke stiffness
                let ki = self.kn[i];
                let kj = self.kn[j];
                self.kn_ij[i][j] = if ki > 0.0 && kj > 0.0 {
                    2.0 * ki * kj / (ki + kj)
                } else {
                    0.0
                };
                let kti = self.kt[i];
                let ktj = self.kt[j];
                self.kt_ij[i][j] = if kti > 0.0 && ktj > 0.0 {
                    2.0 * kti * ktj / (kti + ktj)
                } else {
                    0.0
                };

                // SDS rolling stiffness (harmonic mean)
                let kri = self.rolling_stiffness[i];
                let krj = self.rolling_stiffness[j];
                self.rolling_stiffness_ij[i][j] = if kri > 0.0 && krj > 0.0 {
                    2.0 * kri * krj / (kri + krj)
                } else if kri > 0.0 {
                    kri
                } else {
                    krj
                };

                // SDS rolling damping (geometric mean)
                self.rolling_damping_ij[i][j] =
                    (self.rolling_damping[i].max(0.0) * self.rolling_damping[j].max(0.0)).sqrt();

                // SDS twisting stiffness (harmonic mean)
                let kwi = self.twisting_stiffness[i];
                let kwj = self.twisting_stiffness[j];
                self.twisting_stiffness_ij[i][j] = if kwi > 0.0 && kwj > 0.0 {
                    2.0 * kwi * kwj / (kwi + kwj)
                } else if kwi > 0.0 {
                    kwi
                } else {
                    kwj
                };

                // SDS twisting damping (geometric mean)
                self.twisting_damping_ij[i][j] =
                    (self.twisting_damping[i].max(0.0) * self.twisting_damping[j].max(0.0)).sqrt();

                self.mdr_yield_stress_ij[i][j] =
                    (self.mdr_yield_stress[i].max(0.0) * self.mdr_yield_stress[j].max(0.0)).sqrt();
                self.mdr_psi_b_ij[i][j] = 0.5 * (self.mdr_psi_b[i] + self.mdr_psi_b[j]);
                self.mdr_damping_ij[i][j] =
                    (self.mdr_damping[i].max(0.0) * self.mdr_damping[j].max(0.0)).sqrt();
                self.liquid_bridge_volume_ij[i][j] = (self.liquid_bridge_volume[i].max(0.0)
                    * self.liquid_bridge_volume[j].max(0.0))
                .sqrt();
                self.liquid_surface_tension_ij[i][j] = (self.liquid_surface_tension[i].max(0.0)
                    * self.liquid_surface_tension[j].max(0.0))
                .sqrt();
                self.liquid_contact_angle_ij[i][j] =
                    0.5 * (self.liquid_contact_angle[i] + self.liquid_contact_angle[j]);
                self.liquid_rupture_distance_ij[i][j] = (self.liquid_rupture_distance[i].max(0.0)
                    * self.liquid_rupture_distance[j].max(0.0))
                .sqrt();
            }
        }
    }
}

// ── DemAtom per-atom data ────────────────────────────────────────────────────

/// Per-atom DEM extension data: particle radius, density, inverse inertia, and rotational fields.
///
/// Registered via [`register_atom_data!`] in [`DemAtomPlugin::build`]. The `AtomData` derive
/// macro generates pack/unpack methods for MPI communication based on field attributes:
///
/// - **`#[forward]`** — Sent from owner to ghost atoms each timestep (e.g., `omega`).
/// - **`#[reverse]`** — Accumulated from ghosts back to owner (e.g., `torque`).
/// - **`#[zero]`** — Zeroed before each force computation (e.g., `torque`).
///
/// Fields without attributes are only communicated during atom migration (ownership transfer).
#[derive(AtomData)]
pub struct DemAtom {
    /// Particle radius (m). Set at insertion time; used for contact detection and force calculation.
    pub radius: Vec<f64>,
    /// Particle material density (kg/m³). Used with radius to compute mass and moment of inertia.
    pub density: Vec<f64>,
    /// Inverse moment of inertia (1/(I) where I = 2/5 * m * r²) (1/(kg·m²)).
    /// Precomputed at insertion for efficient torque-to-angular-acceleration conversion.
    pub inv_inertia: Vec<f64>,
    /// Orientation quaternion [w, x, y, z] (unit quaternion). Initialized to [1, 0, 0, 0].
    pub quaternion: Vec<[f64; 4]>,
    /// Angular velocity (rad/s) in [x, y, z] components.
    /// Marked `#[forward]`: communicated from owner to ghost atoms each timestep.
    #[forward]
    pub omega: Vec<[f64; 3]>,
    /// Angular momentum (kg·m²/s) in [x, y, z] components.
    pub ang_mom: Vec<[f64; 3]>,
    /// Torque (N·m) in [x, y, z] components.
    /// Marked `#[reverse]`: ghost contributions are summed back to the owning atom.
    /// Marked `#[zero]`: zeroed before each force computation step.
    #[reverse]
    #[zero]
    pub torque: Vec<[f64; 3]>,
    /// Rigid body ID for clump/multisphere membership. 0.0 = independent particle.
    /// Positive values indicate sub-spheres of the same rigid body (same value = same body).
    #[forward]
    pub body_id: Vec<f64>,
}

impl Default for DemAtom {
    fn default() -> Self {
        Self::new()
    }
}

impl DemAtom {
    /// Creates an empty `DemAtom` with no particles. Particle data is appended during insertion.
    pub fn new() -> Self {
        DemAtom {
            radius: Vec::new(),
            density: Vec::new(),
            inv_inertia: Vec::new(),
            quaternion: Vec::new(),
            omega: Vec::new(),
            ang_mom: Vec::new(),
            torque: Vec::new(),
            body_id: Vec::new(),
        }
    }
}

/// Returns true if atoms `i` and `j` belong to the same rigid body.
#[inline]
pub fn same_body(dem: &DemAtom, i: usize, j: usize) -> bool {
    let bi = dem.body_id[i];
    let bj = dem.body_id[j];
    bi > 0.0 && bj > 0.0 && (bi - bj).abs() < 0.5
}

// ── Plugin ───────────────────────────────────────────────────────────────────

/// Registers [`DemAtom`] extension and [`MaterialTable`] from `[[dem.materials]]` config.
pub struct DemAtomPlugin;

impl Plugin for DemAtomPlugin {
    fn provides(&self) -> Vec<&str> {
        vec!["dem_particles"]
    }

    fn default_config(&self) -> Option<&str> {
        Some(
            r#"# Material definitions for DEM particles
[[dem.materials]]
name = "glass"
youngs_mod = 8.7e9
poisson_ratio = 0.3
restitution = 0.95
friction = 0.4
# rolling_friction = 0.1      # rolling resistance coefficient (default 0.0 = disabled)
# cohesion_energy = 0.05       # SJKR cohesion energy density J/m² (default 0.0 = disabled)
# surface_energy = 0.05        # JKR/DMT surface energy J/m² (default 0.0 = disabled)
# liquid_bridge_volume = 1e-12 # m^3/contact, with liquid_bridge_model = "willett2000"
# liquid_surface_tension = 0.072 # N/m
# liquid_contact_angle = 0.0   # radians
# liquid_rupture_distance = 0.0 # m; 0 = volume-scaled estimate
# mdr_yield_stress = 1.0e5     # Pa, for contact_model = "mdr"
# mdr_psi_b = 0.5              # parsed for LAMMPS MDR parity; bulk branch not yet active
# mdr_damping = 0.0            # MDR normal damping prefactor
# adhesion_model = "jkr"       # "jkr" (default) or "dmt" when surface_energy > 0

# Additional materials can be added:
# [[dem.materials]]
# name = "steel"
# youngs_mod = 200e9
# poisson_ratio = 0.28
# restitution = 0.8
# friction = 0.3"#,
        )
    }

    fn build(&self, app: &mut App) {
        configure_dem_atom(app)
            .unwrap_or_else(|error| panic!("DemAtomPlugin failed to build: {error}"));
    }

    fn try_build(&self, app: &mut App) -> Result<(), AppError> {
        configure_dem_atom(app)
    }
}

fn configure_dem_atom(app: &mut App) -> Result<(), AppError> {
    app.add_plugins(AtomPlugin);

    register_atom_data!(app, DemAtom::new());

    let dem_config = Config::try_load::<DemConfig>(app, "dem")
        .map_err(|error| AppError::message(error.to_string()))?;

    // Ergonomics guard: contact_model = "hooke" silently ignores surface_energy
    // (JKR/DMT adhesion lives only on the Hertz path). Warn loudly instead of
    // silently dropping adhesion. See `hooke_surface_energy_warning`.
    if let Some(msg) = hooke_surface_energy_warning(&dem_config) {
        eprintln!("{}", msg);
    }
    if let Some(msg) = liquid_bridge_model_error(&dem_config) {
        return Err(AppError::message(msg));
    }

    // Build MaterialTable from config at plugin build time
    let mut material_table = MaterialTable::new();

    material_table.contact_model = dem_config.contact_model.clone();
    material_table.adhesion_model = dem_config.adhesion_model.clone();
    material_table.rolling_model = dem_config.rolling_model.clone();
    material_table.twisting_model = dem_config.twisting_model.clone();
    material_table.tangential_model = dem_config.tangential_model.clone();
    material_table.track_orientation = dem_config.track_orientation;
    material_table.limit_damping = dem_config.limit_damping;
    material_table.liquid_bridge_model = dem_config.liquid_bridge_model.clone();

    if let Some(ref materials) = dem_config.materials {
        for mat in materials {
            material_table
                .add(
                    Material::from_config(mat)
                        .map_err(|error| AppError::message(error.to_string()))?,
                )
                .map_err(|error| AppError::message(error.to_string()))?;
        }
        material_table.build_pair_tables();
    }

    app.add_resource(material_table);
    app.add_setup_system(set_dem_ntypes, ScheduleSetupSet::Setup);
    Ok(())
}

/// Setup system that sets `Atom::ntypes` from the number of registered materials.
fn set_dem_ntypes(mut atoms: ResMut<Atom>, material_table: Res<MaterialTable>) {
    if !material_table.names.is_empty() {
        atoms.ntypes = material_table.names.len();
    }
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_material_beta_and_friction() {
        let mut mt = MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0)
            .expect("valid glass material");
        mt.build_pair_tables();

        let e = 0.95_f64;
        // default contact_model is "hertz" → β = α(e)/√5 (Tsuji), matching LAMMPS
        // `damping tsuji` so the two codes integrate the same damping force.
        let expected_beta = hertz_beta_for_cor(e);
        assert!(
            (mt.beta_ij[0][0] - expected_beta).abs() < 1e-12,
            "beta should be {}, got {}",
            expected_beta,
            mt.beta_ij[0][0]
        );
        // β must equal the Tsuji polynomial over √5.
        assert!(
            (mt.beta_ij[0][0] - tsuji_alpha(e) / 5.0_f64.sqrt()).abs() < 1e-12,
            "beta should be α(e)/√5, got {}",
            mt.beta_ij[0][0]
        );
        // The realized COR of a Hertz collision with this β reproduces LAMMPS's
        // Tsuji restitution: it sits slightly ABOVE nominal below e≈0.9 (shared
        // polynomial property). At e=0.95 that offset is small (~0.965).
        let realized = hertz_cor_of_beta(mt.beta_ij[0][0]);
        assert!(
            (0.960..=0.970).contains(&realized),
            "realized COR should be ~0.965 (LAMMPS Tsuji), got {}",
            realized
        );
        assert!(
            (mt.friction_ij[0][0] - 0.4).abs() < 1e-12,
            "friction should be 0.4, got {}",
            mt.friction_ij[0][0]
        );
    }

    #[test]
    fn multi_material_mixing_symmetry() {
        let mut mt = MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0)
            .expect("valid glass material");
        mt.add_material("steel", 200e9, 0.28, 0.8, 0.3, 0.0, 0.0)
            .expect("valid steel material");
        mt.build_pair_tables();

        // Symmetry
        assert!(
            (mt.beta_ij[0][1] - mt.beta_ij[1][0]).abs() < 1e-15,
            "beta_ij should be symmetric"
        );
        assert!(
            (mt.friction_ij[0][1] - mt.friction_ij[1][0]).abs() < 1e-15,
            "friction_ij should be symmetric"
        );

        // Geometric mean mixing for friction
        let expected_friction = (0.4_f64 * 0.3).sqrt();
        assert!(
            (mt.friction_ij[0][1] - expected_friction).abs() < 1e-12,
            "friction_ij should be geometric mean {}, got {}",
            expected_friction,
            mt.friction_ij[0][1]
        );

        // Geometric mean mixing for restitution -> beta (hertz default → Tsuji α/√5)
        let e_mix = (0.95_f64 * 0.8).sqrt();
        let expected_beta = hertz_beta_for_cor(e_mix);
        assert!(
            (mt.beta_ij[0][1] - expected_beta).abs() < 1e-12,
            "beta_ij should use geometric mean restitution"
        );

        // e_eff and g_eff symmetry
        assert!(
            (mt.e_eff_ij[0][1] - mt.e_eff_ij[1][0]).abs() < 1e-6,
            "e_eff_ij should be symmetric"
        );
        assert!(
            (mt.g_eff_ij[0][1] - mt.g_eff_ij[1][0]).abs() < 1e-6,
            "g_eff_ij should be symmetric"
        );
        assert!(mt.e_eff_ij[0][0] > 0.0, "e_eff should be positive");
        assert!(mt.g_eff_ij[0][0] > 0.0, "g_eff should be positive");
    }

    #[test]
    fn typed_materials_match_independent_legacy_mixing_rules() {
        let mut typed = MaterialTable::new();
        let soft = typed
            .add(
                Material::new(
                    "soft",
                    Elastic::new(8.7e9, 0.30, 0.95).with_hooke_stiffness(1.0e6, 5.0e5),
                )
                .with_friction(Friction {
                    sliding: 0.4,
                    rolling: 0.1,
                    twisting: 0.05,
                })
                .with_adhesion(Adhesion::SurfaceEnergy { energy: 0.2 })
                .with_rolling(Rolling::Sds {
                    stiffness: 2.0,
                    damping: 0.3,
                })
                .with_twisting(Twisting::Sds {
                    stiffness: 3.0,
                    damping: 0.4,
                })
                .with_mdr(Mdr {
                    yield_stress: 1.0e6,
                    psi_b: 0.1,
                    damping: 0.02,
                })
                .with_liquid_bridge(LiquidBridge {
                    volume: 1.0e-11,
                    surface_tension: 0.072,
                    contact_angle: 0.2,
                    rupture_distance: 1.0e-4,
                }),
            )
            .unwrap();
        let stiff = typed
            .add(
                Material::new(
                    "stiff",
                    Elastic::new(70e9, 0.22, 0.80).with_hooke_stiffness(2.0e6, 7.0e5),
                )
                .with_friction(Friction {
                    sliding: 0.3,
                    rolling: 0.2,
                    twisting: 0.07,
                })
                .with_adhesion(Adhesion::SurfaceEnergy { energy: 0.5 })
                .with_rolling(Rolling::Sds {
                    stiffness: 4.0,
                    damping: 0.6,
                })
                .with_twisting(Twisting::Sds {
                    stiffness: 5.0,
                    damping: 0.8,
                })
                .with_mdr(Mdr {
                    yield_stress: 2.0e6,
                    psi_b: 0.2,
                    damping: 0.04,
                })
                .with_liquid_bridge(LiquidBridge {
                    volume: 8.0e-12,
                    surface_tension: 0.060,
                    contact_angle: 0.4,
                    rupture_distance: 2.0e-4,
                }),
            )
            .unwrap();
        assert_eq!((soft.raw(), stiff.raw()), (0, 1));
        typed.build_pair_tables();

        // These formulae were transcribed from the pre-PR (origin/main) table
        // builder. They intentionally do not call a deprecated wrapper or the
        // production mixing helpers, making this an independent regression
        // oracle for every generated pair property.
        let geo = |a: f64, b: f64| (a * b).sqrt();
        let harmonic_or_nonzero = |a: f64, b: f64| {
            if a > 0.0 && b > 0.0 {
                2.0 * a * b / (a + b)
            } else if a > 0.0 {
                a
            } else {
                b
            }
        };
        let hertz_beta = |e: f64| {
            let e = e.clamp(1.0e-3, 0.9999);
            (1.2728 - 4.2783 * e + 11.087 * e.powi(2) - 22.348 * e.powi(3) + 27.467 * e.powi(4)
                - 18.022 * e.powi(5)
                + 4.8218 * e.powi(6))
                / 5.0_f64.sqrt()
        };
        let expected = [
            ("beta", typed.beta_ij[0][1], hertz_beta(geo(0.95, 0.80))),
            ("friction", typed.friction_ij[0][1], geo(0.4, 0.3)),
            (
                "rolling_friction",
                typed.rolling_friction_ij[0][1],
                geo(0.1, 0.2),
            ),
            ("cohesion_energy", typed.cohesion_energy_ij[0][1], 0.0),
            (
                "surface_energy",
                typed.surface_energy_ij[0][1],
                geo(0.2, 0.5),
            ),
            (
                "twisting_friction",
                typed.twisting_friction_ij[0][1],
                geo(0.05, 0.07),
            ),
            (
                "e_eff",
                typed.e_eff_ij[0][1],
                1.0 / ((1.0 - 0.30_f64.powi(2)) / 8.7e9 + (1.0 - 0.22_f64.powi(2)) / 70e9),
            ),
            (
                "g_eff",
                typed.g_eff_ij[0][1],
                1.0 / (2.0 * (2.0 - 0.30) * 1.30 / 8.7e9 + 2.0 * (2.0 - 0.22) * 1.22 / 70e9),
            ),
            ("kn", typed.kn_ij[0][1], harmonic_or_nonzero(1.0e6, 2.0e6)),
            ("kt", typed.kt_ij[0][1], harmonic_or_nonzero(5.0e5, 7.0e5)),
            (
                "rolling_stiffness",
                typed.rolling_stiffness_ij[0][1],
                harmonic_or_nonzero(2.0, 4.0),
            ),
            (
                "rolling_damping",
                typed.rolling_damping_ij[0][1],
                geo(0.3, 0.6),
            ),
            (
                "twisting_stiffness",
                typed.twisting_stiffness_ij[0][1],
                harmonic_or_nonzero(3.0, 5.0),
            ),
            (
                "twisting_damping",
                typed.twisting_damping_ij[0][1],
                geo(0.4, 0.8),
            ),
            (
                "mdr_yield_stress",
                typed.mdr_yield_stress_ij[0][1],
                geo(1.0e6, 2.0e6),
            ),
            ("mdr_psi_b", typed.mdr_psi_b_ij[0][1], 0.15),
            ("mdr_damping", typed.mdr_damping_ij[0][1], geo(0.02, 0.04)),
            (
                "liquid_volume",
                typed.liquid_bridge_volume_ij[0][1],
                geo(1.0e-11, 8.0e-12),
            ),
            (
                "liquid_surface_tension",
                typed.liquid_surface_tension_ij[0][1],
                geo(0.072, 0.060),
            ),
            (
                "liquid_contact_angle",
                typed.liquid_contact_angle_ij[0][1],
                0.3,
            ),
            (
                "liquid_rupture_distance",
                typed.liquid_rupture_distance_ij[0][1],
                geo(1.0e-4, 2.0e-4),
            ),
        ];
        for (name, actual, expected) in expected {
            assert!(
                (actual - expected).abs() <= 1.0e-12 * expected.abs().max(1.0),
                "{name}: expected {expected:e}, got {actual:e}"
            );
        }
    }

    #[test]
    fn e_eff_matches_manual_computation() {
        let mut mt = MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0)
            .expect("valid glass material");
        mt.build_pair_tables();

        let nu = 0.3_f64;
        let e = 8.7e9_f64;
        let expected = 1.0 / (2.0 * (1.0 - nu * nu) / e);
        assert!(
            (mt.e_eff_ij[0][0] - expected).abs() < 1.0,
            "e_eff_ij[0][0] should be {}, got {}",
            expected,
            mt.e_eff_ij[0][0]
        );
    }

    #[test]
    fn liquid_bridge_cutoff_padding_tracks_active_rupture_range() {
        let mut mt = MaterialTable::new();
        mt.liquid_bridge_model = "willett2000".to_string();
        let glass = mt
            .add_material_with_liquid_bridge(
                "glass", 1.0e6, 0.25, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0e-11, 0.072, 0.0, 1.5e-4,
            )
            .expect("valid liquid-bridge material");
        mt.add_material("dry", 1.0e6, 0.25, 1.0, 0.0, 0.0, 0.0)
            .expect("valid dry material");
        mt.build_pair_tables();

        assert!((mt.liquid_bridge_cutoff_padding(glass) - 1.5e-4).abs() < 1.0e-15);
        assert_eq!(mt.liquid_bridge_cutoff_padding(1), 0.0);

        mt.liquid_bridge_model = "off".to_string();
        assert_eq!(mt.liquid_bridge_cutoff_padding(glass), 0.0);
    }

    #[test]
    fn conflicting_material_cohesion_is_a_typed_error() {
        let mut table = MaterialTable::new();
        let error = table
            .try_add_material_with_liquid_bridge(
                "bad", 1.0, 0.25, 0.9, 0.4, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            )
            .unwrap_err();
        assert_eq!(
            error,
            MaterialError::ConflictingCohesion {
                name: "bad".to_string()
            }
        );
    }

    #[test]
    fn public_material_constructors_report_conflicting_cohesion() {
        let mut table = MaterialTable::new();
        let error = table
            .add_material_full("bad", 1.0, 0.25, 0.9, 0.4, 0.0, 1.0, 1.0)
            .unwrap_err();
        assert_eq!(
            error,
            MaterialError::ConflictingCohesion {
                name: "bad".to_string()
            }
        );
        assert!(
            table.names.is_empty(),
            "a rejected material must not mutate the table"
        );
    }

    #[test]
    fn typed_material_validation_rejects_invalid_domains_without_mutation() {
        let valid = || Material::new("bad", Elastic::new(1.0e6, 0.25, 0.9));
        let cases = [
            (
                valid().with_friction(Friction {
                    sliding: f64::NAN,
                    ..Friction::default()
                }),
                "friction.sliding",
            ),
            (
                Material::new("bad", Elastic::new(-1.0, 0.25, 0.9)),
                "elastic.youngs_mod",
            ),
            (
                Material::new("bad", Elastic::new(1.0e6, -0.01, 0.9)),
                "elastic.poisson_ratio",
            ),
            (
                Material::new("bad", Elastic::new(1.0e6, 0.25, 1.01)),
                "elastic.restitution",
            ),
            (
                valid().with_rolling(Rolling::Sds {
                    stiffness: -1.0,
                    damping: 0.0,
                }),
                "rolling.sds.stiffness",
            ),
            (
                valid().with_liquid_bridge(LiquidBridge {
                    contact_angle: std::f64::consts::PI + 0.01,
                    ..LiquidBridge::default()
                }),
                "liquid_bridge.contact_angle",
            ),
        ];
        for (material, property) in cases {
            let mut table = MaterialTable::new();
            let error = table.add(material).unwrap_err();
            assert!(
                matches!(error, MaterialError::InvalidProperty { property: actual, .. } if actual == property)
            );
            assert!(
                table.names.is_empty(),
                "invalid input must not mutate any column"
            );
            assert!(table.youngs_mod.is_empty());
            assert!(table.liquid_bridge_volume.is_empty());
        }
    }

    // ── hooke_surface_energy_warning: silent-drop ergonomics guard ────────────

    /// Build a minimal `MaterialConfig` with a given `surface_energy`; all other
    /// force-law parameters are zeroed (they are irrelevant to the diagnostic).
    fn mat(name: &str, surface_energy: f64) -> MaterialConfig {
        MaterialConfig {
            name: name.to_string(),
            youngs_mod: 8.7e9,
            poisson_ratio: 0.3,
            restitution: 0.9,
            friction: 0.4,
            rolling_friction: 0.0,
            cohesion_energy: 0.0,
            surface_energy,
            twisting_friction: 0.0,
            kn: 0.0,
            kt: 0.0,
            rolling_stiffness: 0.0,
            rolling_damping: 0.0,
            twisting_stiffness: 0.0,
            twisting_damping: 0.0,
            mdr_yield_stress: 0.0,
            mdr_psi_b: 0.0,
            mdr_damping: 0.0,
            liquid_bridge_volume: 0.0,
            liquid_surface_tension: 0.0,
            liquid_contact_angle: 0.0,
            liquid_rupture_distance: 0.0,
        }
    }

    fn cfg(contact_model: &str, materials: Vec<MaterialConfig>) -> DemConfig {
        DemConfig {
            contact_model: contact_model.to_string(),
            materials: Some(materials),
            ..DemConfig::default()
        }
    }

    #[test]
    fn warns_for_hooke_with_surface_energy() {
        // hooke + surface_energy > 0 → diagnostic fires, names the field, the
        // offending material, and points at the Hertz path.
        let config = cfg("hooke", vec![mat("glass", 0.05)]);
        let msg = hooke_surface_energy_warning(&config)
            .expect("hooke + surface_energy>0 must produce a diagnostic");
        assert!(
            msg.contains("surface_energy"),
            "must name the ignored field: {msg}"
        );
        assert!(
            msg.contains("glass"),
            "must name the offending material: {msg}"
        );
        assert!(msg.contains("hertz"), "must point at the hertz path: {msg}");
        assert!(
            msg.contains("hooke"),
            "must name the offending contact_model: {msg}"
        );
    }

    #[test]
    fn warns_lists_all_offending_materials() {
        let config = cfg(
            "hooke",
            vec![mat("glass", 0.05), mat("dry", 0.0), mat("wet", 0.02)],
        );
        let msg = hooke_surface_energy_warning(&config).expect("diagnostic must fire");
        assert!(
            msg.contains("glass") && msg.contains("wet"),
            "names all offenders: {msg}"
        );
        assert!(
            !msg.contains("dry"),
            "must not name a zero-surface_energy material: {msg}"
        );
    }

    #[test]
    fn silent_for_hertz_with_surface_energy() {
        // hertz + surface_energy > 0 is the *correct* JKR/DMT configuration — no
        // diagnostic. (Default contact_model is "hertz".)
        let config = cfg("hertz", vec![mat("glass", 0.05)]);
        assert!(
            hooke_surface_energy_warning(&config).is_none(),
            "hertz + surface_energy>0 is valid and must NOT warn"
        );
        let default_model = cfg(&default_contact_model(), vec![mat("glass", 0.05)]);
        assert!(
            hooke_surface_energy_warning(&default_model).is_none(),
            "default (hertz) + surface_energy>0 must NOT warn"
        );
    }

    #[test]
    fn silent_for_hooke_without_surface_energy() {
        // hooke + surface_energy == 0 drops nothing → no diagnostic.
        let config = cfg("hooke", vec![mat("glass", 0.0), mat("steel", 0.0)]);
        assert!(
            hooke_surface_energy_warning(&config).is_none(),
            "hooke + surface_energy=0 drops nothing and must NOT warn"
        );
    }

    #[test]
    fn silent_for_hooke_with_no_materials() {
        let config = cfg("hooke", vec![]);
        assert!(
            hooke_surface_energy_warning(&config).is_none(),
            "hooke with no materials must NOT warn"
        );
    }

    #[test]
    fn accepts_supported_liquid_bridge_model_selectors() {
        let mut config = DemConfig::default();
        for model in ["off", "willett2000"] {
            config.liquid_bridge_model = model.to_string();
            assert!(
                liquid_bridge_model_error(&config).is_none(),
                "{model} must be accepted"
            );
        }
    }

    #[test]
    fn rejects_invalid_liquid_bridge_model_selector() {
        let mut config = DemConfig {
            liquid_bridge_model: "willet2000".to_string(),
            ..DemConfig::default()
        };
        let msg = liquid_bridge_model_error(&config).expect("typo must be rejected");
        assert!(
            msg.contains("liquid_bridge_model"),
            "must name the bad selector field: {msg}"
        );
        assert!(
            msg.contains("willet2000"),
            "must echo the invalid value: {msg}"
        );
        assert!(
            msg.contains("willett2000") && msg.contains("off"),
            "must list supported values: {msg}"
        );

        config.liquid_bridge_model = "Willett2000".to_string();
        assert!(
            liquid_bridge_model_error(&config).is_some(),
            "case-mismatched selectors must not silently disable bridge physics"
        );
    }
}
