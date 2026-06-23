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
//! contact_model = "hertz"  # "hertz" (default) or "hooke"
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

pub mod insert;
pub mod radius;

pub use insert::*;
pub use radius::*;

use std::f64::consts::PI;

use grass_app::prelude::*;
use soil_derive::AtomData;
use grass_scheduler::prelude::*;
use serde::Deserialize;

use soil_core::{register_atom_data, Atom, AtomData, AtomPlugin, Config, ScheduleSetupSet};

// ── Shared physics constants ────────────────────────────────────────────────

pub const SQRT_5_6: f64 = 0.9128709291752768;

// ── Exact coefficient-of-restitution damping (Hertz) ─────────────────────────

/// COR of a head-on Hertz collision with the damping DIRT applies
/// (`f_diss = 2β√(5/6)√(Sₙ mᵣ) vₙ`, `Sₙ = 2E*√(R*δ)`), as a function of β.
///
/// Computed by integrating the dimensionless 1D collision (E*=R*=mᵣ=1, v₀=1):
/// `δ̈ = −(4/3)δ^{3/2} − 2β√(5/6)√2 · δ^{1/4} · δ̇`. The Tsuji scaling makes this
/// velocity-independent, so one integration fixes the β↔COR map. RK4, fixed dt.
fn hertz_cor_of_beta(beta: f64) -> f64 {
    if beta <= 0.0 {
        return 1.0;
    }
    let c = 2.0 * beta * SQRT_5_6 * std::f64::consts::SQRT_2; // damping prefactor
    // acceleration of the overlap coordinate (only while in contact, δ>0)
    let acc = |d: f64, v: f64| -> f64 {
        if d <= 0.0 { 0.0 } else { -(4.0 / 3.0) * d.powf(1.5) - c * d.powf(0.25) * v }
    };
    let dt = 1.0e-4;
    // δ = overlap (grows on approach): start at contact with δ̇ = +1 (approaching).
    let (mut d, mut v) = (0.0_f64, 1.0_f64);
    for _ in 0..2_000_000 {
        // RK4 on (d, v) with d' = v, v' = acc(d, v)
        let (k1d, k1v) = (v, acc(d, v));
        let (k2d, k2v) = (v + 0.5 * dt * k1v, acc(d + 0.5 * dt * k1d, v + 0.5 * dt * k1v));
        let (k3d, k3v) = (v + 0.5 * dt * k2v, acc(d + 0.5 * dt * k2d, v + 0.5 * dt * k2v));
        let (k4d, k4v) = (v + dt * k3v, acc(d + dt * k3d, v + dt * k3v));
        d += dt / 6.0 * (k1d + 2.0 * k2d + 2.0 * k3d + k4d);
        v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        if d <= 0.0 && v < 0.0 {
            return v.abs(); // separated (overlap back to 0, receding); COR = |v_out| (v_in = 1)
        }
    }
    v.abs()
}

/// Invert [`hertz_cor_of_beta`]: the damping ratio β giving an exact restitution
/// `e_target` for the Hertz contact (DIRT's analogue of LAMMPS `damping
/// coeff_restitution`). COR(β) decreases monotonically from 1 (β=0), so bisect.
/// This makes the **input restitution the realized COR**, so DIRT shear/cooling
/// land on the same kinetic-theory / cross-code line as LAMMPS/LIGGGHTS rather
/// than on a shifted "realized-e" line.
pub fn hertz_beta_for_cor(e_target: f64) -> f64 {
    if e_target >= 0.9999 {
        return 0.0;
    }
    let e = e_target.clamp(1.0e-3, 0.9999);
    let (mut lo, mut hi) = (0.0_f64, 5.0_f64); // β range; COR(5) ≈ 0
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if hertz_cor_of_beta(mid) > e {
            lo = mid; // not enough damping → raise β
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
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
}

fn default_adhesion_model() -> String {
    "jkr".to_string()
}

fn default_rolling_model() -> String {
    "constant".to_string()
}

fn default_twisting_model() -> String {
    "constant".to_string()
}

/// TOML `[dem]` — top-level DEM configuration containing material definitions.
#[derive(Deserialize, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct DemConfig {
    /// List of material definitions, each corresponding to a `[[dem.materials]]` block.
    pub materials: Option<Vec<MaterialConfig>>,
    /// Contact model: "hertz" (default) or "hooke".
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
/// `restitution` is stored as the **target coefficient of restitution (COR)**,
/// not a damping ratio. In phase 2, `beta_ij[i][j]` is computed by inverting the
/// exact head-on Hertz collision via [`hertz_beta_for_cor`] — a bisection on the
/// monotone COR(β) curve of the head-on Hertz model. This makes the **input
/// restitution the realized COR** of a binary collision, so DIRT's
/// shear/cooling results land on the same kinetic-theory line as
/// LAMMPS/LIGGGHTS (`damping coeff_restitution`) rather than on a shifted
/// "realized-e" line. A plain damping-ratio mapping would not have that
/// property.
///
/// # Config-error convention
///
/// `add_material*` validates physically inconsistent input (e.g. both
/// `cohesion_energy` and `surface_energy` set) by printing an `ERROR:` line to
/// stderr and calling `std::process::exit(1)` — it does **not** return a
/// `Result`. This is deliberate: a malformed material table is a config bug that
/// should stop the run immediately and identically on every MPI rank, rather
/// than propagate a partially-built table.
///
/// # Hooke vs. Hertz adhesion asymmetry
///
/// JKR/DMT adhesion (`surface_energy`) is only honored under the Hertz contact
/// model. Under `contact_model = "hooke"` the surface-energy term is silently
/// ignored — only SJKR-style cohesion (`cohesion_energy`) is applied. See
/// `dirt_granular` for the per-branch parameter reference.
pub struct MaterialTable {
    pub names: Vec<String>,
    pub youngs_mod: Vec<f64>,
    pub poisson_ratio: Vec<f64>,
    pub friction: Vec<f64>,
    pub restitution: Vec<f64>,
    pub rolling_friction: Vec<f64>,
    pub twisting_friction: Vec<f64>,
    pub cohesion_energy: Vec<f64>,
    pub surface_energy: Vec<f64>,
    pub beta_ij: Vec<Vec<f64>>,
    pub friction_ij: Vec<Vec<f64>>,
    pub rolling_friction_ij: Vec<Vec<f64>>,
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
    /// Track per-sphere orientation (quaternion). Default `false`; see [`DemConfig::track_orientation`].
    pub track_orientation: bool,
    /// Per-material rolling spring stiffness (SDS model).
    pub rolling_stiffness: Vec<f64>,
    /// Per-material rolling damping coefficient (SDS model).
    pub rolling_damping: Vec<f64>,
    /// Per-material twisting spring stiffness (SDS model).
    pub twisting_stiffness: Vec<f64>,
    /// Per-material twisting damping coefficient (SDS model).
    pub twisting_damping: Vec<f64>,
    /// Per-pair rolling stiffness (harmonic mean).
    pub rolling_stiffness_ij: Vec<Vec<f64>>,
    /// Per-pair rolling damping (geometric mean).
    pub rolling_damping_ij: Vec<Vec<f64>>,
    /// Per-pair twisting stiffness (harmonic mean).
    pub twisting_stiffness_ij: Vec<Vec<f64>>,
    /// Per-pair twisting damping (geometric mean).
    pub twisting_damping_ij: Vec<Vec<f64>>,
}

impl Default for MaterialTable {
    fn default() -> Self {
        Self::new()
    }
}

impl MaterialTable {
    /// Creates an empty `MaterialTable` with default contact/adhesion/rolling/twisting models.
    ///
    /// # Building a table by hand
    ///
    /// The full two-phase pattern — register materials, then build the pair
    /// tables once before any contact force is evaluated:
    ///
    /// ```
    /// use dirt_atom::MaterialTable;
    ///
    /// let mut mat = MaterialTable::new();
    ///
    /// // Phase 1 — register materials. Returns the material index.
    /// let glass = mat.add_material(
    ///     "glass",
    ///     8.7e9, // Young's modulus E [Pa]
    ///     0.3,   // Poisson's ratio ν
    ///     0.95,  // restitution (target COR)
    ///     0.5,   // sliding friction
    ///     0.0,   // rolling friction
    ///     0.0,   // cohesion energy
    /// );
    /// assert_eq!(glass, 0);
    /// assert!(mat.beta_ij.is_empty()); // pair tables still empty in phase 1
    ///
    /// // Phase 2 — build the per-pair mixing tables. Required before contact eval.
    /// mat.build_pair_tables();
    ///
    /// // restitution 0.95 inverts to a small (but non-zero) Hertz damping ratio.
    /// let beta = mat.beta_ij[glass as usize][glass as usize];
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
            track_orientation: false,
            rolling_stiffness: Vec::new(),
            rolling_damping: Vec::new(),
            twisting_stiffness: Vec::new(),
            twisting_damping: Vec::new(),
            rolling_stiffness_ij: Vec::new(),
            rolling_damping_ij: Vec::new(),
            twisting_stiffness_ij: Vec::new(),
            twisting_damping_ij: Vec::new(),
        }
    }

    /// Adds a material with basic properties. Returns the material index.
    ///
    /// This is a convenience wrapper around [`add_material_full`](Self::add_material_full)
    /// with `surface_energy = 0.0` (no adhesion).
    pub fn add_material(
        &mut self,
        name: &str,
        youngs_mod: f64,
        poisson_ratio: f64,
        restitution: f64,
        friction: f64,
        rolling_friction: f64,
        cohesion_energy: f64,
    ) -> u32 {
        self.add_material_full(name, youngs_mod, poisson_ratio, restitution, friction, rolling_friction, cohesion_energy, 0.0)
    }

    /// Adds a material with basic properties plus surface energy. Returns the material index.
    ///
    /// Wraps [`add_material_extended`](Self::add_material_extended) with twisting/Hooke stiffness = 0.
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
    ) -> u32 {
        self.add_material_extended(
            name, youngs_mod, poisson_ratio, restitution, friction,
            rolling_friction, cohesion_energy, surface_energy, 0.0, 0.0, 0.0,
        )
    }

    /// Adds a material with twisting friction and Hooke stiffnesses. Returns the material index.
    ///
    /// Wraps [`add_material_with_sds`](Self::add_material_with_sds) with SDS parameters = 0.
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
    ) -> u32 {
        self.add_material_with_sds(
            name, youngs_mod, poisson_ratio, restitution, friction,
            rolling_friction, cohesion_energy, surface_energy, twisting_friction,
            kn, kt, 0.0, 0.0, 0.0, 0.0,
        )
    }

    /// Add a material with all fields including SDS rolling/twisting parameters.
    #[allow(clippy::too_many_arguments)]
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
    ) -> u32 {
        if cohesion_energy > 0.0 && surface_energy > 0.0 {
            eprintln!(
                "ERROR: material '{}' has both cohesion_energy and surface_energy > 0. Use only one.",
                name
            );
            std::process::exit(1);
        }
        let idx = self.names.len() as u32;
        self.names.push(name.to_string());
        self.youngs_mod.push(youngs_mod);
        self.poisson_ratio.push(poisson_ratio);
        self.restitution.push(restitution);
        self.friction.push(friction);
        self.rolling_friction.push(rolling_friction);
        self.twisting_friction.push(twisting_friction);
        self.cohesion_energy.push(cohesion_energy);
        self.surface_energy.push(surface_energy);
        self.kn.push(kn);
        self.kt.push(kt);
        self.rolling_stiffness.push(rolling_stiffness);
        self.rolling_damping.push(rolling_damping);
        self.twisting_stiffness.push(twisting_stiffness);
        self.twisting_damping.push(twisting_damping);
        idx
    }

    /// Looks up a material by name, returning its index if found.
    pub fn find_material(&self, name: &str) -> Option<u32> {
        self.names.iter().position(|n| n == name).map(|i| i as u32)
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
        for i in 0..n {
            for j in 0..n {
                // Geometric mean mixing for restitution
                let e_ij = (self.restitution[i] * self.restitution[j]).sqrt();
                let log_e = e_ij.ln();
                // Damping coefficient β from the restitution, so the realized COR
                // EQUALS the input e (DIRT's analogue of LAMMPS `damping
                // coeff_restitution`). This matters because T*∝1/(1−e²) near the
                // elastic limit, so any input-vs-realized gap throws the stress off
                // and DIRT would land on a separate line from LAMMPS/LIGGGHTS/KT.
                //   - Hooke (linear): β = -ln(e)/√(π²+ln²e) is exact for a
                //     constant-stiffness spring-dashpot (velocity-independent).
                //   - Hertz (nonlinear): the old Tsuji *polynomial* fit realized a
                //     COR above nominal (e.g. 0.95→0.965). Replace it with a
                //     numerically EXACT inversion of the Hertz-collision COR(β)
                //     (`hertz_beta_for_cor`), so input e = realized COR.
                self.beta_ij[i][j] = if self.contact_model == "hertz" {
                    hertz_beta_for_cor(e_ij)
                } else {
                    -log_e / (PI * PI + log_e * log_e).sqrt()
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
                self.twisting_friction_ij[i][j] =
                    (self.twisting_friction[i].max(0.0) * self.twisting_friction[j].max(0.0)).sqrt();

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
        app.add_plugins(AtomPlugin);

        register_atom_data!(app, DemAtom::new());

        let dem_config = Config::load::<DemConfig>(app, "dem");

        // Build MaterialTable from config at plugin build time
        let mut material_table = MaterialTable::new();

        material_table.contact_model = dem_config.contact_model.clone();
        material_table.adhesion_model = dem_config.adhesion_model.clone();
        material_table.rolling_model = dem_config.rolling_model.clone();
        material_table.twisting_model = dem_config.twisting_model.clone();
        material_table.track_orientation = dem_config.track_orientation;

        if let Some(ref materials) = dem_config.materials {
            for mat in materials {
                material_table.add_material_with_sds(
                    &mat.name,
                    mat.youngs_mod,
                    mat.poisson_ratio,
                    mat.restitution,
                    mat.friction,
                    mat.rolling_friction,
                    mat.cohesion_energy,
                    mat.surface_energy,
                    mat.twisting_friction,
                    mat.kn,
                    mat.kt,
                    mat.rolling_stiffness,
                    mat.rolling_damping,
                    mat.twisting_stiffness,
                    mat.twisting_damping,
                );
            }
            material_table.build_pair_tables();
        }

        app.add_resource(material_table);
        app.add_setup_system(set_dem_ntypes, ScheduleSetupSet::Setup);
    }
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
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0);
        mt.build_pair_tables();

        let e = 0.95_f64;
        // default contact_model is "hertz" → β set by the exact COR inversion, so
        // the realized restitution of a Hertz collision equals the input e.
        let expected_beta = hertz_beta_for_cor(e);
        assert!(
            (mt.beta_ij[0][0] - expected_beta).abs() < 1e-12,
            "beta should be {}, got {}",
            expected_beta,
            mt.beta_ij[0][0]
        );
        // The defining invariant: realized COR(β) ≈ nominal e.
        assert!(
            (hertz_cor_of_beta(mt.beta_ij[0][0]) - e).abs() < 2e-3,
            "realized COR should equal input e={}, got {}",
            e,
            hertz_cor_of_beta(mt.beta_ij[0][0])
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
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0);
        mt.add_material("steel", 200e9, 0.28, 0.8, 0.3, 0.0, 0.0);
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

        // Geometric mean mixing for restitution -> beta (hertz default → exact COR)
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
    fn e_eff_matches_manual_computation() {
        let mut mt = MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0);
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
}
