//! Bonded Particle Model (BPM) forces for DEM simulations.
//!
//! This crate provides [`DemBondPlugin`], which adds elastic bond forces
//! between particle pairs using a beam-theory (cylindrical bond) model.
//! Each bond resists four independent deformation channels:
//!
//! ## Deformation channels
//!
//! | Channel                   | Stiffness (beam form) | Physical meaning                             |
//! |---------------------------|-----------------------|----------------------------------------------|
//! | Normal (extension/compression) | `K_n   = E · A / L`   | stretching/compressing the bond along **n̂** |
//! | Shear                     | `K_t   = G · A / L`   | sliding perpendicular to **n̂**              |
//! | Twist (torsion)           | `K_tor = G · J / L`   | rotating about **n̂**                        |
//! | Bending                   | `K_bend = E · I / L`  | relative rotation perpendicular to **n̂**    |
//!
//! where for a solid cylindrical bond of radius `r_b`:
//! `A = π r_b²`, `J = ½ π r_b⁴` (polar second moment), and
//! `I = ¼ π r_b⁴ = ½ J` (second moment for bending).
//! `L = r₀` is the equilibrium bond length.
//!
//! ## Force and moment equations
//!
//! **Normal force** along unit bond axis **n̂** (from *i* to *j*):
//!
//! ```text
//! F_n = (K_n · δ  +  γ_n · v_n) · n̂,     δ = |r_ij| − r₀
//! ```
//!
//! **Shear force** (history-dependent, Δs re-projected ⊥ to current n̂ each step):
//!
//! ```text
//! F_t  = K_t · Δs  +  γ_t · v_t
//! ```
//!
//! applied as **+F_t on atom i (lower tag) and −F_t on atom j** — so when the
//! higher-tag atom slides below the lower-tag atom (v_t · ẑ < 0), the lower-tag
//! atom is pulled downward and the higher-tag atom is pulled back up,
//! damping the relative transverse motion. Shear is evaluated at the
//! bond mid-point; the resulting force produces a lever-arm torque on both
//! particles (`τ_shear = (L/2) n̂ × F_t`).
//!
//! **Twist moment** (along n̂, Δθ component along n̂):
//!
//! ```text
//! M_tor  = K_tor · (Δθ · n̂) n̂  +  γ_tor · (ω_rel · n̂) n̂
//! ```
//!
//! **Bending moment** (⊥ to n̂, Δθ with n̂ component removed):
//!
//! ```text
//! M_bend = K_bend · (Δθ − (Δθ · n̂) n̂)  +  γ_bend · (ω_rel − (ω_rel · n̂) n̂)
//! ```
//!
//! Both moments applied as **+M on atom i, −M on atom j**, which damps
//! relative rotation (matches LIGGGHTS/Fortran BPM convention).
//!
//! ## Damping
//!
//! Per-channel viscous damping is derived from a **critical-damping ratio** `β ∈ [0, 1]`:
//!
//! ```text
//! γ   = 2 β √( m* · K_eff )      for F_n, F_t
//! γ_M = 2 β √( I* · K_eff )      for M_tor, M_bend
//! ```
//!
//! using the reduced mass `m* = m_i m_j / (m_i + m_j)` and reduced MOI
//! `I* = I_i I_j / (I_i + I_j)` of the bonded pair. Each channel accepts an
//! optional raw-`γ` override that bypasses the β-based calculation.
//!
//! ## Breakage (beam-stress criterion)
//!
//! A bond breaks when either combined stress at the extreme fibre exceeds its limit:
//!
//! ```text
//! σ = F_n / A  +  2 |M_bend| r_b / J     →  break if σ > σ_max   (tensile)
//! τ = |F_t| / A  +  |M_tor| r_b / J      →  break if τ > τ_max   (shear)
//! ```
//!
//! ## Bond geometry
//!
//! The bond radius for each pair is `r_b = bond_radius_ratio · min(R_i, R_j)`.
//! Setting `bond_radius_ratio = 1.0` gives bonds as wide as the smaller particle.
//!
//! ## Configuration mode
//!
//! You can parametrise stiffness two ways:
//!
//! - **Material mode:** give `youngs_modulus` *E* and `shear_modulus` *G*;
//!   stiffnesses are derived per-bond from beam theory. This is the paper-standard mode.
//! - **Direct mode:** give `normal_stiffness`, `shear_stiffness`, `twist_stiffness`,
//!   `bending_stiffness` directly (units N/m or N·m/rad). Used when E/G are not set.
//!
//! If both are set, material-mode wins for whichever channels E/G apply to.
//!
//! ## TOML example
//!
//! ```toml
//! [bonds]
//! auto_bond = true
//! bond_tolerance = 1.001
//! bond_radius_ratio = 1.0
//!
//! # Material mode
//! youngs_modulus = 1.0e9      # E (Pa)
//! shear_modulus  = 4.0e8      # G (Pa)
//!
//! # Damping ratios (critical = 1.0)
//! beta_normal  = 0.05
//! beta_shear   = 0.05
//! beta_twist   = 0.05
//! beta_bending = 0.05
//!
//! # Breakage (combined-stress with Weibull tensile threshold):
//! seed = 0                    # Weibull-sampler RNG seed
//! [bonds.breakage]
//! kind    = "combined_stress"
//! tensile = { kind = "weibull", mean = 5.0e7, m = 8.0, l_calib = 0.020 }
//! shear   = { kind = "constant", value = 3.0e7 }
//!
//! # Plasticity — bending and axial channels independently configurable.
//! # Bending: Guo 2018 elastic-perfectly-plastic at M^p = (4/3) σ_0 r_b³.
//! [bonds.plasticity.bending]
//! kind         = "guo_bending"
//! yield_stress = 1.23e8
//!
//! # Axial: piecewise-linear hardening — slope K_e until 1% strain, then
//! # K_e/2 until 2%, then K_e/10 until 3%, then flat (perfectly plastic).
//! [bonds.plasticity.axial]
//! kind                = "piecewise"
//! breakpoint_strains  = [0.01, 0.02, 0.03]
//! slope_multipliers   = [0.5, 0.1, 0.0]
//! ```

use std::any::Any;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufRead, BufReader};

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use serde::Deserialize;

use dirt_atom::DemAtom;
use soil_core::{Atom, AtomData, AtomDataRegistry, BondEntry, BondStore, CommResource, Config, Domain, ParticleSimScheduleSet, ScheduleSetupSet, VirialStress, VirialStressPlugin};
use soil_print::Thermo;

pub mod breakage;
pub mod plasticity;
use breakage::{BreakageConfig, BreakageCriterion, Unbreakable};
use plasticity::{BondPlasticityModel, PlasticityConfig};

// ── BondConfig ──────────────────────────────────────────────────────────────

/// Deserialized TOML `[bonds]` configuration section.
///
/// See the [module-level docs](crate) for the full parameter reference.
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct BondConfig {
    /// When `true`, automatically bond all particle pairs whose centre-to-centre
    /// distance is within `bond_tolerance × (R_i + R_j)` during setup.
    #[serde(default)]
    pub auto_bond: bool,
    /// Multiplier on sum-of-radii for auto-bond eligibility. Default: `1.001`.
    #[serde(default = "default_bond_tolerance")]
    pub bond_tolerance: f64,
    /// Bond radius as a multiple of `min(R_i, R_j)`. Default: `1.0`.
    #[serde(default = "default_bond_radius_ratio")]
    pub bond_radius_ratio: f64,
    /// Multiplier applied to the maximum bond length when extending MPI
    /// `ghost_cutoff` at setup. Must cover bonded-pair reach (1×) plus
    /// shared-neighbour 1-3 exclusion (2×) plus a safety margin for stretch.
    /// Default: `2.5` — enough for 1-3 exclusion + 25 % bond stretch.
    /// Set to `0.0` to disable the extension (single-process / MPI-1×1×1 only).
    #[serde(default = "default_ghost_cutoff_multiplier")]
    pub ghost_cutoff_multiplier: f64,

    // ── Material-mode inputs ────────────────────────────────────────────────
    /// Young's modulus *E* (Pa). If set, normal & bending stiffness derive from `E`.
    #[serde(default)]
    pub youngs_modulus: Option<f64>,
    /// Shear modulus *G* (Pa). If set, shear & twist stiffness derive from `G`.
    #[serde(default)]
    pub shear_modulus: Option<f64>,

    // ── Direct stiffness overrides ──────────────────────────────────────────
    /// Direct normal stiffness (N/m). Used when `youngs_modulus` is not set.
    #[serde(default)]
    pub normal_stiffness: f64,
    /// Direct shear stiffness (N/m). Used when `shear_modulus` is not set.
    #[serde(default)]
    pub shear_stiffness: f64,
    /// Direct twist (torsion) stiffness (N·m/rad). Used when `shear_modulus` is not set.
    #[serde(default)]
    pub twist_stiffness: f64,
    /// Direct bending stiffness (N·m/rad). Used when `youngs_modulus` is not set.
    #[serde(default)]
    pub bending_stiffness: f64,

    // ── Critical-damping ratios ─────────────────────────────────────────────
    /// Normal damping ratio β ∈ [0,1]. Critical damping = 1.0. Default: `0.0`.
    #[serde(default)]
    pub beta_normal: f64,
    /// Shear damping ratio β ∈ [0,1]. Default: `0.0`.
    #[serde(default)]
    pub beta_shear: f64,
    /// Twist damping ratio β ∈ [0,1]. Default: `0.0`.
    #[serde(default)]
    pub beta_twist: f64,
    /// Bending damping ratio β ∈ [0,1]. Default: `0.0`.
    #[serde(default)]
    pub beta_bending: f64,

    // ── Raw damping overrides (optional) ────────────────────────────────────
    /// Raw normal damping *γ_n* (N·s/m). Overrides `beta_normal` when set.
    #[serde(default)]
    pub normal_damping: Option<f64>,
    /// Raw shear damping *γ_t* (N·s/m). Overrides `beta_shear` when set.
    #[serde(default)]
    pub shear_damping: Option<f64>,
    /// Raw twist damping *γ_tor* (N·m·s/rad). Overrides `beta_twist` when set.
    #[serde(default)]
    pub twist_damping: Option<f64>,
    /// Raw bending damping *γ_bend* (N·m·s/rad). Overrides `beta_bending` when set.
    #[serde(default)]
    pub bending_damping: Option<f64>,

    // ── Breakage ────────────────────────────────────────────────────────────
    /// Breakage criterion configuration. See [`breakage::BreakageConfig`] for
    /// the variant menu and TOML syntax. Default: `None` (bonds never break).
    #[serde(default)]
    pub breakage: Option<BreakageConfig>,
    /// Seed for the per-bond Weibull threshold RNG. Sampling is deterministic
    /// given the seed and the bond-creation order. Defaults to `0`.
    #[serde(default)]
    pub seed: u64,

    // ── Plasticity ──────────────────────────────────────────────────────────
    /// Plasticity model configuration. See [`plasticity::PlasticityConfig`].
    /// Default: `None` (every channel is purely elastic).
    #[serde(default)]
    pub plasticity: Option<PlasticityConfig>,

    /// Path to a LAMMPS data file containing a `Bonds` section.
    pub file: Option<String>,
    /// File format identifier. Only `"lammps_data"` is supported.
    pub format: Option<String>,
}

fn default_bond_tolerance() -> f64 { 1.001 }
fn default_bond_radius_ratio() -> f64 { 1.0 }
fn default_ghost_cutoff_multiplier() -> f64 { 2.5 }

impl Default for BondConfig {
    fn default() -> Self {
        BondConfig {
            auto_bond: false,
            bond_tolerance: 1.001,
            bond_radius_ratio: 1.0,
            ghost_cutoff_multiplier: 2.5,
            youngs_modulus: None,
            shear_modulus: None,
            normal_stiffness: 0.0,
            shear_stiffness: 0.0,
            twist_stiffness: 0.0,
            bending_stiffness: 0.0,
            beta_normal: 0.0,
            beta_shear: 0.0,
            beta_twist: 0.0,
            beta_bending: 0.0,
            normal_damping: None,
            shear_damping: None,
            twist_damping: None,
            bending_damping: None,
            breakage: None,
            seed: 0,
            plasticity: None,
            file: None,
            format: None,
        }
    }
}

// ── BondHistoryStore ────────────────────────────────────────────────────────

/// Per-bond history: accumulated shear displacement, total rotation angle,
/// per-bond sampled failure thresholds, and plastic-state vectors.
///
/// `delta_t` is re-projected perpendicular to the current bond axis each step.
/// `delta_theta` is a running integral of `ω_rel · dt`; its split into twist
/// (along n̂) and bending (⊥ n̂) components is done on-the-fly each step.
/// `thresholds` are drawn once at bond creation from
/// [`breakage::ThresholdDistribution`] and never resampled; their meaning
/// depends on the active criterion (see [`breakage::BondThresholds`]).
/// `theta_p_bend` and `eps_p_axial` are the kinematic plastic anchors for
/// the bending and axial channels respectively (see [`plasticity`]); zero
/// when the corresponding channel is purely elastic. `theta_max_bend` and
/// `eps_max_axial` are the largest kinematic strain magnitudes the bond
/// has ever reached on each channel, used to evaluate the monotonic
/// hardening envelope under kinematic hardening.
#[derive(Clone, Debug)]
pub struct BondHistoryEntry {
    /// Global tag of the bonded partner.
    pub partner_tag: u32,
    /// Accumulated tangential displacement **Δs** (m), ⊥ to current bond axis.
    pub delta_t: [f64; 3],
    /// Accumulated relative rotation angle **Δθ** (rad) — full vector; the
    /// along-n̂ and ⊥-n̂ components are extracted each step.
    pub delta_theta: [f64; 3],
    /// Per-bond sampled failure thresholds, indexed by criterion (up to four).
    pub thresholds: [f64; 4],
    /// Plastic-bending kinematic anchor, ⊥ to bond axis. Zero ⇒ no plastic flow yet.
    pub theta_p_bend: [f64; 3],
    /// Plastic-axial kinematic anchor (signed scalar). Zero ⇒ no plastic flow yet.
    pub eps_p_axial: f64,
    /// Largest `|θ_bend|` ever reached on this bond (≥ 0).
    pub theta_max_bend: f64,
    /// Largest `|ε_axial|` ever reached on this bond (≥ 0).
    pub eps_max_axial: f64,
}

/// Per-atom list of [`BondHistoryEntry`], kept in sync with [`BondStore`].
pub struct BondHistoryStore {
    /// Outer index = local atom index; inner vec = one entry per bond on that atom.
    pub history: Vec<Vec<BondHistoryEntry>>,
}

impl BondHistoryStore {
    /// Creates an empty bond history store.
    pub fn new() -> Self {
        BondHistoryStore { history: Vec::new() }
    }
}

impl AtomData for BondHistoryStore {
    fn as_any(&self) -> &dyn Any { self }
    fn as_any_mut(&mut self) -> &mut dyn Any { self }

    fn truncate(&mut self, n: usize) {
        self.history.resize_with(n, Vec::new);
        self.history.truncate(n);
    }

    fn swap_remove(&mut self, i: usize) {
        if i < self.history.len() {
            self.history.swap_remove(i);
        }
    }

    fn apply_permutation(&mut self, perm: &[usize], n: usize) {
        let new_history: Vec<Vec<BondHistoryEntry>> =
            perm.iter().map(|&p| self.history[p].clone()).collect();
        self.history[..n].clone_from_slice(&new_history);
    }

    /// Wire format: `[count, (partner_tag, dt[3], dθ[3], thr[4], θ_p_bend[3],
    /// ε_p_axial, θ_max_bend, ε_max_axial) × count]` — `1 + 17 × count` f64s.
    fn pack(&self, i: usize, buf: &mut Vec<f64>) {
        if i < self.history.len() {
            let list = &self.history[i];
            buf.push(list.len() as f64);
            for e in list {
                buf.push(e.partner_tag as f64);
                buf.push(e.delta_t[0]); buf.push(e.delta_t[1]); buf.push(e.delta_t[2]);
                buf.push(e.delta_theta[0]); buf.push(e.delta_theta[1]); buf.push(e.delta_theta[2]);
                buf.push(e.thresholds[0]); buf.push(e.thresholds[1]);
                buf.push(e.thresholds[2]); buf.push(e.thresholds[3]);
                buf.push(e.theta_p_bend[0]); buf.push(e.theta_p_bend[1]); buf.push(e.theta_p_bend[2]);
                buf.push(e.eps_p_axial);
                buf.push(e.theta_max_bend);
                buf.push(e.eps_max_axial);
            }
        } else {
            buf.push(0.0);
        }
    }

    fn unpack(&mut self, buf: &[f64]) -> usize {
        let count = buf[0] as usize;
        let mut list = Vec::with_capacity(count);
        let mut pos = 1;
        for _ in 0..count {
            let partner_tag = buf[pos] as u32;
            let delta_t = [buf[pos + 1], buf[pos + 2], buf[pos + 3]];
            let delta_theta = [buf[pos + 4], buf[pos + 5], buf[pos + 6]];
            let thresholds = [buf[pos + 7], buf[pos + 8], buf[pos + 9], buf[pos + 10]];
            let theta_p_bend = [buf[pos + 11], buf[pos + 12], buf[pos + 13]];
            let eps_p_axial = buf[pos + 14];
            let theta_max_bend = buf[pos + 15];
            let eps_max_axial = buf[pos + 16];
            list.push(BondHistoryEntry {
                partner_tag, delta_t, delta_theta, thresholds,
                theta_p_bend, eps_p_axial,
                theta_max_bend, eps_max_axial,
            });
            pos += 17;
        }
        self.history.push(list);
        pos
    }
}

// ── BondBreakage (active criterion) ─────────────────────────────────────────

/// Active breakage criterion plus the global seed used to derive per-bond
/// threshold draws. Built once at setup from [`BondConfig::breakage`] and
/// `BondConfig::seed`, consumed by [`init_bond_history`] and [`bond_force`].
///
/// Per-bond `u`-sample draws are **deterministic in `(min(tag_a, tag_b),
/// max(tag_a, tag_b), seed)`** via [`breakage::per_bond_uniform_samples`]
/// — independent of which rank owns the bond, in which order bonds are
/// visited, and how the simulation is decomposed across MPI processes.
pub struct BondBreakage {
    /// Run-time criterion trait object. `Unbreakable` when no `[bonds.breakage]`
    /// section is present.
    pub criterion: Box<dyn BreakageCriterion>,
    /// Global seed for the per-bond threshold sampler. Two simulations with
    /// the same seed and the same bond topology produce identical per-bond
    /// thresholds even when run with different MPI decompositions.
    pub seed: u64,
}

impl Default for BondBreakage {
    fn default() -> Self {
        BondBreakage {
            criterion: Box::new(Unbreakable),
            seed: 0,
        }
    }
}

/// Setup system: builds [`BondBreakage`] from [`BondConfig::breakage`] and
/// `BondConfig::seed`. Runs once during `PostSetup`, after [`Config::load`]
/// has populated `BondConfig`.
pub fn init_breakage(bond_config: Res<BondConfig>, mut breakage: ResMut<BondBreakage>) {
    breakage.criterion = match &bond_config.breakage {
        Some(cfg) => cfg.build(),
        None => Box::new(Unbreakable),
    };
    breakage.seed = bond_config.seed;
}

// ── BondPlasticity (active plasticity model) ────────────────────────────────

/// Active plasticity model. Built once at setup from [`BondConfig::plasticity`]
/// and consumed by [`bond_force`].
#[derive(Default)]
pub struct BondPlasticity {
    /// Run-time plasticity state. Defaults to no plasticity.
    pub model: BondPlasticityModel,
}

/// Setup system: builds [`BondPlasticity`] from [`BondConfig::plasticity`].
/// Runs once during `PostSetup`.
pub fn init_plasticity(bond_config: Res<BondConfig>, mut plasticity: ResMut<BondPlasticity>) {
    plasticity.model = BondPlasticityModel::from_config(
        bond_config.plasticity.as_ref(),
        bond_config.youngs_modulus,
    );
}

// ── BondMetrics ─────────────────────────────────────────────────────────────

/// Accumulated per-step bond metrics, exposed via the thermo system.
#[derive(Default)]
pub struct BondMetrics {
    /// Sum of `δ / r₀` (axial strain) over all active bonds this step.
    pub strain_sum: f64,
    /// Number of active bonds evaluated this step.
    pub bond_count: usize,
    /// Number of bonds broken during this step.
    pub bonds_broken_this_step: usize,
    /// Cumulative number of bonds broken since the start of the simulation.
    pub total_bonds_broken: usize,
    /// Number of bonds skipped this step because the partner atom was not
    /// present as a local or ghost. Non-zero means `ghost_cutoff` is too
    /// small — bump `ghost_cutoff_multiplier` in `[bonds]`.
    pub missing_partner_skips: usize,
    /// Whether a rank-0 warning has already been printed; prevents flooding.
    pub warned_missing_partner: bool,
}

// ── Plugin ──────────────────────────────────────────────────────────────────

/// Plugin that enables BPM bond forces between DEM particles.
pub struct DemBondPlugin;

impl Plugin for DemBondPlugin {
    fn default_config(&self) -> Option<&str> {
        Some(
            r#"[bonds]
# auto_bond = false
# bond_tolerance = 1.001
# bond_radius_ratio = 1.0
# ghost_cutoff_multiplier = 2.5   # MPI: extends ghost skin to cover bond + 1-3 reach
#
# Material mode (paper-standard beam theory):
# youngs_modulus = 1.0e9
# shear_modulus  = 4.0e8
#
# Direct stiffness overrides (used when E/G are not set):
# normal_stiffness  = 0.0   # N/m
# shear_stiffness   = 0.0   # N/m
# twist_stiffness   = 0.0   # N·m/rad
# bending_stiffness = 0.0   # N·m/rad
#
# Damping ratios (critical = 1.0):
# beta_normal  = 0.0
# beta_shear   = 0.0
# beta_twist   = 0.0
# beta_bending = 0.0
#
# Raw damping overrides (optional):
# normal_damping  = 0.0
# shear_damping   = 0.0
# twist_damping   = 0.0
# bending_damping = 0.0
#
# Per-bond Weibull RNG seed (deterministic given the seed):
# seed = 0
#
# Breakage criterion (omit to leave bonds unbreakable). See
# `breakage::BreakageConfig` for the variant menu. Example — stress-based
# combined (Guo / Potyondy-Cundall):
# [bonds.breakage]
# kind    = "combined_stress"
# tensile = { kind = "constant", value = 5.0e7 }
# shear   = { kind = "constant", value = 3.0e7 }
#
# Plasticity (omit for purely elastic bonds). Both `bending` and `axial`
# channels are independently optional. Examples:
# [bonds.plasticity.bending]
# kind         = "guo_bending"
# yield_stress = 1.23e8
# [bonds.plasticity.axial]
# kind                = "piecewise"
# breakpoint_strains  = [0.01, 0.02, 0.03]
# slope_multipliers   = [0.5, 0.1, 0.0]"#,
        )
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(soil_core::BondPlugin);
        app.add_plugins(VirialStressPlugin);
        Config::load::<BondConfig>(app, "bonds");
        app.add_resource(BondMetrics::default());
        app.add_resource(BondBreakage::default());
        app.add_resource(BondPlasticity::default());
        soil_core::register_atom_data!(app, BondHistoryStore::new());
        app.add_setup_system(
            auto_bond_touching
                .label("auto_bond_touching")
                .run_if(first_stage_only()),
            ScheduleSetupSet::PostSetup,
        );
        app.add_setup_system(
            load_bonds_from_file
                .label("load_bonds_from_file")
                .run_if(first_stage_only()),
            ScheduleSetupSet::PostSetup,
        );
        // ghost_cutoff must cover bond length (bonded partners across ranks)
        // AND 2× bond length (shared-neighbour 1-3 exclusion at rank boundaries).
        // Must run AFTER bonds are created and BEFORE neighbor_setup bakes the
        // ghost_cutoff into the bin grid / borders skin.
        app.add_setup_system(
            extend_ghost_cutoff_for_bonds
                .label("extend_ghost_cutoff_for_bonds")
                .after("auto_bond_touching")
                .after("load_bonds_from_file")
                .before("neighbor_setup")
                .run_if(first_stage_only()),
            ScheduleSetupSet::PostSetup,
        );
        app.add_setup_system(
            init_breakage
                .label("init_breakage")
                .run_if(first_stage_only()),
            ScheduleSetupSet::PostSetup,
        );
        app.add_setup_system(
            init_plasticity
                .label("init_plasticity")
                .run_if(first_stage_only()),
            ScheduleSetupSet::PostSetup,
        );
        app.add_setup_system(
            init_bond_history
                .after("init_breakage")
                .after("auto_bond_touching")
                .after("load_bonds_from_file")
                .run_if(first_stage_only()),
            ScheduleSetupSet::PostSetup,
        );
        app.add_update_system(zero_bond_metrics, ParticleSimScheduleSet::PreForce);
        app.add_update_system(bond_force.label("dem_bond_force"), ParticleSimScheduleSet::Force);
        app.add_update_system(output_bond_metrics, ParticleSimScheduleSet::PostForce);
    }
}

// ── Setup systems ───────────────────────────────────────────────────────────

/// Auto-bond initially touching particles at setup time. Uses periodic
/// minimum-image distances on axes flagged periodic in `Domain`, so a
/// single-rank periodic chain whose endpoints are within bond reach across
/// the wrap (e.g. a closed-loop fibre) gets its seam bond like any other.
pub fn auto_bond_touching(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    bond_config: Res<BondConfig>,
    comm: Res<CommResource>,
    domain: Res<soil_core::Domain>,
) {
    if !bond_config.auto_bond { return; }

    let dem = registry.expect::<DemAtom>("auto_bond_touching");
    let mut bond_store = registry.expect_mut::<BondStore>("auto_bond_touching");

    let nlocal = atoms.nlocal as usize;
    while bond_store.bonds.len() < nlocal {
        bond_store.bonds.push(Vec::new());
    }

    let tol = bond_config.bond_tolerance;
    let pflags = domain.periodic_flags();
    let box_size = domain.size;
    let mut bond_count = 0u64;

    for i in 0..nlocal {
        for j in (i + 1)..nlocal {
            let mut d = [
                atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64,
                atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64,
                atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64,
            ];
            for k in 0..3 {
                if pflags[k] {
                    d[k] -= box_size[k] * (d[k] / box_size[k]).round();
                }
            }
            let dist = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]).sqrt();
            let sum_r = dem.radius[i] + dem.radius[j];
            if dist <= sum_r * tol {
                bond_store.bonds[i].push(BondEntry {
                    partner_tag: atoms.tag[j], bond_type: 0, r0: dist,
                });
                bond_store.bonds[j].push(BondEntry {
                    partner_tag: atoms.tag[i], bond_type: 0, r0: dist,
                });
                bond_count += 1;
            }
        }
    }

    if comm.rank() == 0 {
        println!("DemBond: auto-bonded {} pairs", bond_count);
    }
}

/// Load bonds from a LAMMPS data file's `Bonds` section. r₀ is computed
/// using periodic minimum-image distance on axes flagged periodic in
/// `Domain`, so explicit bonds that wrap across a periodic seam get the
/// correct rest length (= shortest periodic distance, not the full
/// across-the-box distance).
pub fn load_bonds_from_file(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    bond_config: Res<BondConfig>,
    comm: Res<CommResource>,
    domain: Res<soil_core::Domain>,
) {
    let file_path = match bond_config.file.as_deref() {
        Some(p) => p,
        None => return,
    };
    let format = bond_config.format.as_deref().unwrap_or("lammps_data");
    if format != "lammps_data" {
        eprintln!("ERROR: Unsupported bond file format '{}'. Supported: lammps_data", format);
        std::process::exit(1);
    }

    let file = File::open(file_path).unwrap_or_else(|e| {
        eprintln!("ERROR: Failed to open bond file '{}': {}", file_path, e);
        std::process::exit(1);
    });
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines()
        .map(|l| l.expect("failed to read line from bond file"))
        .collect();

    let section_headers = [
        "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
        "Masses", "Pair Coeffs",
    ];
    let is_section_header = |line: &str| -> bool {
        let t = line.trim();
        section_headers.iter().any(|h| t.starts_with(h))
    };

    let mut bonds_start = None;
    for (i, line) in lines.iter().enumerate() {
        if line.trim().starts_with("Bonds") {
            bonds_start = Some(i + 1);
        }
    }
    let bonds_start = match bonds_start {
        Some(s) => s,
        None => {
            if comm.rank() == 0 {
                println!("DemBond: no Bonds section found in '{}', skipping", file_path);
            }
            return;
        }
    };

    let nlocal = atoms.nlocal as usize;
    let mut tag_to_local: HashMap<u32, usize> = HashMap::with_capacity(nlocal);
    for i in 0..nlocal {
        tag_to_local.insert(atoms.tag[i], i);
    }

    let mut bond_store = registry.expect_mut::<BondStore>("load_bonds_from_file");
    while bond_store.bonds.len() < nlocal {
        bond_store.bonds.push(Vec::new());
    }

    let mut bond_count = 0u64;

    for i in bonds_start..lines.len() {
        let t = lines[i].trim();
        if t.is_empty() { continue; }
        if is_section_header(t) { break; }
        if t.starts_with('#') { continue; }

        let fields: Vec<&str> = t.split_whitespace().collect();
        if fields.len() < 4 { continue; }

        let bond_type: u32 = fields[1].parse().unwrap_or(0);
        let tag1: u32 = fields[2].parse().expect("failed to parse atom tag1 in Bonds section");
        let tag2: u32 = fields[3].parse().expect("failed to parse atom tag2 in Bonds section");

        let idx1 = match tag_to_local.get(&tag1) { Some(&i) => i, None => continue };
        let idx2 = match tag_to_local.get(&tag2) { Some(&i) => i, None => continue };

        let mut d = [
            atoms.pos[idx2][0] as f64 - atoms.pos[idx1][0] as f64,
            atoms.pos[idx2][1] as f64 - atoms.pos[idx1][1] as f64,
            atoms.pos[idx2][2] as f64 - atoms.pos[idx1][2] as f64,
        ];
        let pflags = domain.periodic_flags();
        let box_size = domain.size;
        for k in 0..3 {
            if pflags[k] {
                d[k] -= box_size[k] * (d[k] / box_size[k]).round();
            }
        }
        let r0 = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]).sqrt();

        bond_store.bonds[idx1].push(BondEntry { partner_tag: tag2, bond_type, r0 });
        bond_store.bonds[idx2].push(BondEntry { partner_tag: tag1, bond_type, r0 });
        bond_count += 1;
    }

    if comm.rank() == 0 {
        println!(
            "DemBond: loaded {} bonds from LAMMPS data file '{}'",
            bond_count, file_path
        );
    }
}

/// Extends `Domain::ghost_cutoff` so bonded partners (and their 1-3
/// shared neighbours) remain visible as ghosts when atoms migrate across
/// MPI rank boundaries.
///
/// Without this, a bond spanning a rank boundary silently "disappears"
/// from the force loop: `bond_force` can't resolve the partner tag into a
/// local/ghost index, so the bond term is skipped and the atom drifts free.
///
/// Runs in `ScheduleSetupSet::PostSetup`, ordered **after**
/// `auto_bond_touching` / `load_bonds_from_file` (so bonds exist) and
/// **before** `neighbor_setup` (which locks `ghost_cutoff` into the bin
/// grid and border skin).
pub fn extend_ghost_cutoff_for_bonds(
    registry: Res<AtomDataRegistry>,
    bond_config: Res<BondConfig>,
    mut domain: ResMut<Domain>,
    comm: Res<CommResource>,
) {
    if bond_config.ghost_cutoff_multiplier <= 0.0 {
        return;
    }

    // Global max r0 across all ranks (bonds currently only exist on the
    // rank(s) that auto-bonded or loaded them at setup).
    let local_max_r0 = {
        let bond_store = match registry.get::<BondStore>() {
            Some(bs) => bs,
            None => return,
        };
        let mut m = 0.0f64;
        for list in &bond_store.bonds {
            for b in list {
                if b.r0 > m { m = b.r0; }
            }
        }
        m
    };
    // all_reduce_max via negated min (only min is in the CommBackend trait).
    let global_max_r0 = -comm.all_reduce_min_f64(-local_max_r0);
    if global_max_r0 <= 0.0 {
        return;
    }

    let required = global_max_r0 * bond_config.ghost_cutoff_multiplier;
    if required > domain.ghost_cutoff {
        let old = domain.ghost_cutoff;
        domain.ghost_cutoff = required;
        if comm.rank() == 0 {
            println!(
                "DemBond: extended ghost_cutoff {:.6} → {:.6} (max r₀ = {:.6} × multiplier {:.2})",
                old, required, global_max_r0, bond_config.ghost_cutoff_multiplier
            );
        }
    }
}

/// Seed [`BondHistoryStore`] entries for every bond that does not already
/// have one. Each new entry's failure thresholds are drawn from
/// [`BondBreakage::criterion`] via [`breakage::per_bond_uniform_samples`],
/// which uses the bond's tag pair plus `BondBreakage::seed` — **not** a
/// rank-local RNG. So MPI-decomposition order does not affect which
/// thresholds end up on which bond; the same physical bond gets the same
/// thresholds on every rank that ever sees it.
///
/// The atom's local index is read from [`Atom::tag`] so a bond's two tags
/// can be looked up without needing the global topology to be present on
/// the local rank.
pub fn init_bond_history(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    breakage: Res<BondBreakage>,
) {
    let bond_store = registry.get::<BondStore>();
    let bonds = match bond_store {
        Some(ref b) => b,
        None => return,
    };

    let mut history = registry.expect_mut::<BondHistoryStore>("init_bond_history");
    while history.history.len() < bonds.bonds.len() {
        history.history.push(Vec::new());
    }
    let nlocal = atoms.nlocal as usize;
    for i in 0..bonds.bonds.len().min(nlocal) {
        let tag_a = atoms.tag[i];
        for bond in &bonds.bonds[i] {
            let has = history.history[i].iter().any(|h| h.partner_tag == bond.partner_tag);
            if !has {
                // MPI-stable: same (tag pair, seed) → same `u` on every rank.
                let u = breakage::per_bond_uniform_samples(
                    tag_a, bond.partner_tag, breakage.seed,
                );
                let thr = breakage.criterion.sample(bond.r0, u);
                history.history[i].push(BondHistoryEntry {
                    partner_tag: bond.partner_tag,
                    delta_t: [0.0; 3],
                    delta_theta: [0.0; 3],
                    thresholds: thr.t,
                    theta_p_bend: [0.0; 3],
                    eps_p_axial: 0.0,
                    theta_max_bend: 0.0,
                    eps_max_axial: 0.0,
                });
            }
        }
    }
}

// ── Force system ────────────────────────────────────────────────────────────

/// Computes BPM bond forces and moments for all local atoms.
pub fn bond_force(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    bond_config: Res<BondConfig>,
    breakage: Res<BondBreakage>,
    plasticity: Res<BondPlasticity>,
    mut metrics: ResMut<BondMetrics>,
    mut virial: Option<ResMut<VirialStress>>,
    domain: Res<soil_core::Domain>,
) {
    let bond_store = registry.get::<BondStore>();
    let bonds = match bond_store {
        Some(ref b) => b,
        None => return,
    };

    let pflags = domain.periodic_flags();
    let box_size = domain.size;

    let nlocal = atoms.nlocal as usize;
    if bonds.bonds.len() < nlocal { return; }

    let ratio = bond_config.bond_radius_ratio;
    if ratio <= 0.0 { return; }

    // Material-mode (E, G) or direct stiffness fallback.
    let e_mod = bond_config.youngs_modulus;
    let g_mod = bond_config.shear_modulus;
    let k_n_direct = bond_config.normal_stiffness;
    let k_t_direct = bond_config.shear_stiffness;
    let k_tor_direct = bond_config.twist_stiffness;
    let k_bend_direct = bond_config.bending_stiffness;

    let beta_n = bond_config.beta_normal;
    let beta_t = bond_config.beta_shear;
    let beta_tor = bond_config.beta_twist;
    let beta_bend = bond_config.beta_bending;

    let dt = atoms.dt;

    // Always need DemAtom (radius, omega, torque, inv_inertia); always need history
    // for shear and rotation channels.
    let mut dem = registry.expect_mut::<DemAtom>("bond_force");
    let mut hist = registry.expect_mut::<BondHistoryStore>("bond_force");

    // Tag → index lookup (local + ghost)
    let mut tag_to_index: HashMap<u32, usize> = HashMap::with_capacity(atoms.len());
    for idx in 0..atoms.len() {
        tag_to_index.insert(atoms.tag[idx], idx);
    }

    let mut bonds_to_break: Vec<(u32, u32)> = Vec::new();

    for i in 0..nlocal {
        for b_idx in 0..bonds.bonds[i].len() {
            let bond = &bonds.bonds[i][b_idx];
            let j = match tag_to_index.get(&bond.partner_tag) {
                Some(&idx) => idx,
                None => {
                    metrics.missing_partner_skips += 1;
                    continue;
                }
            };
            // Process each bond once: lower tag owns the computation.
            if atoms.tag[i] >= bond.partner_tag { continue; }

            // Minimum-image distance on periodic axes. The tag_to_index
            // lookup above can resolve a partner to a ghost copy at a
            // wrap-image position when the same tag appears as both local
            // and ghost; using minimum-image normalises that so within-chain
            // bonds always see the short distance and wrap bonds always
            // see the wrapped distance regardless of which atom copy the
            // map happens to land on.
            let mut dxv = [
                atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64,
                atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64,
                atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64,
            ];
            for k in 0..3 {
                if pflags[k] {
                    dxv[k] -= box_size[k] * (dxv[k] / box_size[k]).round();
                }
            }
            let dx = dxv[0]; let dy = dxv[1]; let dz = dxv[2];
            let dist = (dx*dx + dy*dy + dz*dz).sqrt();
            if dist < 1e-20 { continue; }
            let nhat = [dx / dist, dy / dist, dz / dist];

            // Bond geometry (cylindrical beam).
            let r_b = ratio * dem.radius[i].min(dem.radius[j]);
            let area = PI * r_b * r_b;
            let jpol = 0.5 * PI * r_b.powi(4);        // polar 2nd moment of area
            let iben = 0.5 * jpol;                    // bending 2nd moment (½ J)
            let len = bond.r0;

            // Stiffnesses — material mode wins when E/G provided.
            let k_n = match e_mod { Some(e) => e * area / len, None => k_n_direct };
            let k_t = match g_mod { Some(g) => g * area / len, None => k_t_direct };
            let k_tor = match g_mod { Some(g) => g * jpol / len, None => k_tor_direct };
            let k_bend = match e_mod { Some(e) => e * iben / len, None => k_bend_direct };

            // Reduced mass / reduced MOI for damping.
            let m_i = atoms.mass[i] as f64;
            let m_j = atoms.mass[j] as f64;
            let m_red = if m_i + m_j > 0.0 { m_i * m_j / (m_i + m_j) } else { 0.0 };
            let moi_i = if dem.inv_inertia[i] > 0.0 { 1.0 / dem.inv_inertia[i] } else { 0.0 };
            let moi_j = if dem.inv_inertia[j] > 0.0 { 1.0 / dem.inv_inertia[j] } else { 0.0 };
            let moi_red = if moi_i + moi_j > 0.0 { moi_i * moi_j / (moi_i + moi_j) } else { 0.0 };

            // Damping: raw override if provided, else critical-damping formula.
            let gamma_n = bond_config.normal_damping
                .unwrap_or_else(|| 2.0 * beta_n * (m_red * k_n.max(0.0)).sqrt());
            let gamma_t = bond_config.shear_damping
                .unwrap_or_else(|| 2.0 * beta_t * (m_red * k_t.max(0.0)).sqrt());
            let gamma_tor = bond_config.twist_damping
                .unwrap_or_else(|| 2.0 * beta_tor * (moi_red * k_tor.max(0.0)).sqrt());
            let gamma_bend = bond_config.bending_damping
                .unwrap_or_else(|| 2.0 * beta_bend * (moi_red * k_bend.max(0.0)).sqrt());

            // Kinematics at contact mid-point (lever arm = L/2 n̂).
            let half_l = 0.5 * len;
            let r1 = [half_l * nhat[0], half_l * nhat[1], half_l * nhat[2]]; // from i → contact
            // ω × r
            let w_i = dem.omega[i];
            let w_j = dem.omega[j];
            let v_i_c = [
                atoms.vel[i][0] as f64 + w_i[1]*r1[2] - w_i[2]*r1[1],
                atoms.vel[i][1] as f64 + w_i[2]*r1[0] - w_i[0]*r1[2],
                atoms.vel[i][2] as f64 + w_i[0]*r1[1] - w_i[1]*r1[0],
            ];
            // r2 = -r1 for j → contact
            let v_j_c = [
                atoms.vel[j][0] as f64 - (w_j[1]*r1[2] - w_j[2]*r1[1]),
                atoms.vel[j][1] as f64 - (w_j[2]*r1[0] - w_j[0]*r1[2]),
                atoms.vel[j][2] as f64 - (w_j[0]*r1[1] - w_j[1]*r1[0]),
            ];
            let v_rel = [v_j_c[0] - v_i_c[0], v_j_c[1] - v_i_c[1], v_j_c[2] - v_i_c[2]];
            let v_n_s = v_rel[0]*nhat[0] + v_rel[1]*nhat[1] + v_rel[2]*nhat[2];
            let v_n = [v_n_s*nhat[0], v_n_s*nhat[1], v_n_s*nhat[2]];
            let v_t = [v_rel[0] - v_n[0], v_rel[1] - v_n[1], v_rel[2] - v_n[2]];

            // Axial elongation (kinematic).
            let delta = dist - bond.r0;

            // ── Locate (or create) the history entry for this bond ──
            //
            // Defensive fallback: `init_bond_history` should have seeded every
            // entry (with sampled thresholds) at setup. If an entry is missing
            // here we create one with zero thresholds — which trips any
            // non-`Unbreakable` criterion immediately. That fail-loud behaviour
            // surfaces missed bond-creation paths instead of silently making
            // bonds unbreakable.
            while hist.history.len() <= i { hist.history.push(Vec::new()); }
            let h_idx = match hist.history[i].iter().position(|h| h.partner_tag == bond.partner_tag) {
                Some(idx) => idx,
                None => {
                    hist.history[i].push(BondHistoryEntry {
                        partner_tag: bond.partner_tag,
                        delta_t: [0.0; 3],
                        delta_theta: [0.0; 3],
                        thresholds: [0.0; 4],
                        theta_p_bend: [0.0; 3],
                        eps_p_axial: 0.0,
                        theta_max_bend: 0.0,
                        eps_max_axial: 0.0,
                    });
                    hist.history[i].len() - 1
                }
            };

            // Normal force conservative part: elastic `K_n · δ` by default,
            // or piecewise-linear plasticity return-map on the axial channel
            // when configured. Damping `γ_n · v_n_s` is added on top.
            let (f_n_cons_mag, eps_p_axial_new, eps_max_axial_new) =
                if let Some(axial_model) = plasticity.model.axial.as_ref() {
                    let eps_axial = if bond.r0 > 0.0 { delta / bond.r0 } else { 0.0 };
                    plasticity::update_axial(
                        eps_axial,
                        hist.history[i][h_idx].eps_p_axial,
                        hist.history[i][h_idx].eps_max_axial,
                        k_n,
                        bond.r0,
                        axial_model,
                    )
                } else {
                    (
                        k_n * delta,
                        hist.history[i][h_idx].eps_p_axial,
                        hist.history[i][h_idx].eps_max_axial,
                    )
                };
            hist.history[i][h_idx].eps_p_axial = eps_p_axial_new;
            hist.history[i][h_idx].eps_max_axial = eps_max_axial_new;
            let f_n_mag = f_n_cons_mag + gamma_n * v_n_s;
            let f_n = [f_n_mag * nhat[0], f_n_mag * nhat[1], f_n_mag * nhat[2]];

            // Shear: re-project Δs ⊥ to new n̂, then integrate.
            {
                let h = &mut hist.history[i][h_idx];
                let s_n = h.delta_t[0]*nhat[0] + h.delta_t[1]*nhat[1] + h.delta_t[2]*nhat[2];
                h.delta_t[0] -= s_n * nhat[0];
                h.delta_t[1] -= s_n * nhat[1];
                h.delta_t[2] -= s_n * nhat[2];
                h.delta_t[0] += v_t[0] * dt;
                h.delta_t[1] += v_t[1] * dt;
                h.delta_t[2] += v_t[2] * dt;
            }
            let ds = hist.history[i][h_idx].delta_t;
            // Bond-internal shear force (same sign convention as Fortran:
            // grows with +Δs and +v_t). Applied as +f_t on atom i (lower tag)
            // and −f_t on atom j: when atom j slides below atom i, +Δs is
            // negative, so f_t points downward on atom i (pulls it toward
            // atom j) and upward on atom j (pulls it back into alignment).
            let f_t = [
                k_t * ds[0] + gamma_t * v_t[0],
                k_t * ds[1] + gamma_t * v_t[1],
                k_t * ds[2] + gamma_t * v_t[2],
            ];

            // Rotation kinematics
            let w_rel = [w_j[0] - w_i[0], w_j[1] - w_i[1], w_j[2] - w_i[2]];
            let w_rel_n_s = w_rel[0]*nhat[0] + w_rel[1]*nhat[1] + w_rel[2]*nhat[2];
            let w_n = [w_rel_n_s*nhat[0], w_rel_n_s*nhat[1], w_rel_n_s*nhat[2]];
            let w_t = [w_rel[0] - w_n[0], w_rel[1] - w_n[1], w_rel[2] - w_n[2]];

            // Update Δθ and split into twist (along n̂) and bending (⊥ n̂) parts.
            {
                let h = &mut hist.history[i][h_idx];
                h.delta_theta[0] += w_rel[0] * dt;
                h.delta_theta[1] += w_rel[1] * dt;
                h.delta_theta[2] += w_rel[2] * dt;
            }
            let dth = hist.history[i][h_idx].delta_theta;
            let dth_n_s = dth[0]*nhat[0] + dth[1]*nhat[1] + dth[2]*nhat[2];
            let dth_twist = [dth_n_s*nhat[0], dth_n_s*nhat[1], dth_n_s*nhat[2]];
            let dth_bend  = [dth[0] - dth_twist[0], dth[1] - dth_twist[1], dth[2] - dth_twist[2]];

            // Bond-internal twist and bending moments (positive sign: the
            // magnitudes grow with ω_rel and Δθ_rel). Applied as +m on atom i
            // (lower tag) and −m on atom j (higher tag), matching the Fortran
            // reference. This damps relative rotation: atom j receives a
            // torque −γ·ω_rel that opposes its rotation, while atom i receives
            // +γ·ω_rel that accelerates it toward the same rotation — in both
            // cases reducing the *relative* angular velocity.
            let m_tor = [
                k_tor * dth_twist[0] + gamma_tor * w_n[0],
                k_tor * dth_twist[1] + gamma_tor * w_n[1],
                k_tor * dth_twist[2] + gamma_tor * w_n[2],
            ];
            // Bending: conservative part is elastic by default, or follows
            // the configured plasticity return-map (Guo / piecewise) when
            // bending plasticity is configured. Damping is added on top.
            let (m_bend_cons, theta_p_bend_new, theta_max_bend_new) =
                if let Some(bending_model) = plasticity.model.bending.as_ref() {
                    plasticity::update_bending(
                        dth_bend,
                        hist.history[i][h_idx].theta_p_bend,
                        hist.history[i][h_idx].theta_max_bend,
                        k_bend,
                        bending_model,
                        r_b,
                        bond.r0,
                    )
                } else {
                    let m_elastic = [
                        k_bend * dth_bend[0],
                        k_bend * dth_bend[1],
                        k_bend * dth_bend[2],
                    ];
                    (
                        m_elastic,
                        hist.history[i][h_idx].theta_p_bend,
                        hist.history[i][h_idx].theta_max_bend,
                    )
                };
            hist.history[i][h_idx].theta_p_bend = theta_p_bend_new;
            hist.history[i][h_idx].theta_max_bend = theta_max_bend_new;
            let m_bend = [
                m_bend_cons[0] + gamma_bend * w_t[0],
                m_bend_cons[1] + gamma_bend * w_t[1],
                m_bend_cons[2] + gamma_bend * w_t[2],
            ];

            // Breakage: build the geom/loads/kinematics snapshots and ask the
            // active criterion whether this bond has failed.
            let m_bend_mag = (m_bend[0]*m_bend[0] + m_bend[1]*m_bend[1] + m_bend[2]*m_bend[2]).sqrt();
            let m_tor_mag  = (m_tor[0]*m_tor[0]  + m_tor[1]*m_tor[1]  + m_tor[2]*m_tor[2]).sqrt();
            let f_t_mag    = (f_t[0]*f_t[0] + f_t[1]*f_t[1] + f_t[2]*f_t[2]).sqrt();
            let ds_mag     = (ds[0]*ds[0] + ds[1]*ds[1] + ds[2]*ds[2]).sqrt();
            let dth_bend_mag =
                (dth_bend[0]*dth_bend[0] + dth_bend[1]*dth_bend[1] + dth_bend[2]*dth_bend[2]).sqrt();
            let l_now = dist.max(f64::MIN_POSITIVE);
            let geom = breakage::BondGeom { r_b, area, iben, jpol, l0: bond.r0 };
            let loads = breakage::BondLoads { f_n: f_n_mag, f_t_mag, m_bend_mag, m_tor_mag };
            let kin = breakage::BondKinematics {
                eps_axial:   delta / bond.r0,
                gamma_shear: ds_mag / l_now,
                kappa_bend:  dth_bend_mag / l_now,
                kappa_tor:   dth_n_s.abs() / l_now,
            };
            let thr = breakage::BondThresholds { t: hist.history[i][h_idx].thresholds };
            if breakage.criterion.check(&geom, &loads, &kin, &thr).is_some() {
                bonds_to_break.push((atoms.tag[i], bond.partner_tag));
                continue;
            }

            // ── Apply forces ──
            let f_total = [f_n[0] + f_t[0], f_n[1] + f_t[1], f_n[2] + f_t[2]];
            atoms.force[i][0] += f_total[0] as soil_core::Accum;
            atoms.force[i][1] += f_total[1] as soil_core::Accum;
            atoms.force[i][2] += f_total[2] as soil_core::Accum;
            atoms.force[j][0] -= f_total[0] as soil_core::Accum;
            atoms.force[j][1] -= f_total[1] as soil_core::Accum;
            atoms.force[j][2] -= f_total[2] as soil_core::Accum;

            if let Some(ref mut v) = virial {
                if v.active {
                    v.add_pair(dx, dy, dz, f_total[0], f_total[1], f_total[2]);
                }
            }

            // Torque from shear at lever arm (both particles get r1 × f_t).
            let tau_shear = [
                r1[1]*f_t[2] - r1[2]*f_t[1],
                r1[2]*f_t[0] - r1[0]*f_t[2],
                r1[0]*f_t[1] - r1[1]*f_t[0],
            ];

            // Total torque: +M on i, −M on j; shear torque same sign on both.
            let m_total = [m_tor[0] + m_bend[0], m_tor[1] + m_bend[1], m_tor[2] + m_bend[2]];
            dem.torque[i][0] += tau_shear[0] + m_total[0];
            dem.torque[i][1] += tau_shear[1] + m_total[1];
            dem.torque[i][2] += tau_shear[2] + m_total[2];
            dem.torque[j][0] += tau_shear[0] - m_total[0];
            dem.torque[j][1] += tau_shear[1] - m_total[1];
            dem.torque[j][2] += tau_shear[2] - m_total[2];

            metrics.strain_sum += delta / bond.r0;
            metrics.bond_count += 1;
        }
    }

    if !bonds_to_break.is_empty() {
        drop(bond_store);
        drop(dem);
        drop(hist);

        let mut bond_store = registry.expect_mut::<BondStore>("bond_force_break");
        let mut history_store = registry.expect_mut::<BondHistoryStore>("bond_force_break");

        for (tag_a, tag_b) in &bonds_to_break {
            for idx in 0..atoms.len() {
                if atoms.tag[idx] == *tag_a || atoms.tag[idx] == *tag_b {
                    let partner = if atoms.tag[idx] == *tag_a { *tag_b } else { *tag_a };
                    if idx < bond_store.bonds.len() {
                        bond_store.bonds[idx].retain(|b| b.partner_tag != partner);
                    }
                    if idx < history_store.history.len() {
                        history_store.history[idx].retain(|h| h.partner_tag != partner);
                    }
                }
            }
        }

        metrics.bonds_broken_this_step += bonds_to_break.len();
        metrics.total_bonds_broken += bonds_to_break.len();
    }
}

/// Reset per-step bond metrics to zero before the force computation pass.
pub fn zero_bond_metrics(mut metrics: ResMut<BondMetrics>) {
    metrics.strain_sum = 0.0;
    metrics.bond_count = 0;
    metrics.bonds_broken_this_step = 0;
    metrics.missing_partner_skips = 0;
}

/// Write bond metrics to thermo output after force computation.
///
/// Also emits a one-shot rank-0 warning if any bond this step could not find
/// its partner (ghost_cutoff too small to span the rank boundary).
pub fn output_bond_metrics(
    mut metrics: ResMut<BondMetrics>,
    comm: Res<CommResource>,
    mut thermo: Option<ResMut<Thermo>>,
) {
    let strain_sum = comm.all_reduce_sum_f64(metrics.strain_sum);
    let bond_count = comm.all_reduce_sum_f64(metrics.bond_count as f64);
    let missing_global = comm.all_reduce_sum_f64(metrics.missing_partner_skips as f64);

    if missing_global > 0.0 && !metrics.warned_missing_partner && comm.rank() == 0 {
        eprintln!(
            "WARNING: DemBond skipped {} bond(s) this step because the partner \
             was not present as a local/ghost atom on the owning rank. \
             Ghost cutoff is too small to span a bond across a rank boundary. \
             Increase [bonds].ghost_cutoff_multiplier (default 2.5) or reduce \
             MPI decomposition along the bonded direction.",
            missing_global as usize
        );
        metrics.warned_missing_partner = true;
    }

    if let Some(ref mut thermo) = thermo {
        if bond_count > 0.0 {
            thermo.set("bond_strain", strain_sum / bond_count);
        } else {
            thermo.set("bond_strain", 0.0);
        }
        thermo.set("bonds_broken", metrics.total_bonds_broken as f64);
        thermo.set("bond_missing", missing_global);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dirt_atom::DemAtom;
    use soil_core::{Atom, AtomDataRegistry, BondEntry, BondStore, CommResource, SingleProcessComm, toml};
    use dirt_test_utils::push_dem_test_atom;

    fn make_bond_config() -> BondConfig {
        BondConfig {
            normal_stiffness: 1e7,
            ..BondConfig::default()
        }
    }

    /// Build a 2-atom pair simulation app. `vel1`/`omega1` apply to atom index 1.
    fn build_pair_app_with(
        radius: f64,
        sep: f64,
        cfg: BondConfig,
        vel1: [f64; 3],
        omega1: [f64; 3],
    ) -> App {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        atom.dt = 1e-6;
        push_dem_test_atom(&mut atom, &mut dem, 1, [0.0, 0.0, 0.0], radius);
        push_dem_test_atom(&mut atom, &mut dem, 2, [sep, 0.0, 0.0], radius);
        atom.vel[1] = vel1;
        dem.omega[1] = omega1;
        atom.nlocal = 2; atom.natoms = 2;

        let mut bond_store = BondStore::new();
        bond_store.bonds.push(vec![BondEntry { partner_tag: 2, bond_type: 0, r0: 0.002 }]);
        bond_store.bonds.push(vec![BondEntry { partner_tag: 1, bond_type: 0, r0: 0.002 }]);

        // History starts empty — `init_bond_history` will seed it with
        // per-bond thresholds sampled from the active criterion at setup.
        let history = BondHistoryStore::new();

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(bond_store);
        registry.register(history);

        let mut domain = soil_core::Domain::new();
        domain.size = [10.0, 10.0, 10.0];

        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(cfg);
        app.add_resource(BondMetrics::default());
        app.add_resource(BondBreakage::default());
        app.add_resource(BondPlasticity::default());
        app.add_resource(CommResource(Box::new(SingleProcessComm::new())));
        app.add_resource(Thermo::new());
        app.add_resource(domain);
        app.add_setup_system(init_breakage.label("init_breakage"), ScheduleSetupSet::PostSetup);
        app.add_setup_system(init_plasticity.label("init_plasticity"), ScheduleSetupSet::PostSetup);
        app.add_setup_system(
            init_bond_history.after("init_breakage"),
            ScheduleSetupSet::PostSetup,
        );
        app.add_update_system(bond_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        // `app.run()` only steps the update loop; we need to run setup
        // explicitly here so `init_breakage` and `init_bond_history` seed
        // per-bond thresholds before bond_force fires.
        app.setup();
        app
    }

    fn build_pair_app(radius: f64, sep: f64, cfg: BondConfig) -> App {
        build_pair_app_with(radius, sep, cfg, [0.0; 3], [0.0; 3])
    }

    #[test]
    fn auto_bond_creates_symmetric_bonds() {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;
        push_dem_test_atom(&mut atom, &mut dem, 1, [0.0, 0.0, 0.0], radius);
        push_dem_test_atom(&mut atom, &mut dem, 2, [0.002, 0.0, 0.0], radius);
        atom.nlocal = 2; atom.natoms = 2;

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(BondStore::new());
        registry.register(BondHistoryStore::new());

        let mut domain = soil_core::Domain::new();
        domain.size = [10.0, 10.0, 10.0];

        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(BondConfig { auto_bond: true, ..make_bond_config() });
        app.add_resource(CommResource(Box::new(SingleProcessComm::new())));
        app.add_resource(SchedulerManager::default());
        app.add_resource(domain);
        app.add_setup_system(auto_bond_touching, ScheduleSetupSet::PostSetup);
        app.organize_systems();
        app.setup();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let bonds = registry.expect::<BondStore>("test");
        assert_eq!(bonds.bonds[0].len(), 1);
        assert_eq!(bonds.bonds[1].len(), 1);
        assert_eq!(bonds.bonds[0][0].partner_tag, 2);
        assert_eq!(bonds.bonds[1][0].partner_tag, 1);
    }

    #[test]
    fn auto_bond_skips_separated_atoms() {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;
        push_dem_test_atom(&mut atom, &mut dem, 1, [0.0, 0.0, 0.0], radius);
        push_dem_test_atom(&mut atom, &mut dem, 2, [0.01, 0.0, 0.0], radius);
        atom.nlocal = 2; atom.natoms = 2;

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(BondStore::new());
        registry.register(BondHistoryStore::new());

        let mut domain = soil_core::Domain::new();
        domain.size = [10.0, 10.0, 10.0];

        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(BondConfig { auto_bond: true, ..make_bond_config() });
        app.add_resource(CommResource(Box::new(SingleProcessComm::new())));
        app.add_resource(SchedulerManager::default());
        app.add_resource(domain);
        app.add_setup_system(auto_bond_touching, ScheduleSetupSet::PostSetup);
        app.organize_systems();
        app.setup();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let bonds = registry.expect::<BondStore>("test");
        assert_eq!(bonds.bonds[0].len(), 0);
        assert_eq!(bonds.bonds[1].len(), 0);
    }

    #[test]
    fn bond_force_attracts_stretched_pair() {
        let app = build_pair_app(0.001, 0.0025, make_bond_config());
        let mut app = app;
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(atom.force[0][0] > 0.0, "stretched bond attracts atom 0");
        assert!(atom.force[1][0] < 0.0, "stretched bond attracts atom 1");
        assert!((atom.force[0][0] + atom.force[1][0]).abs() < 1e-6);
    }

    #[test]
    fn bond_force_repels_compressed_pair() {
        let mut app = build_pair_app(0.001, 0.0015, make_bond_config());
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(atom.force[0][0] < 0.0, "compressed bond repels atom 0");
        assert!(atom.force[1][0] > 0.0, "compressed bond repels atom 1");
    }

    #[test]
    fn bond_force_zero_at_equilibrium() {
        let mut app = build_pair_app(0.001, 0.002, make_bond_config());
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(atom.force[0][0].abs() < 1e-10);
        assert!(atom.force[1][0].abs() < 1e-10);
    }

    #[test]
    fn tangential_bond_force_perpendicular() {
        let cfg = BondConfig {
            shear_stiffness: 5e6,
            ..make_bond_config()
        };
        let mut app = build_pair_app_with(0.001, 0.002, cfg, [0.0, 0.1, 0.0], [0.0; 3]);
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(atom.force[0][1].abs() > 0.0, "tangential force on atom 0");
        assert!(
            (atom.force[0][1] + atom.force[1][1]).abs() < 1e-6,
            "Newton's 3rd law for tangential"
        );
    }

    #[test]
    fn twist_moment_opposes_relative_twist() {
        // Relative ω along bond axis (x) should give a twist moment along x only.
        let cfg = BondConfig {
            twist_stiffness: 1e4,
            ..make_bond_config()
        };
        let mut app = build_pair_app_with(0.001, 0.002, cfg, [0.0; 3], [100.0, 0.0, 0.0]);
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // Atom 1 has +ω_x; the bond opposes its relative rotation by applying
        // a −x torque on atom 1 (slow it down) and a +x torque on atom 0
        // (speed it up in the same direction → damps *relative* rotation).
        assert!(dem.torque[1][0] < 0.0, "twist on atom 1 opposes +ω_x, got {}", dem.torque[1][0]);
        assert!(dem.torque[0][0] > 0.0, "twist on atom 0 is opposite of atom 1");
        // y/z components should be ~0 for pure twist
        assert!(dem.torque[0][1].abs() < 1e-10);
        assert!(dem.torque[0][2].abs() < 1e-10);
    }

    #[test]
    fn bending_moment_opposes_relative_bending() {
        // Relative ω perpendicular to bond axis is pure bending.
        let cfg = BondConfig {
            bending_stiffness: 1e4,
            ..make_bond_config()
        };
        let mut app = build_pair_app_with(0.001, 0.002, cfg, [0.0; 3], [0.0, 100.0, 0.0]);
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // Atom 1 has +ω_y (perpendicular to bond axis = bending).
        // Atom 1 gets a −y torque opposing its rotation; atom 0 gets +y
        // to rotate in sync, damping the relative rotation.
        assert!(dem.torque[1][1] < 0.0, "bending on atom 1 opposes +ω_y, got {}", dem.torque[1][1]);
        assert!(dem.torque[0][1] > 0.0, "bending on atom 0 is opposite of atom 1");
        // No twist moment for pure perpendicular ω
        assert!(dem.torque[0][0].abs() < 1e-10);
    }

    #[test]
    fn twist_and_bending_are_independent() {
        // Supplying only twist_stiffness with perpendicular ω should give zero moment.
        let cfg = BondConfig { twist_stiffness: 1e4, ..make_bond_config() };
        let mut app = build_pair_app_with(0.001, 0.002, cfg, [0.0; 3], [0.0, 100.0, 0.0]);
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        assert!(
            dem.torque[0][1].abs() < 1e-10,
            "perpendicular ω must produce no moment when only twist_stiffness is set"
        );
    }

    #[test]
    fn material_mode_derives_normal_stiffness_from_e_a_over_l() {
        // With E = 1 GPa, r_b = R = 0.001, L = r0 = 0.002 → K_n = E·A/L = 1e9·π·1e-6/2e-3 = π/2 · 1e6
        // Stretch by 1e-5 → force should be K_n · 1e-5.
        let e = 1e9;
        let cfg = BondConfig {
            youngs_modulus: Some(e),
            bond_radius_ratio: 1.0,
            ..BondConfig::default()
        };
        let r = 0.001;
        let l = 0.002;
        let delta = 1e-5;
        let mut app = build_pair_app(r, l + delta, cfg);
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let expected_k_n = e * PI * r * r / l;
        let expected_force = expected_k_n * delta;
        assert!(
            (atom.force[0][0] - expected_force).abs() / expected_force < 1e-6,
            "F_n got {}, expected {}",
            atom.force[0][0],
            expected_force
        );
    }

    fn combined_stress_break(value_tensile: f64, value_shear: Option<f64>) -> BreakageConfig {
        BreakageConfig::CombinedStress {
            tensile: breakage::ThresholdDistribution::Constant { value: value_tensile },
            shear:   value_shear.map(|v| breakage::ThresholdDistribution::Constant { value: v }),
        }
    }

    #[test]
    fn bond_breaks_on_tensile_stress() {
        // Large stretch + moderate σ_max → should break.
        let cfg = BondConfig {
            normal_stiffness: 1e10,    // huge → easy to exceed σ_max
            breakage: Some(combined_stress_break(1e5, None)),
            ..BondConfig::default()
        };
        let mut app = build_pair_app(0.001, 0.003, cfg); // 50% stretch
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let bonds = registry.expect::<BondStore>("test");
        assert_eq!(bonds.bonds[0].len(), 0, "bond should break on tensile stress");
        assert_eq!(bonds.bonds[1].len(), 0);

        let metrics = app.get_resource_ref::<BondMetrics>().unwrap();
        assert_eq!(metrics.total_bonds_broken, 1);
    }

    #[test]
    fn bond_no_break_below_tensile_stress() {
        let cfg = BondConfig {
            normal_stiffness: 1e7,                          // modest stiffness
            breakage: Some(combined_stress_break(1e12, None)), // huge → never breaks
            ..BondConfig::default()
        };
        let mut app = build_pair_app(0.001, 0.0021, cfg); // 5% stretch
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let bonds = registry.expect::<BondStore>("test");
        assert_eq!(bonds.bonds[0].len(), 1);
    }

    #[test]
    fn bond_breaks_on_shear_stress() {
        // Large tangential velocity integrated one step → large Δs → large F_t → large τ.
        let cfg = BondConfig {
            shear_stiffness: 1e10,
            breakage: Some(combined_stress_break(f64::INFINITY, Some(1e5))),
            ..BondConfig::default()
        };
        let mut app = build_pair_app_with(0.001, 0.002, cfg, [0.0, 100.0, 0.0], [0.0; 3]);
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let bonds = registry.expect::<BondStore>("test");
        assert_eq!(bonds.bonds[0].len(), 0, "bond should break on shear stress");
    }

    #[test]
    fn bond_history_pack_unpack_round_trip() {
        let mut store = BondHistoryStore::new();
        store.history.push(vec![
            BondHistoryEntry {
                partner_tag: 5,
                delta_t: [0.1, 0.2, 0.3],
                delta_theta: [0.4, 0.5, 0.6],
                thresholds: [7.0, 8.0, 9.0, 10.0],
                theta_p_bend: [0.15, 0.25, 0.35],
                eps_p_axial: 0.0125,
                theta_max_bend: 0.42,
                eps_max_axial: 0.018,
            },
            BondHistoryEntry {
                partner_tag: 10,
                delta_t: [1.0, 2.0, 3.0],
                delta_theta: [4.0, 5.0, 6.0],
                thresholds: [11.0, 12.0, 13.0, 14.0],
                theta_p_bend: [1.5, 2.5, 3.5],
                eps_p_axial: -0.005,
                theta_max_bend: 3.7,
                eps_max_axial: 0.05,
            },
        ]);

        let mut buf = Vec::new();
        store.pack(0, &mut buf);
        assert_eq!(buf.len(), 1 + 17 * 2);

        let mut store2 = BondHistoryStore::new();
        let consumed = store2.unpack(&buf);
        assert_eq!(consumed, buf.len());
        assert_eq!(store2.history[0].len(), 2);
        assert_eq!(store2.history[0][0].partner_tag, 5);
        assert!((store2.history[0][0].delta_t[0] - 0.1).abs() < 1e-15);
        assert!((store2.history[0][1].delta_theta[2] - 6.0).abs() < 1e-15);
        assert!((store2.history[0][0].thresholds[3] - 10.0).abs() < 1e-15);
        assert!((store2.history[0][1].thresholds[0] - 11.0).abs() < 1e-15);
        assert!((store2.history[0][0].theta_p_bend[2] - 0.35).abs() < 1e-15);
        assert!((store2.history[0][1].theta_p_bend[0] - 1.5).abs() < 1e-15);
        assert!((store2.history[0][0].eps_p_axial - 0.0125).abs() < 1e-15);
        assert!((store2.history[0][1].eps_p_axial - (-0.005)).abs() < 1e-15);
        assert!((store2.history[0][0].theta_max_bend - 0.42).abs() < 1e-15);
        assert!((store2.history[0][1].eps_max_axial - 0.05).abs() < 1e-15);
    }

    #[test]
    fn bond_config_deserialization() {
        let toml_str = r#"
youngs_modulus = 1e9
shear_modulus  = 4e8
bond_radius_ratio = 0.8
beta_normal  = 0.05
beta_shear   = 0.05
beta_twist   = 0.05
beta_bending = 0.05
seed = 42

[breakage]
kind    = "combined_stress"
tensile = { kind = "constant", value = 5.0e7 }
shear   = { kind = "constant", value = 3.0e7 }
"#;
        let cfg: BondConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(cfg.youngs_modulus, Some(1e9));
        assert_eq!(cfg.shear_modulus, Some(4e8));
        assert!((cfg.bond_radius_ratio - 0.8).abs() < 1e-12);
        assert_eq!(cfg.beta_normal, 0.05);
        assert_eq!(cfg.seed, 42);
        assert!(matches!(cfg.breakage, Some(BreakageConfig::CombinedStress { .. })));
    }

    #[test]
    fn bond_config_with_file_fields() {
        let toml_str = r#"
normal_stiffness = 1e7
file = "data.lammps"
format = "lammps_data"
"#;
        let cfg: BondConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(cfg.file.as_deref(), Some("data.lammps"));
        assert_eq!(cfg.format.as_deref(), Some("lammps_data"));
    }

    #[test]
    fn extend_ghost_cutoff_respects_max_bond_r0() {
        // Registry with a BondStore carrying a known max r0; Domain starts with
        // a tiny ghost_cutoff. The system should bump it to r0 × multiplier.
        let mut app = App::new();

        let mut bond_store = BondStore::new();
        bond_store.bonds.push(vec![
            BondEntry { partner_tag: 2, bond_type: 0, r0: 0.002 },
        ]);
        bond_store.bonds.push(vec![
            BondEntry { partner_tag: 1, bond_type: 0, r0: 0.002 },
            BondEntry { partner_tag: 3, bond_type: 0, r0: 0.005 }, // max
        ]);
        bond_store.bonds.push(vec![
            BondEntry { partner_tag: 2, bond_type: 0, r0: 0.005 },
        ]);

        let mut registry = AtomDataRegistry::new();
        registry.register(bond_store);

        let mut domain = soil_core::Domain::new();
        domain.ghost_cutoff = 0.001; // deliberately too small

        app.add_resource(registry);
        app.add_resource(domain);
        app.add_resource(BondConfig {
            ghost_cutoff_multiplier: 2.5,
            ..BondConfig::default()
        });
        app.add_resource(CommResource(Box::new(SingleProcessComm::new())));
        app.add_update_system(extend_ghost_cutoff_for_bonds, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let domain = app.get_resource_ref::<soil_core::Domain>().unwrap();
        // max r0 = 0.005, multiplier = 2.5 → required = 0.0125
        assert!(
            (domain.ghost_cutoff - 0.0125).abs() < 1e-12,
            "expected ghost_cutoff ≈ 0.0125, got {}",
            domain.ghost_cutoff
        );
    }

    #[test]
    fn extend_ghost_cutoff_disabled_when_multiplier_zero() {
        let mut app = App::new();

        let mut bond_store = BondStore::new();
        bond_store.bonds.push(vec![BondEntry { partner_tag: 2, bond_type: 0, r0: 0.005 }]);

        let mut registry = AtomDataRegistry::new();
        registry.register(bond_store);

        let mut domain = soil_core::Domain::new();
        domain.ghost_cutoff = 0.001;

        app.add_resource(registry);
        app.add_resource(domain);
        app.add_resource(BondConfig { ghost_cutoff_multiplier: 0.0, ..BondConfig::default() });
        app.add_resource(CommResource(Box::new(SingleProcessComm::new())));
        app.add_update_system(extend_ghost_cutoff_for_bonds, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let domain = app.get_resource_ref::<soil_core::Domain>().unwrap();
        assert_eq!(domain.ghost_cutoff, 0.001);
    }

    #[test]
    fn extend_ghost_cutoff_never_shrinks() {
        let mut app = App::new();

        let mut bond_store = BondStore::new();
        bond_store.bonds.push(vec![BondEntry { partner_tag: 2, bond_type: 0, r0: 0.002 }]);

        let mut registry = AtomDataRegistry::new();
        registry.register(bond_store);

        let mut domain = soil_core::Domain::new();
        domain.ghost_cutoff = 0.05; // already larger than 0.002 × 2.5 = 0.005

        app.add_resource(registry);
        app.add_resource(domain);
        app.add_resource(BondConfig { ghost_cutoff_multiplier: 2.5, ..BondConfig::default() });
        app.add_resource(CommResource(Box::new(SingleProcessComm::new())));
        app.add_update_system(extend_ghost_cutoff_for_bonds, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let domain = app.get_resource_ref::<soil_core::Domain>().unwrap();
        assert_eq!(domain.ghost_cutoff, 0.05, "must not shrink an already-larger cutoff");
    }

    #[test]
    fn bond_config_without_file_fields() {
        let toml_str = r#"
auto_bond = true
normal_stiffness = 1e7
"#;
        let cfg: BondConfig = toml::from_str(toml_str).unwrap();
        assert!(cfg.file.is_none());
        assert!(cfg.auto_bond);
    }
}
