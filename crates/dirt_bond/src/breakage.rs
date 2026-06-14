//! Bond breakage criteria.
//!
//! This module defines a [`BreakageCriterion`] trait and a menu of concrete
//! implementations. Every criterion takes the same three snapshots
//! ([`BondGeom`], [`BondLoads`], [`BondKinematics`]) plus a per-bond
//! [`BondThresholds`] record (sampled once at bond creation from a
//! [`ThresholdDistribution`]) and returns `Some(BreakMode)` when the bond has
//! failed.
//!
//! ## Criterion menu
//!
//! Three families × three value-types, eight combinations plus an
//! `Unbreakable` no-op.
//!
//! | Criterion              | Tensile branch                                  | Shear branch                                     |
//! |------------------------|-------------------------------------------------|--------------------------------------------------|
//! | `AxialForce`           | `|F_n| > t[0]`                                  | `|F_t| > t[1]`                                   |
//! | `AxialStress`          | `|F_n|/A > t[0]`                                | `|F_t|/A > t[1]`                                 |
//! | `AxialStrain`          | `ε_axial > t[0]`                                | `|Δs|/L > t[1]`                                  |
//! | `CombinedStress`       | `|F_n|/A + r_b·|M_bend|/I > t[0]`               | `|F_t|/A + r_b·|M_tor|/I_p > t[1]`               |
//! | `CombinedStrain`       | `ε_axial + r_b·|κ_bend| > t[0]`                 | `|Δs|/L + r_b·|κ_tor| > t[1]`                    |
//! | `InteractionLinearForce`  | single envelope: `Σ |L_i|/t[i] ≥ 1` over 4 force/moment channels                  |
//! | `InteractionLinearStress` | single envelope: `Σ stress_i/t[i] ≥ 1` over 4 stress channels                     |
//! | `InteractionLinearStrain` | single envelope: `Σ strain_i/t[i] ≥ 1` over 4 strain channels                     |
//!
//! ## Channel indexing
//!
//! For [`BondThresholds`] the four entries `t[0..4]` are:
//!
//! | Index | Axial / Combined family    | InteractionLinear family       |
//! |-------|----------------------------|--------------------------------|
//! | `t[0]` | tensile threshold         | axial channel threshold        |
//! | `t[1]` | shear threshold (or `∞`)  | shear channel threshold        |
//! | `t[2]` | unused                    | bending channel threshold      |
//! | `t[3]` | unused                    | twist channel threshold        |
//!
//! ## Sources
//!
//! * `CombinedStress` is the Potyondy-Cundall (2004) / Brown (2014) / Guo (2017)
//!   extreme-fibre beam-stress criterion.
//! * `CombinedStrain` is the strain analog derived in
//!   `lit_review_fiber_dem_breakage/BPM_BREAKAGE_MIGRATION.md` Eq. 1.7-1.
//! * `InteractionLinear*` follows the linear damage-accumulation envelope of
//!   LAMMPS `bond_style bpm/rotational` (Clemmer, Monti & Lechman,
//!   *Soft Matter* **20**, 1702, 2024 — `B = Σ |X_i|/X_i,c ≥ 1`).
//! * The Weibull size-effect scaling and per-bond inverse-CDF sampler follow
//!   Weibull (1939) and the migration doc Eqs. 1.5-1, 1.6-1.

use serde::Deserialize;

// ── MPI-stable per-bond uniform sampling ────────────────────────────────────

/// Draw the four uniform samples `u ∈ (0, 1)` for a single bond, **deterministically
/// from the bond's tag pair and the global seed alone** — no rank-local RNG state.
///
/// This is the MPI-stable replacement for sampling from a per-rank `StdRng`:
/// the same bond (identified by an unordered tag pair) always gets the same
/// four samples regardless of which rank ends up owning it, in which order
/// bonds are visited, or which MPI decomposition was chosen. Two ranks
/// independently visiting the same bond pair will compute identical
/// thresholds, so a simulation re-run with a different decomposition is
/// bit-reproducible in its breakage pattern.
///
/// Implementation: canonicalise the tag pair (lo, hi), mix into a 64-bit
/// seed via SplitMix64, drive a [`SmallRng`](rand::rngs::SmallRng) from that
/// seed, draw four `u ∈ [1e-15, 1 - 1e-15]` samples. SmallRng is overkill
/// for four numbers but matches the rest of the crate's RNG conventions and
/// is fast at this point in the run (bond creation only).
pub fn per_bond_uniform_samples(tag_a: u32, tag_b: u32, seed: u64) -> [f64; 4] {
    use rand::{Rng, SeedableRng};
    use rand::rngs::SmallRng;

    let (lo, hi) = if tag_a <= tag_b { (tag_a, tag_b) } else { (tag_b, tag_a) };
    // Mix seed, lo, hi into a unique 64-bit bond seed via two SplitMix64 rounds.
    // The bit-mixing is good enough that adjacent tag pairs land in
    // uncorrelated parts of the [0,1) domain, which is what the Weibull
    // sampler downstream needs to behave like independent draws.
    let mut bond_seed = seed;
    bond_seed = splitmix64(bond_seed ^ (lo as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    bond_seed = splitmix64(bond_seed ^ (hi as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9));

    let mut rng = SmallRng::seed_from_u64(bond_seed);
    [
        rng.random_range(1.0e-15..1.0 - 1.0e-15),
        rng.random_range(1.0e-15..1.0 - 1.0e-15),
        rng.random_range(1.0e-15..1.0 - 1.0e-15),
        rng.random_range(1.0e-15..1.0 - 1.0e-15),
    ]
}

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

// ── Geometry / loads / kinematics passed to every criterion ─────────────────

/// Snapshot of bond geometry needed by any criterion.
#[derive(Clone, Copy, Debug)]
pub struct BondGeom {
    /// Bond radius `r_b` (m).
    pub r_b: f64,
    /// Cross-sectional area `A = π r_b²` (m²).
    pub area: f64,
    /// Second moment of area for bending `I = π r_b⁴ / 4` (m⁴).
    pub iben: f64,
    /// Polar second moment `J = π r_b⁴ / 2` (m⁴) — for torsion.
    pub jpol: f64,
    /// Initial bond length `L₀` (m).
    pub l0: f64,
}

/// Bond loads (forces and moments) at the current step.
#[derive(Clone, Copy, Debug)]
pub struct BondLoads {
    /// Signed axial force (+ tension, − compression).
    pub f_n: f64,
    /// Magnitude of shear force vector `|F_t|`.
    pub f_t_mag: f64,
    /// Magnitude of bending moment vector `|M_bend|`.
    pub m_bend_mag: f64,
    /// Magnitude of twist (torsion) moment along bond axis `|M_tor|`.
    pub m_tor_mag: f64,
}

/// Bond kinematics — purely positional / orientational, no `E` entering.
#[derive(Clone, Copy, Debug)]
pub struct BondKinematics {
    /// Axial strain `(L − L₀)/L₀`.
    pub eps_axial: f64,
    /// Shear strain magnitude `|Δs|/L`.
    pub gamma_shear: f64,
    /// Bending-curvature magnitude `|Δθ_bend|/L`.
    pub kappa_bend: f64,
    /// Twist-rate magnitude `|Δθ_tor|/L`.
    pub kappa_tor: f64,
}

/// Per-bond failure thresholds. Up to four entries; meaning is criterion-dependent.
#[derive(Clone, Copy, Debug, Default)]
pub struct BondThresholds {
    /// Threshold values; entries unused by a given criterion are ignored.
    pub t: [f64; 4],
}

/// Which failure mode tripped.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BreakMode {
    /// Tensile (axial / bending) branch of a two-branch criterion.
    Tensile,
    /// Shear (tangential / torsion) branch of a two-branch criterion.
    Shear,
    /// Combined-envelope failure (`InteractionLinear*`).
    Interaction,
}

// ── Threshold distribution (Constant or Weibull) ────────────────────────────

/// Distribution from which a per-bond threshold is sampled at bond creation.
///
/// `Constant` returns the same value for every bond. `Weibull` returns a draw
/// from a length-scaled 2-parameter Weibull
///
/// ```text
///     value = (mean / Γ(1 + 1/m))
///           · (L_calib / max(L_bond, L_min))^{1/m}
///           · (−ln(1 − u))^{1/m}
/// ```
///
/// where `u ∈ (0, 1)` is supplied by the caller. The factor `L_calib / L_eff`
/// implements the weakest-link size effect; `L_min` floors the size effect at
/// a material-flaw spacing to avoid unphysical strengthening of very short
/// bonds (see migration doc §2.12).
#[derive(Deserialize, Clone, Debug)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ThresholdDistribution {
    /// Single deterministic value used for every bond.
    Constant {
        /// Threshold value (units depend on criterion: N, Pa, or dimensionless).
        value: f64,
    },
    /// Length-scaled two-parameter Weibull distribution.
    Weibull {
        /// Experimental mean breaking value `\bar{x}_exp` at gauge length `l_calib`.
        mean: f64,
        /// Weibull modulus / shape parameter `m`.
        m: f64,
        /// Gauge length `L_calib` (m) at which `mean` was measured.
        l_calib: f64,
        /// Lower bound on the effective length used in the size-effect factor
        /// (m). Defaults to `0.0` (no floor); set to material flaw spacing to
        /// avoid unphysically strong short bonds.
        #[serde(default)]
        l_min: f64,
    },
    /// Crack-band-style length-rescaled deterministic threshold (Bažant 1976,
    /// Hillerborg-Modéer-Petersson 1976). The threshold's post-yield extent
    /// rescales inversely with bond length so the per-bond plastic + brittle
    /// energy budget × bond length is invariant:
    ///
    /// ```text
    ///     threshold(L_bond) = eps_yield + (value_ref - eps_yield) · L_ref / L_bond
    /// ```
    ///
    /// For breakage criteria that operate on total strain (`axial_strain`,
    /// `combined_strain`, `interaction_linear_strain`), this is the correct
    /// regularization that flattens the apparent strength-vs-mesh-refinement
    /// curve. For non-strain criteria (force, stress) set `eps_yield = 0`
    /// to scale the entire threshold value.
    CrackBand {
        /// Threshold value at the calibration bond length `l_ref`.
        value_ref: f64,
        /// Reference bond length (m) at which `value_ref` was calibrated.
        l_ref: f64,
        /// Elastic yield strain (or other "elastic anchor"). The portion of
        /// the threshold **below** `eps_yield` is treated as a true material
        /// constant; only the part **above** rescales with bond length. For
        /// non-strain criteria, set this to `0.0`. Default: `0.0`.
        #[serde(default)]
        eps_yield: f64,
        /// Lower bound on the bond length used in the rescaling (m). Avoids
        /// unphysical strengthening of very short bonds. Default: `0.0`.
        #[serde(default)]
        l_min: f64,
    },
}

impl ThresholdDistribution {
    /// Draw one threshold for a bond of length `l_bond` given a uniform sample `u ∈ (0,1)`.
    pub fn sample(&self, l_bond: f64, u: f64) -> f64 {
        match *self {
            Self::Constant { value } => value,
            Self::Weibull { mean, m, l_calib, l_min } => {
                let l_eff = l_bond.max(l_min).max(f64::MIN_POSITIVE);
                let size_factor = (l_calib / l_eff).powf(1.0 / m);
                let u_clamped = u.clamp(1e-15, 1.0 - 1e-15);
                let scale = mean / gamma_lanczos(1.0 + 1.0 / m);
                scale * size_factor * (-((1.0 - u_clamped).ln())).powf(1.0 / m)
            }
            Self::CrackBand { value_ref, l_ref, eps_yield, l_min } => {
                let l_eff = l_bond.max(l_min).max(f64::MIN_POSITIVE);
                eps_yield + (value_ref - eps_yield) * (l_ref / l_eff)
            }
        }
    }
}

// ── Trait + implementations ─────────────────────────────────────────────────

/// A bond breakage criterion. Implementations are pure functions of the
/// criterion's configuration; per-bond state lives in [`BondThresholds`] and
/// the run-time snapshots ([`BondGeom`], [`BondLoads`], [`BondKinematics`]).
pub trait BreakageCriterion: Send + Sync + std::fmt::Debug {
    /// Number of distributions the criterion samples per bond (0..=4).
    fn num_thresholds(&self) -> usize;

    /// Sample one bond's thresholds. `u[i]` must be a uniform draw in (0,1).
    fn sample(&self, l_bond: f64, u: [f64; 4]) -> BondThresholds;

    /// Check the bond against its thresholds. Returns the failure mode, if any.
    fn check(
        &self,
        geom: &BondGeom,
        loads: &BondLoads,
        kin: &BondKinematics,
        thr: &BondThresholds,
    ) -> Option<BreakMode>;
}

/// A bond that never breaks. Useful as a default and for purely-elastic tests.
#[derive(Clone, Debug, Default)]
pub struct Unbreakable;

impl BreakageCriterion for Unbreakable {
    fn num_thresholds(&self) -> usize { 0 }
    fn sample(&self, _: f64, _: [f64; 4]) -> BondThresholds { BondThresholds::default() }
    fn check(&self, _: &BondGeom, _: &BondLoads, _: &BondKinematics, _: &BondThresholds) -> Option<BreakMode> { None }
}

// Helper macro: implement the two-branch criteria (Axial*, Combined*) where
// the tensile and shear value-extractors differ but the structure is identical.
macro_rules! impl_two_branch {
    ($name:ident, $tensile_expr:expr, $shear_expr:expr) => {
        impl BreakageCriterion for $name {
            fn num_thresholds(&self) -> usize { 2 }
            fn sample(&self, l_bond: f64, u: [f64; 4]) -> BondThresholds {
                let t0 = self.tensile.sample(l_bond, u[0]);
                let t1 = match &self.shear {
                    Some(d) => d.sample(l_bond, u[1]),
                    None => f64::INFINITY,
                };
                BondThresholds { t: [t0, t1, 0.0, 0.0] }
            }
            fn check(
                &self,
                geom: &BondGeom,
                loads: &BondLoads,
                kin: &BondKinematics,
                thr: &BondThresholds,
            ) -> Option<BreakMode> {
                let tensile_val: f64 = $tensile_expr(geom, loads, kin);
                if tensile_val > thr.t[0] { return Some(BreakMode::Tensile); }
                let shear_val: f64 = $shear_expr(geom, loads, kin);
                if shear_val > thr.t[1] { return Some(BreakMode::Shear); }
                None
            }
        }
    };
}

/// Tensile/shear failure when the raw axial / shear force exceeds its limit.
/// Bending and torsion moments are ignored.
#[derive(Clone, Debug)]
pub struct AxialForce {
    /// Tensile threshold distribution (units: N).
    pub tensile: ThresholdDistribution,
    /// Optional shear-force threshold distribution. `None` disables shear failure.
    pub shear: Option<ThresholdDistribution>,
}
impl_two_branch!(
    AxialForce,
    |_g: &BondGeom, l: &BondLoads, _k: &BondKinematics| l.f_n.max(0.0),
    |_g: &BondGeom, l: &BondLoads, _k: &BondKinematics| l.f_t_mag
);

/// Tensile/shear failure based on the axial component of stress (force / area).
/// Bending and torsion contributions are ignored.
#[derive(Clone, Debug)]
pub struct AxialStress {
    /// Tensile-stress threshold distribution (units: Pa).
    pub tensile: ThresholdDistribution,
    /// Optional shear-stress threshold distribution (units: Pa).
    pub shear: Option<ThresholdDistribution>,
}
impl_two_branch!(
    AxialStress,
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics|
        if g.area > 0.0 { l.f_n.max(0.0) / g.area } else { 0.0 },
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics|
        if g.area > 0.0 { l.f_t_mag / g.area } else { 0.0 }
);

/// Tensile/shear failure based on kinematic strain only — no moment / curvature.
/// `AxialStrain` is the criterion used by Clemmer & Robbins PRL 2022 and the
/// LAMMPS `bpm/spring` bond style.
#[derive(Clone, Debug)]
pub struct AxialStrain {
    /// Tensile-strain threshold (dimensionless).
    pub tensile: ThresholdDistribution,
    /// Optional shear-strain threshold (dimensionless).
    pub shear: Option<ThresholdDistribution>,
}
impl_two_branch!(
    AxialStrain,
    |_g: &BondGeom, _l: &BondLoads, k: &BondKinematics| k.eps_axial.max(0.0),
    |_g: &BondGeom, _l: &BondLoads, k: &BondKinematics| k.gamma_shear
);

/// Extreme-fibre beam-stress criterion: tensile failure sums axial-stress and
/// bending-stress at the outer fibre; shear failure sums shear-stress and
/// torsion-stress at the same point. This is the Potyondy & Cundall (2004) /
/// Brown (2014) / Guo (2017 Eq. 16-17) form. Matches the current `dirt_bond`
/// behaviour with `sigma_max` / `tau_max` set to constant thresholds.
#[derive(Clone, Debug)]
pub struct CombinedStress {
    /// Tensile-stress threshold (Pa) — applied to `|F_n|/A + r_b·|M_b|/I`.
    pub tensile: ThresholdDistribution,
    /// Optional shear-stress threshold (Pa) — applied to `|F_t|/A + r_b·|M_tor|/I_p`.
    pub shear: Option<ThresholdDistribution>,
}
impl_two_branch!(
    CombinedStress,
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics| {
        let axial = if g.area > 0.0 { l.f_n.max(0.0) / g.area } else { 0.0 };
        let bend  = if g.iben > 0.0 { g.r_b * l.m_bend_mag / g.iben } else { 0.0 };
        axial + bend
    },
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics| {
        let shear = if g.area > 0.0 { l.f_t_mag / g.area } else { 0.0 };
        let tor   = if g.jpol > 0.0 { g.r_b * l.m_tor_mag / g.jpol } else { 0.0 };
        shear + tor
    }
);

/// Extreme-fibre strain criterion (migration doc Eq. 1.7-1). Tensile branch
/// sums axial strain and bending strain `ε_axial + r_b·|κ_bend|`; shear branch
/// sums shear strain and torsion-rate strain `|Δs|/L + r_b·|κ_tor|`. Both
/// branches use **purely kinematic** quantities — no `E`, no force — so the
/// criterion behaves correctly past plastic yield.
#[derive(Clone, Debug)]
pub struct CombinedStrain {
    /// Tensile-strain threshold (dimensionless) — applied to `ε_axial + r_b·|κ_bend|`.
    pub tensile: ThresholdDistribution,
    /// Optional shear-strain threshold (dimensionless) — applied to `|Δs|/L + r_b·|κ_tor|`.
    pub shear: Option<ThresholdDistribution>,
}
impl_two_branch!(
    CombinedStrain,
    |g: &BondGeom, _l: &BondLoads, k: &BondKinematics|
        k.eps_axial.max(0.0) + g.r_b * k.kappa_bend,
    |g: &BondGeom, _l: &BondLoads, k: &BondKinematics|
        k.gamma_shear + g.r_b * k.kappa_tor
);

// ── InteractionLinear family ────────────────────────────────────────────────

/// Linear damage-accumulation criterion in force/moment space. Channels are
/// (axial-force, shear-force, bending-moment, twist-moment). Failure when
///
/// ```text
///     |F_n|/t[0]  +  |F_t|/t[1]  +  |M_bend|/t[2]  +  |M_tor|/t[3]  ≥  1
/// ```
///
/// Set any channel's distribution to `None` to leave that channel out of the
/// sum. Matches the LAMMPS `bond_style bpm/rotational` form (Clemmer-Monti-
/// Lechman 2024) when all four channels are active.
#[derive(Clone, Debug)]
pub struct InteractionLinearForce {
    /// Axial-force-threshold distribution (N).
    pub axial:   Option<ThresholdDistribution>,
    /// Shear-force-threshold distribution (N).
    pub shear:   Option<ThresholdDistribution>,
    /// Bending-moment-threshold distribution (N·m).
    pub bending: Option<ThresholdDistribution>,
    /// Twist-moment-threshold distribution (N·m).
    pub twist:   Option<ThresholdDistribution>,
}

/// Linear damage-accumulation criterion in stress space. Channels are
/// (axial-stress `|F_n|/A`, shear-stress `|F_t|/A`, bending-stress at extreme
/// fibre `r_b·|M_bend|/I`, torsion-stress at extreme fibre `r_b·|M_tor|/I_p`).
/// Failure when their normalized sum reaches 1.
#[derive(Clone, Debug)]
pub struct InteractionLinearStress {
    /// Axial-stress-threshold distribution (Pa).
    pub axial:   Option<ThresholdDistribution>,
    /// Shear-stress-threshold distribution (Pa).
    pub shear:   Option<ThresholdDistribution>,
    /// Bending-stress-threshold distribution (Pa).
    pub bending: Option<ThresholdDistribution>,
    /// Torsion-stress-threshold distribution (Pa).
    pub twist:   Option<ThresholdDistribution>,
}

/// Linear damage-accumulation criterion in strain space. Channels are
/// (axial-strain `ε_axial`, shear-strain `|Δs|/L`, bending-strain at extreme
/// fibre `r_b·|κ_bend|`, twist-rate strain at extreme fibre `r_b·|κ_tor|`).
/// Failure when their normalized sum reaches 1.
#[derive(Clone, Debug)]
pub struct InteractionLinearStrain {
    /// Axial-strain-threshold distribution (dimensionless).
    pub axial:   Option<ThresholdDistribution>,
    /// Shear-strain-threshold distribution (dimensionless).
    pub shear:   Option<ThresholdDistribution>,
    /// Bending-strain-threshold distribution (dimensionless).
    pub bending: Option<ThresholdDistribution>,
    /// Twist-strain-threshold distribution (dimensionless).
    pub twist:   Option<ThresholdDistribution>,
}

// Helper macro: share the sampling and check logic across the three
// InteractionLinear variants. The only thing that differs between them is the
// per-channel value extractor passed in.
macro_rules! impl_interaction_linear {
    ($name:ident, $axial:expr, $shear:expr, $bending:expr, $twist:expr) => {
        impl BreakageCriterion for $name {
            fn num_thresholds(&self) -> usize { 4 }
            fn sample(&self, l_bond: f64, u: [f64; 4]) -> BondThresholds {
                let s = |d: &Option<ThresholdDistribution>, ui: f64|
                    d.as_ref().map(|x| x.sample(l_bond, ui)).unwrap_or(f64::INFINITY);
                BondThresholds { t: [
                    s(&self.axial,   u[0]),
                    s(&self.shear,   u[1]),
                    s(&self.bending, u[2]),
                    s(&self.twist,   u[3]),
                ] }
            }
            fn check(
                &self,
                geom: &BondGeom,
                loads: &BondLoads,
                kin: &BondKinematics,
                thr: &BondThresholds,
            ) -> Option<BreakMode> {
                let mut sum = 0.0;
                let v_axial:   f64 = $axial(geom, loads, kin);
                let v_shear:   f64 = $shear(geom, loads, kin);
                let v_bending: f64 = $bending(geom, loads, kin);
                let v_twist:   f64 = $twist(geom, loads, kin);
                // Channels with infinite thresholds (i.e. `None` in config) drop
                // out cleanly since x/∞ = 0.
                sum += v_axial   / thr.t[0];
                sum += v_shear   / thr.t[1];
                sum += v_bending / thr.t[2];
                sum += v_twist   / thr.t[3];
                if sum >= 1.0 { Some(BreakMode::Interaction) } else { None }
            }
        }
    };
}

impl_interaction_linear!(
    InteractionLinearForce,
    |_g: &BondGeom, l: &BondLoads, _k: &BondKinematics| l.f_n.max(0.0).abs(),
    |_g: &BondGeom, l: &BondLoads, _k: &BondKinematics| l.f_t_mag,
    |_g: &BondGeom, l: &BondLoads, _k: &BondKinematics| l.m_bend_mag,
    |_g: &BondGeom, l: &BondLoads, _k: &BondKinematics| l.m_tor_mag
);

impl_interaction_linear!(
    InteractionLinearStress,
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics|
        if g.area > 0.0 { l.f_n.max(0.0) / g.area } else { 0.0 },
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics|
        if g.area > 0.0 { l.f_t_mag / g.area } else { 0.0 },
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics|
        if g.iben > 0.0 { g.r_b * l.m_bend_mag / g.iben } else { 0.0 },
    |g: &BondGeom, l: &BondLoads, _k: &BondKinematics|
        if g.jpol > 0.0 { g.r_b * l.m_tor_mag / g.jpol } else { 0.0 }
);

impl_interaction_linear!(
    InteractionLinearStrain,
    |_g: &BondGeom, _l: &BondLoads, k: &BondKinematics| k.eps_axial.max(0.0),
    |_g: &BondGeom, _l: &BondLoads, k: &BondKinematics| k.gamma_shear,
    |g: &BondGeom, _l: &BondLoads, k: &BondKinematics| g.r_b * k.kappa_bend,
    |g: &BondGeom, _l: &BondLoads, k: &BondKinematics| g.r_b * k.kappa_tor
);

// ── Config enum (TOML-deserializable) ───────────────────────────────────────

/// Deserializable configuration for the `[bonds.breakage]` table. Each variant
/// matches one of the criteria above; the variant tag in TOML is the
/// `snake_case` form of the criterion name.
///
/// ```toml
/// # Stress-based combined (Guo / Potyondy-Cundall — the dirt_bond default):
/// [bonds.breakage]
/// kind = "combined_stress"
/// tensile = { kind = "constant", value = 5.0e7 }
/// shear   = { kind = "constant", value = 3.0e7 }
///
/// # Strain-based combined with a Weibull tensile threshold:
/// [bonds.breakage]
/// kind = "combined_strain"
/// tensile = { kind = "weibull", mean = 0.018, m = 5.3, l_calib = 0.020 }
/// shear   = { kind = "constant", value = 0.05 }
///
/// # Clemmer bpm/rotational — linear damage envelope over four channels:
/// [bonds.breakage]
/// kind  = "interaction_linear_stress"
/// axial   = { kind = "constant", value = 5.0e7 }
/// shear   = { kind = "constant", value = 3.0e7 }
/// bending = { kind = "constant", value = 5.0e7 }
/// twist   = { kind = "constant", value = 3.0e7 }
/// ```
#[derive(Deserialize, Clone, Debug)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum BreakageConfig {
    /// Never breaks any bond.
    Unbreakable,
    /// See [`AxialForce`].
    AxialForce {
        tensile: ThresholdDistribution,
        #[serde(default)]
        shear: Option<ThresholdDistribution>,
    },
    /// See [`AxialStress`].
    AxialStress {
        tensile: ThresholdDistribution,
        #[serde(default)]
        shear: Option<ThresholdDistribution>,
    },
    /// See [`AxialStrain`].
    AxialStrain {
        tensile: ThresholdDistribution,
        #[serde(default)]
        shear: Option<ThresholdDistribution>,
    },
    /// See [`CombinedStress`].
    CombinedStress {
        tensile: ThresholdDistribution,
        #[serde(default)]
        shear: Option<ThresholdDistribution>,
    },
    /// See [`CombinedStrain`].
    CombinedStrain {
        tensile: ThresholdDistribution,
        #[serde(default)]
        shear: Option<ThresholdDistribution>,
    },
    /// See [`InteractionLinearForce`].
    InteractionLinearForce {
        #[serde(default)] axial:   Option<ThresholdDistribution>,
        #[serde(default)] shear:   Option<ThresholdDistribution>,
        #[serde(default)] bending: Option<ThresholdDistribution>,
        #[serde(default)] twist:   Option<ThresholdDistribution>,
    },
    /// See [`InteractionLinearStress`].
    InteractionLinearStress {
        #[serde(default)] axial:   Option<ThresholdDistribution>,
        #[serde(default)] shear:   Option<ThresholdDistribution>,
        #[serde(default)] bending: Option<ThresholdDistribution>,
        #[serde(default)] twist:   Option<ThresholdDistribution>,
    },
    /// See [`InteractionLinearStrain`].
    InteractionLinearStrain {
        #[serde(default)] axial:   Option<ThresholdDistribution>,
        #[serde(default)] shear:   Option<ThresholdDistribution>,
        #[serde(default)] bending: Option<ThresholdDistribution>,
        #[serde(default)] twist:   Option<ThresholdDistribution>,
    },
}

impl BreakageConfig {
    /// Construct the run-time criterion trait object.
    pub fn build(&self) -> Box<dyn BreakageCriterion> {
        match self {
            Self::Unbreakable => Box::new(Unbreakable),
            Self::AxialForce { tensile, shear } =>
                Box::new(AxialForce { tensile: tensile.clone(), shear: shear.clone() }),
            Self::AxialStress { tensile, shear } =>
                Box::new(AxialStress { tensile: tensile.clone(), shear: shear.clone() }),
            Self::AxialStrain { tensile, shear } =>
                Box::new(AxialStrain { tensile: tensile.clone(), shear: shear.clone() }),
            Self::CombinedStress { tensile, shear } =>
                Box::new(CombinedStress { tensile: tensile.clone(), shear: shear.clone() }),
            Self::CombinedStrain { tensile, shear } =>
                Box::new(CombinedStrain { tensile: tensile.clone(), shear: shear.clone() }),
            Self::InteractionLinearForce { axial, shear, bending, twist } =>
                Box::new(InteractionLinearForce {
                    axial: axial.clone(), shear: shear.clone(),
                    bending: bending.clone(), twist: twist.clone(),
                }),
            Self::InteractionLinearStress { axial, shear, bending, twist } =>
                Box::new(InteractionLinearStress {
                    axial: axial.clone(), shear: shear.clone(),
                    bending: bending.clone(), twist: twist.clone(),
                }),
            Self::InteractionLinearStrain { axial, shear, bending, twist } =>
                Box::new(InteractionLinearStrain {
                    axial: axial.clone(), shear: shear.clone(),
                    bending: bending.clone(), twist: twist.clone(),
                }),
        }
    }
}

// ── Gamma function (Lanczos g=7, n=9) ───────────────────────────────────────

/// `Γ(z)` for `z > 0` via Lanczos approximation (~15-digit accuracy).
fn gamma_lanczos(z: f64) -> f64 {
    const G: f64 = 7.0;
    const P: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_59,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_571_6e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z < 0.5 {
        std::f64::consts::PI / ((std::f64::consts::PI * z).sin() * gamma_lanczos(1.0 - z))
    } else {
        let z = z - 1.0;
        let mut x = P[0];
        for (i, &p) in P.iter().enumerate().skip(1) {
            x += p / (z + i as f64);
        }
        let t = z + G + 0.5;
        (2.0 * std::f64::consts::PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x
    }
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Reference geometry: r_b = 1 mm cylindrical bond, L₀ = 2 mm.
    fn geom() -> BondGeom {
        let r_b: f64 = 1.0e-3;
        BondGeom {
            r_b,
            area: std::f64::consts::PI * r_b * r_b,
            iben: 0.25 * std::f64::consts::PI * r_b.powi(4),
            jpol: 0.5  * std::f64::consts::PI * r_b.powi(4),
            l0: 2.0e-3,
        }
    }

    fn zero_loads() -> BondLoads {
        BondLoads { f_n: 0.0, f_t_mag: 0.0, m_bend_mag: 0.0, m_tor_mag: 0.0 }
    }
    fn zero_kin() -> BondKinematics {
        BondKinematics { eps_axial: 0.0, gamma_shear: 0.0, kappa_bend: 0.0, kappa_tor: 0.0 }
    }

    // ── MPI-stable per-bond uniform sampler ─────────────────────────────

    #[test]
    fn per_bond_samples_are_deterministic_in_pair_and_seed() {
        // Same (tag_a, tag_b, seed) → same four samples. Calling twice gives
        // identical output (no internal mutable state).
        let a = per_bond_uniform_samples(3, 17, 42);
        let b = per_bond_uniform_samples(3, 17, 42);
        assert_eq!(a, b);
        // All four samples sit strictly inside (0, 1).
        for u in a {
            assert!(u > 0.0 && u < 1.0, "sample {} outside (0,1)", u);
        }
    }

    #[test]
    fn per_bond_samples_are_canonical_in_tag_order() {
        // (tag_a, tag_b) and (tag_b, tag_a) describe the same bond — sampler
        // must produce identical output regardless of which end is "i" and
        // which is "j". This is what guarantees MPI stability: which atom
        // is local vs ghost depends on decomposition, but the canonical
        // (min, max) tag pair does not.
        let a = per_bond_uniform_samples(7, 99, 12345);
        let b = per_bond_uniform_samples(99, 7, 12345);
        assert_eq!(a, b);
    }

    #[test]
    fn per_bond_samples_differ_across_distinct_pairs() {
        // Adjacent tag pairs land in uncorrelated parts of [0,1) — this is
        // what makes the downstream Weibull draws look independent across
        // neighbouring bonds in a chain.
        let p01 = per_bond_uniform_samples(0, 1, 1);
        let p12 = per_bond_uniform_samples(1, 2, 1);
        let p23 = per_bond_uniform_samples(2, 3, 1);
        assert_ne!(p01, p12);
        assert_ne!(p12, p23);
        assert_ne!(p01, p23);
    }

    #[test]
    fn per_bond_samples_change_with_seed() {
        // Same bond pair under different seeds → different draws (no degenerate
        // "all seeds collapse to one realisation" failure mode).
        let s1 = per_bond_uniform_samples(5, 10, 1);
        let s2 = per_bond_uniform_samples(5, 10, 2);
        let s3 = per_bond_uniform_samples(5, 10, 999_999);
        assert_ne!(s1, s2);
        assert_ne!(s2, s3);
        assert_ne!(s1, s3);
    }

    #[test]
    fn per_bond_samples_roughly_uniform_over_many_pairs() {
        // Pull samples across many bond pairs, bin into 10 deciles, and
        // check each decile gets a reasonable share. This is a smoke test
        // that the SplitMix64 mixing is "uniform enough" — we don't claim
        // statistical rigour, just no obvious bias.
        let n_pairs = 4000;
        let mut bins = [0usize; 10];
        for pair_idx in 0..n_pairs {
            let a = pair_idx as u32;
            let b = (pair_idx + 1) as u32;
            let u = per_bond_uniform_samples(a, b, 7);
            // Use only the first sample to look at the dimension we sampled.
            let decile = (u[0] * 10.0).floor() as usize;
            let decile = decile.min(9);
            bins[decile] += 1;
        }
        let expected = n_pairs / 10;
        let tolerance = (n_pairs as f64 * 0.10) as usize;   // 10 % envelope
        for (i, &count) in bins.iter().enumerate() {
            assert!(
                count.abs_diff(expected) <= tolerance,
                "decile {i}: count {count} too far from expected {expected}"
            );
        }
    }

    #[test]
    fn per_bond_samples_replicate_under_simulated_mpi_partition() {
        // Simulate two MPI decompositions of the same 5-bond chain. In the
        // first decomposition the bonds are visited in order (0,1), (1,2),
        // ..., (4,5). In the second they are visited in the reverse order.
        // The sampler must give the same per-bond `u` regardless of which
        // rank visits the pair first.
        let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)];
        let seed = 0xDEADBEEFu64;
        let forward: Vec<[f64; 4]> =
            bonds.iter().map(|(a, b)| per_bond_uniform_samples(*a, *b, seed)).collect();
        let reverse: Vec<[f64; 4]> =
            bonds.iter().rev().map(|(a, b)| per_bond_uniform_samples(*a, *b, seed)).collect();
        // forward[i] should equal reverse[N-1-i] for the same bond.
        for i in 0..bonds.len() {
            assert_eq!(forward[i], reverse[bonds.len() - 1 - i]);
        }
    }

    #[test]
    fn gamma_reference_values() {
        // Γ(1) = 1, Γ(2) = 1, Γ(1.5) = √π/2, Γ(5) = 24.
        assert!((gamma_lanczos(1.0) - 1.0).abs() < 1e-12);
        assert!((gamma_lanczos(2.0) - 1.0).abs() < 1e-12);
        assert!((gamma_lanczos(1.5) - (std::f64::consts::PI.sqrt() / 2.0)).abs() < 1e-12);
        assert!((gamma_lanczos(5.0) - 24.0).abs() < 1e-9);
    }

    #[test]
    fn constant_distribution_passes_through() {
        let d = ThresholdDistribution::Constant { value: 1.234 };
        for u in [0.01, 0.5, 0.99] {
            assert_eq!(d.sample(1.0e-3, u), 1.234);
        }
    }

    #[test]
    fn weibull_size_effect_reduces_threshold_for_longer_bond() {
        // Two bonds of different lengths; longer bond should sample a smaller
        // threshold at the same `u` (weakest-link size effect).
        let d = ThresholdDistribution::Weibull {
            mean: 1.0e9, m: 5.0, l_calib: 1.0e-3, l_min: 0.0,
        };
        let short = d.sample(1.0e-3, 0.5);
        let long  = d.sample(10.0e-3, 0.5);
        assert!(long < short, "longer bond ({:.3e}) should be weaker than shorter ({:.3e})", long, short);
        // L scales by 10, m = 5, ratio should be 10^{-1/5} ≈ 0.6310.
        let ratio = long / short;
        assert!((ratio - 10f64.powf(-0.2)).abs() < 1e-12);
    }

    #[test]
    fn weibull_mean_recovered_at_uniform_l() {
        // At L_bond = L_calib and u = (1 - 1/e), the size factor is 1 and
        // (-ln(1-u))^{1/m} = 1 — so the sampled value equals the
        // characteristic strength `mean / Γ(1+1/m)`.
        let d = ThresholdDistribution::Weibull {
            mean: 5.0e7, m: 5.0, l_calib: 2.0e-3, l_min: 0.0,
        };
        let u = 1.0 - (-1.0_f64).exp();   // u such that -ln(1-u) = 1
        let v = d.sample(2.0e-3, u);
        let expected = 5.0e7 / gamma_lanczos(1.0 + 1.0 / 5.0);
        assert!((v - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn crack_band_threshold_rescales_with_bond_length() {
        // value_ref = 0.06 at l_ref = 2 mm with eps_yield = 0.02.
        // At l_bond = l_ref/2 = 1 mm: ε_break = 0.02 + 0.04·2  = 0.10.
        // At l_bond = l_ref   = 2 mm: ε_break = 0.02 + 0.04·1  = 0.06.
        // At l_bond = 2·l_ref = 4 mm: ε_break = 0.02 + 0.04·0.5 = 0.04.
        let d = ThresholdDistribution::CrackBand {
            value_ref: 0.06, l_ref: 2.0e-3, eps_yield: 0.02, l_min: 0.0,
        };
        let u = 0.5; // unused for deterministic CrackBand
        assert!((d.sample(1.0e-3, u) - 0.10).abs() < 1e-12);
        assert!((d.sample(2.0e-3, u) - 0.06).abs() < 1e-12);
        assert!((d.sample(4.0e-3, u) - 0.04).abs() < 1e-12);
    }

    #[test]
    fn crack_band_threshold_eps_yield_zero_scales_full_value() {
        // With eps_yield = 0 the whole threshold scales as l_ref / l_bond,
        // useful for force / stress criteria where there's no elastic anchor.
        let d = ThresholdDistribution::CrackBand {
            value_ref: 1.0e8, l_ref: 1.0e-3, eps_yield: 0.0, l_min: 0.0,
        };
        let u = 0.5;
        assert!((d.sample(0.5e-3, u) - 2.0e8).abs() / 2.0e8 < 1e-12);
        assert!((d.sample(2.0e-3, u) - 5.0e7).abs() / 5.0e7 < 1e-12);
    }

    #[test]
    fn crack_band_threshold_l_min_floor() {
        let d = ThresholdDistribution::CrackBand {
            value_ref: 0.06, l_ref: 2.0e-3, eps_yield: 0.02, l_min: 1.0e-3,
        };
        let u = 0.5;
        // Below the floor, the threshold is clamped to the at-floor value.
        assert!((d.sample(1.0e-9, u) - d.sample(1.0e-3, u)).abs() < 1e-12);
    }

    #[test]
    fn weibull_l_min_floor() {
        // A bond shorter than `l_min` should be treated as if it were `l_min`.
        let d = ThresholdDistribution::Weibull {
            mean: 1.0e9, m: 5.0, l_calib: 1.0e-3, l_min: 5.0e-3,
        };
        let very_short = d.sample(1.0e-9, 0.5);
        let at_floor   = d.sample(5.0e-3, 0.5);
        assert_eq!(very_short, at_floor);
    }

    #[test]
    fn unbreakable_never_breaks() {
        let c = Unbreakable;
        let thr = BondThresholds::default();
        let g = geom();
        let l = BondLoads { f_n: 1.0e30, f_t_mag: 1.0e30, m_bend_mag: 1.0e30, m_tor_mag: 1.0e30 };
        let k = BondKinematics { eps_axial: 10.0, gamma_shear: 10.0, kappa_bend: 1.0e6, kappa_tor: 1.0e6 };
        assert!(c.check(&g, &l, &k, &thr).is_none());
    }

    #[test]
    fn axial_force_tensile_break() {
        let c = AxialForce {
            tensile: ThresholdDistribution::Constant { value: 100.0 },
            shear:   None,
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        let l = BondLoads { f_n: 150.0, ..zero_loads() };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), Some(BreakMode::Tensile));
        // Compression of equal magnitude must not trip the tensile branch.
        let l = BondLoads { f_n: -150.0, ..zero_loads() };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), None);
    }

    #[test]
    fn axial_stress_threshold() {
        let c = AxialStress {
            tensile: ThresholdDistribution::Constant { value: 1.0e6 },
            shear:   Some(ThresholdDistribution::Constant { value: 5.0e5 }),
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        // σ = 0.5e6 < 1e6 — no break.
        let l = BondLoads { f_n: 0.5e6 * g.area, ..zero_loads() };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), None);
        // σ = 2e6 > 1e6 — tensile break.
        let l = BondLoads { f_n: 2.0e6 * g.area, ..zero_loads() };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), Some(BreakMode::Tensile));
        // shear stress = 1e6 > 5e5 — shear break.
        let l = BondLoads { f_t_mag: 1.0e6 * g.area, ..zero_loads() };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), Some(BreakMode::Shear));
    }

    #[test]
    fn axial_strain_threshold() {
        let c = AxialStrain {
            tensile: ThresholdDistribution::Constant { value: 0.02 },
            shear:   None,
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        let kin_under = BondKinematics { eps_axial: 0.015, ..zero_kin() };
        let kin_over  = BondKinematics { eps_axial: 0.025, ..zero_kin() };
        assert_eq!(c.check(&g, &zero_loads(), &kin_under, &thr), None);
        assert_eq!(c.check(&g, &zero_loads(), &kin_over,  &thr), Some(BreakMode::Tensile));
    }

    #[test]
    fn combined_stress_matches_guo_eq16() {
        // σ = F_n/A + r_b·|M_bend|/I  — set each contribution at half the
        // threshold and confirm the sum trips while either alone does not.
        let c = CombinedStress {
            tensile: ThresholdDistribution::Constant { value: 1.0e7 },
            shear:   None,
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        // Axial half: F_n = 0.5e7·A → axial-stress = 0.5e7.
        let l_axial_only = BondLoads { f_n: 0.5e7 * g.area, ..zero_loads() };
        assert_eq!(c.check(&g, &l_axial_only, &zero_kin(), &thr), None);
        // Bending half: r_b·M_b/I = 0.5e7  ⇒  M_b = 0.5e7·I/r_b.
        let l_bend_only = BondLoads { m_bend_mag: 0.5e7 * g.iben / g.r_b, ..zero_loads() };
        assert_eq!(c.check(&g, &l_bend_only, &zero_kin(), &thr), None);
        // Both at half — sum = 1e7 ≥ threshold (not strictly >, so tweak above).
        let l_both = BondLoads {
            f_n: 0.51e7 * g.area,
            m_bend_mag: 0.5e7 * g.iben / g.r_b,
            ..zero_loads()
        };
        assert_eq!(c.check(&g, &l_both, &zero_kin(), &thr), Some(BreakMode::Tensile));
    }

    #[test]
    fn combined_strain_matches_migration_doc_eq17() {
        // ε_T,max = ε_axial + r_b·|κ_bend|
        let c = CombinedStrain {
            tensile: ThresholdDistribution::Constant { value: 0.02 },
            shear:   None,
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        let half_axial   = BondKinematics { eps_axial: 0.011, ..zero_kin() };
        let half_bend    = BondKinematics { kappa_bend: 0.011 / g.r_b, ..zero_kin() };
        let combined    = BondKinematics {
            eps_axial: 0.011,
            kappa_bend: 0.011 / g.r_b,
            ..zero_kin()
        };
        assert_eq!(c.check(&g, &zero_loads(), &half_axial, &thr), None);
        assert_eq!(c.check(&g, &zero_loads(), &half_bend,  &thr), None);
        assert_eq!(c.check(&g, &zero_loads(), &combined,   &thr), Some(BreakMode::Tensile));
    }

    #[test]
    fn interaction_linear_force_sums_to_one() {
        // Four equal-strength channels, each loaded to 0.3 of its threshold.
        // Sum = 1.2 ≥ 1 — should break.
        let c = InteractionLinearForce {
            axial:   Some(ThresholdDistribution::Constant { value: 1.0 }),
            shear:   Some(ThresholdDistribution::Constant { value: 1.0 }),
            bending: Some(ThresholdDistribution::Constant { value: 1.0 }),
            twist:   Some(ThresholdDistribution::Constant { value: 1.0 }),
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        let l_below = BondLoads { f_n: 0.2, f_t_mag: 0.2, m_bend_mag: 0.2, m_tor_mag: 0.2 };
        assert_eq!(c.check(&g, &l_below, &zero_kin(), &thr), None);
        let l_above = BondLoads { f_n: 0.3, f_t_mag: 0.3, m_bend_mag: 0.3, m_tor_mag: 0.3 };
        assert_eq!(c.check(&g, &l_above, &zero_kin(), &thr), Some(BreakMode::Interaction));
    }

    #[test]
    fn interaction_linear_disabled_channel_drops_out() {
        // Only the axial channel is active; loading the others to huge values
        // must not contribute to the sum.
        let c = InteractionLinearForce {
            axial:   Some(ThresholdDistribution::Constant { value: 10.0 }),
            shear:   None,
            bending: None,
            twist:   None,
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        let l = BondLoads { f_n: 5.0, f_t_mag: 1.0e6, m_bend_mag: 1.0e6, m_tor_mag: 1.0e6 };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), None);
        let l = BondLoads { f_n: 11.0, ..l };
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), Some(BreakMode::Interaction));
    }

    #[test]
    fn interaction_linear_stress_recovers_clemmer_bpm_rotational() {
        // Each channel at exactly one quarter of its threshold → sum = 1, break.
        let c = InteractionLinearStress {
            axial:   Some(ThresholdDistribution::Constant { value: 4.0e6 }),
            shear:   Some(ThresholdDistribution::Constant { value: 4.0e6 }),
            bending: Some(ThresholdDistribution::Constant { value: 4.0e6 }),
            twist:   Some(ThresholdDistribution::Constant { value: 4.0e6 }),
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        let l = BondLoads {
            f_n:        1.0e6 * g.area,
            f_t_mag:    1.0e6 * g.area,
            m_bend_mag: 1.0e6 * g.iben / g.r_b,
            m_tor_mag:  1.0e6 * g.jpol / g.r_b,
        };
        // sum = 4 · 1e6 / 4e6 = 1.0 → break (sum ≥ 1).
        assert_eq!(c.check(&g, &l, &zero_kin(), &thr), Some(BreakMode::Interaction));
    }

    #[test]
    fn interaction_linear_strain_uses_kinematics() {
        // ε_axial + |Δs|/L + r_b·|κ_bend| + r_b·|κ_tor| against unit thresholds.
        let c = InteractionLinearStrain {
            axial:   Some(ThresholdDistribution::Constant { value: 0.01 }),
            shear:   Some(ThresholdDistribution::Constant { value: 0.01 }),
            bending: Some(ThresholdDistribution::Constant { value: 0.01 }),
            twist:   Some(ThresholdDistribution::Constant { value: 0.01 }),
        };
        let thr = c.sample(geom().l0, [0.5; 4]);
        let g = geom();
        // Each contribution = 0.003 → sum = 0.012 / 0.01 × 4 / 4 = 1.2 > 1.
        let k = BondKinematics {
            eps_axial:   0.003,
            gamma_shear: 0.003,
            kappa_bend:  0.003 / g.r_b,
            kappa_tor:   0.003 / g.r_b,
        };
        assert_eq!(c.check(&g, &zero_loads(), &k, &thr), Some(BreakMode::Interaction));
    }
}
