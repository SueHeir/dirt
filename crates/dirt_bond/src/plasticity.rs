//! Bond plasticity models.
//!
//! Two channels are supported in Phase 1:
//!
//! * **Bending** (`[bonds.plasticity.bending]`) — see [`BendingPlasticityConfig`].
//! * **Axial**   (`[bonds.plasticity.axial]`)   — see [`AxialPlasticityConfig`].
//!
//! Each channel is independently configurable, can be left elastic (omit the
//! sub-table), and shares a common machinery:
//!
//! * a piecewise-linear envelope over **strain magnitude**, with breakpoints
//!   in dimensionless strain and slopes expressed as multipliers on the
//!   elastic stiffness
//! * kinematic-hardening return-map with elastic-slope-`K_e` unloading
//! * per-bond plastic anchor stored in [`super::BondHistoryEntry`]
//!
//! For the bending channel, breakpoints are in **extreme-fibre strain**
//! `ε = r_b · θ_bend / l_b` so they are geometry-independent material
//! properties. For the axial channel, breakpoints are in **axial strain**
//! `ε_axial = (L − L₀)/L₀` directly.
//!
//! ## Update rule (both channels share this structure)
//!
//! Let `x` be the kinematic strain measure (a vector for bending, a scalar
//! for axial), `x_p` the per-bond plastic anchor, and `x_max ≥ 0` the
//! largest `|x|` ever reached on this bond. With `x_e = x − x_p` the elastic
//! strain, `k_e` the elastic stiffness in (force-or-moment)/strain units,
//! and `F_env(·)` the **monotonic kinematic-strain envelope**:
//!
//! 1. update `x_max ← max(x_max, |x|)`
//! 2. trial magnitude       `|F_trial| = k_e · |x_e|`
//! 3. yield ceiling         `F_cap     = F_env(x_max)`
//! 4. if `|F_trial| ≤ F_cap`: elastic — `F = k_e · x_e`, anchor unchanged
//! 5. else (plastic flow):  `F = (F_cap / |F_trial|) · k_e · x_e`,
//!                          anchor advances along `x_e / |x_e|` so that
//!                          `|x_e_new| = F_cap / k_e`.
//!
//! Why `x_max`: the envelope as specified by the user is the
//! stress-vs-kinematic-strain curve under **monotonic loading from origin**.
//! Past first yield the elastic strain grows more slowly than the kinematic
//! strain, so evaluating the envelope in elastic-strain space would give the
//! wrong shape (it would degenerate to "elastic up to first yield + flat
//! cap"). Evaluating at the largest historical kinematic-strain magnitude
//! reproduces the user-specified curve exactly under monotonic loading and
//! provides kinematic-hardening hysteresis under cycling (the cap value
//! freezes at `F_env(x_max)` until the bond is pushed past its prior
//! extreme excursion).
//!
//! ## Bending — Guo 2018 simplified
//!
//! The `BendingPlasticityConfig::GuoBending` variant caps `|M_bend|` at the
//! fully-plastic moment `M^p = (4/3) σ_0 r_b³` from Guo et al. 2018
//! (*Chem. Eng. Sci.* **175**, 118–129, Eq. 31). Single material parameter
//! `σ_0`. The trilinear Guo Eq. 32 shape (elastic → elasto-plastic slope
//! `K_ep = K_e/2` → fully plastic) is recoverable via `Piecewise` with
//! `breakpoint_strains = [σ_0/E_b, ε_p]`, `slope_multipliers = [0.5, 0.0]`.

use serde::Deserialize;

// ── Config ──────────────────────────────────────────────────────────────────

/// Top-level configuration for the `[bonds.plasticity]` table. Each channel
/// (`bending`, `axial`) is independently optional. Omit the field to keep
/// that channel purely elastic.
///
/// ```toml
/// # Bending: Guo elastic-perfectly-plastic. Axial: piecewise hardening.
/// [bonds.plasticity.bending]
/// kind         = "guo_bending"
/// yield_stress = 1.23e8
///
/// [bonds.plasticity.axial]
/// kind                = "piecewise"
/// breakpoint_strains  = [0.01, 0.02, 0.03]
/// slope_multipliers   = [0.5, 0.1, 0.0]
/// ```
#[derive(Deserialize, Clone, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct PlasticityConfig {
    /// Bending-plasticity sub-config. `None` ⇒ bending is purely elastic.
    #[serde(default)]
    pub bending: Option<BendingPlasticityConfig>,
    /// Axial-plasticity sub-config. `None` ⇒ axial is purely elastic.
    #[serde(default)]
    pub axial:   Option<AxialPlasticityConfig>,
}

/// Bending-channel plasticity configuration.
#[derive(Deserialize, Clone, Debug)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum BendingPlasticityConfig {
    /// Guo 2018 elastic-perfectly-plastic: cap `|M_bend|` at `(4/3) σ_0 r_b³`.
    GuoBending {
        /// Material yield stress σ_0 (Pa).
        yield_stress: f64,
    },
    /// Guo 2018 **trilinear** bending envelope (Eq. 32): elastic segment with
    /// slope `K_e = E_b·I/l_b` up to `M^e = σ_0·I/r_b`, then an elasto-plastic
    /// segment with slope `K_ep = K_e/2` up to the fully-plastic cap
    /// `M^p = (4/3)·σ_0·r_b³`, then perfectly plastic at `M^p`.
    ///
    /// At config-build time this variant is expanded to a two-breakpoint
    /// [`Piecewise`] in extreme-fibre-strain space with breakpoints
    /// `[σ_0/E_b, (32 − 3π)/(3π) · σ_0/E_b]` (Guo Eq. 33 in ε-space) and
    /// `slope_multipliers = [0.5, 0.0]`. The runtime envelope and return-map
    /// are the same machinery used by `Piecewise`.
    ///
    /// Requires `[bonds].youngs_modulus` so `ε_e = σ_0/E_b` can be computed.
    GuoTrilinear {
        /// Material yield stress σ_0 (Pa).
        yield_stress: f64,
    },
    /// Piecewise-linear bending envelope.
    Piecewise {
        /// Breakpoints in **extreme-fibre strain** `ε = r_b · θ_bend / l_b`
        /// (dimensionless), strictly ascending. The implicit first segment
        /// is elastic.
        breakpoint_strains: Vec<f64>,
        /// Slope multipliers (per segment past each breakpoint) relative to
        /// `K_e`. Length must equal `breakpoint_strains.len()`.
        slope_multipliers: Vec<f64>,
        /// Optional crack-band-style length calibration (Bažant 1976). When
        /// set, post-yield breakpoints rescale at runtime as
        /// `ε_eff[i] = ε_y + (ε[i] − ε_y) · length_calibration / l_b` so
        /// the per-bond plastic dissipation × bond length stays invariant
        /// under mesh refinement. `breakpoint_strains[0]` (the elastic
        /// yield) is treated as a true material constant and is **not**
        /// rescaled. Omit (or set to `None`) to recover the unregularized
        /// behavior. Units: m.
        #[serde(default)]
        length_calibration: Option<f64>,
    },
}

/// Axial-channel plasticity configuration.
#[derive(Deserialize, Clone, Debug)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum AxialPlasticityConfig {
    /// Piecewise-linear axial envelope.
    Piecewise {
        /// Breakpoints in **axial strain** `ε_axial = (L − L₀)/L₀`
        /// (dimensionless), strictly ascending.
        breakpoint_strains: Vec<f64>,
        /// Slope multipliers (per segment past each breakpoint) relative to
        /// the elastic axial stiffness-per-strain `E_b · A = K_n · L₀`.
        /// Length must equal `breakpoint_strains.len()`.
        slope_multipliers: Vec<f64>,
        /// Optional crack-band-style length calibration (Bažant 1976). See
        /// the bending variant for the rescaling formula and rationale; the
        /// axial channel applies the same rule. Units: m.
        #[serde(default)]
        length_calibration: Option<f64>,
    },
}

// ── Runtime models ──────────────────────────────────────────────────────────

/// Run-time bending-plasticity description, built once at setup from
/// [`BendingPlasticityConfig`].
#[derive(Clone, Debug)]
pub enum BendingPlasticity {
    GuoBending {
        yield_stress: f64,
    },
    Piecewise {
        breakpoint_strains: Vec<f64>,
        slope_multipliers: Vec<f64>,
        length_calibration: Option<f64>,
    },
}

impl BendingPlasticity {
    /// Evaluate the bending-moment envelope `|M_env|` at elastic-strain
    /// magnitude `|θ_e|`. `k_e = E_b·I/l_b` is the elastic bending stiffness
    /// for this bond; `r_b`, `l_b` are the bond geometry.
    pub fn envelope(&self, theta_e_mag: f64, k_e: f64, r_b: f64, l_b: f64) -> f64 {
        match self {
            Self::GuoBending { yield_stress } => {
                let m_p = (4.0 / 3.0) * yield_stress * r_b.powi(3);
                (k_e * theta_e_mag).min(m_p)
            }
            Self::Piecewise { breakpoint_strains, slope_multipliers, length_calibration } => {
                // ε = θ · (r_b / l_b) → θ_break = ε_break · l_b / r_b.
                let scale = if r_b > 0.0 { l_b / r_b } else { 0.0 };
                let length_rescale = length_rescale_factor(*length_calibration, l_b);
                evaluate_piecewise(theta_e_mag, k_e, breakpoint_strains, slope_multipliers, scale, length_rescale)
            }
        }
    }
}

/// Run-time axial-plasticity description, built once at setup from
/// [`AxialPlasticityConfig`].
#[derive(Clone, Debug)]
pub enum AxialPlasticity {
    Piecewise {
        breakpoint_strains: Vec<f64>,
        slope_multipliers: Vec<f64>,
        length_calibration: Option<f64>,
    },
}

impl AxialPlasticity {
    /// Evaluate the axial-force envelope `|F_env|` at elastic-axial-strain
    /// magnitude `|ε_e|`. `k_eff = K_n · L₀ = E_b · A` is the
    /// stiffness-per-strain for this bond. `l_bond = L₀` is the bond rest
    /// length, used only when `length_calibration` is set (crack-band
    /// rescaling); otherwise it is ignored.
    pub fn envelope(&self, eps_e_mag: f64, k_eff: f64, l_bond: f64) -> f64 {
        match self {
            Self::Piecewise { breakpoint_strains, slope_multipliers, length_calibration } => {
                // For axial the strain breakpoints already live in ε-space,
                // so the conversion factor is unity.
                let length_rescale = length_rescale_factor(*length_calibration, l_bond);
                evaluate_piecewise(eps_e_mag, k_eff, breakpoint_strains, slope_multipliers, 1.0, length_rescale)
            }
        }
    }
}

/// Crack-band rescaling factor `α = l_calib / l_bond` (Bažant 1976,
/// Hillerborg-Modéer-Petersson 1976). Returns `1.0` (no rescaling) when
/// `length_calibration` is `None` or either length is non-positive. Callers
/// pass this into [`evaluate_piecewise`]; post-yield strain extents get
/// multiplied by `α` while the elastic yield breakpoint is preserved.
#[inline]
fn length_rescale_factor(length_calibration: Option<f64>, l_bond: f64) -> f64 {
    match length_calibration {
        Some(l_calib) if l_calib > 0.0 && l_bond > 0.0 => l_calib / l_bond,
        _ => 1.0,
    }
}

/// Shared piecewise-linear envelope evaluator. Strain breakpoints are
/// `breakpoint_strains[i]` in dimensionless strain; the conversion to the
/// caller's input domain is `x_break = strain · scale`. `slopes` are
/// multipliers relative to the elastic stiffness `k_e`; the implicit first
/// segment (before the first breakpoint) has multiplier `1.0`.
///
/// `length_rescale = l_calib / l_bond` is the crack-band rescaling factor.
/// Post-yield strain extents `(ε[i] − ε[0])` are multiplied by
/// `length_rescale`, **and post-yield slopes are divided by `length_rescale`**.
/// Both rescalings together preserve the envelope in plastic-displacement
/// space (`u = (ε − ε[0]) · l_bond`), so the per-bond plastic dissipation
/// `∫ F du` is invariant in `l_bond`. `ε[0]` (elastic yield) and the elastic
/// slope `k_e` are treated as material constants and are **not** rescaled.
/// Pass `1.0` for no rescaling — the default for callers without
/// `length_calibration`.
fn evaluate_piecewise(
    x_mag: f64,
    k_e: f64,
    breakpoint_strains: &[f64],
    slope_multipliers: &[f64],
    scale: f64,
    length_rescale: f64,
) -> f64 {
    let eps_yield = match breakpoint_strains.first() {
        Some(&e) => e,
        None => return k_e * x_mag,
    };
    let inv_rescale = if length_rescale > 0.0 { 1.0 / length_rescale } else { 1.0 };
    let mut f_acc = 0.0;
    let mut x_prev = 0.0;
    let mut slope = k_e;
    for (i, &eps_break_ref) in breakpoint_strains.iter().enumerate() {
        // First breakpoint (elastic yield) is a material constant; subsequent
        // breakpoints rescale by `length_rescale`. With `length_rescale = 1.0`
        // this collapses to the original strain-space envelope.
        let eps_break_eff = if i == 0 {
            eps_break_ref
        } else {
            eps_yield + (eps_break_ref - eps_yield) * length_rescale
        };
        let x_break = eps_break_eff * scale;
        if x_mag <= x_break {
            return f_acc + slope * (x_mag - x_prev);
        }
        f_acc += slope * (x_break - x_prev);
        x_prev = x_break;
        // Post-yield slope: rescale by 1/length_rescale so the stress at the
        // rescaled breakpoint matches the reference (length-invariant σ in
        // plastic-displacement space).
        slope = k_e * slope_multipliers[i] * inv_rescale;
    }
    f_acc + slope * (x_mag - x_prev)
}

/// Run-time plasticity state for one bond crate, built once at setup.
#[derive(Clone, Debug, Default)]
pub struct BondPlasticityModel {
    /// Active bending-plasticity model, or `None` if bending is purely elastic.
    pub bending: Option<BendingPlasticity>,
    /// Active axial-plasticity model, or `None` if axial is purely elastic.
    pub axial:   Option<AxialPlasticity>,
}

impl BondPlasticityModel {
    /// Build the run-time model from configuration. `youngs_modulus` is the
    /// material `E_b` from `[bonds]`; it's only required when bending
    /// plasticity uses [`BendingPlasticityConfig::GuoTrilinear`] (which
    /// derives breakpoints from `σ_0 / E_b`). Validates breakpoint ordering
    /// and length agreement; panics on invalid config.
    pub fn from_config(cfg: Option<&PlasticityConfig>, youngs_modulus: Option<f64>) -> Self {
        let mut out = Self::default();
        let cfg = match cfg { Some(c) => c, None => return out };
        out.bending = cfg.bending.as_ref().map(|b| match b {
            BendingPlasticityConfig::GuoBending { yield_stress } =>
                BendingPlasticity::GuoBending { yield_stress: *yield_stress },
            BendingPlasticityConfig::GuoTrilinear { yield_stress } => {
                let e_b = youngs_modulus.expect(
                    "[bonds] youngs_modulus must be set when using \
                     [bonds.plasticity.bending] kind = \"guo_trilinear\""
                );
                let eps_e = yield_stress / e_b;
                // Guo Eq. 33 in ε-space: ε_p / ε_e = (32 − 3π) / (3π).
                let pi = std::f64::consts::PI;
                let eps_p = eps_e * (32.0 - 3.0 * pi) / (3.0 * pi);
                BendingPlasticity::Piecewise {
                    breakpoint_strains: vec![eps_e, eps_p],
                    slope_multipliers: vec![0.5, 0.0],
                    length_calibration: None,
                }
            }
            BendingPlasticityConfig::Piecewise { breakpoint_strains, slope_multipliers, length_calibration } => {
                validate_piecewise("bending", breakpoint_strains, slope_multipliers);
                BendingPlasticity::Piecewise {
                    breakpoint_strains: breakpoint_strains.clone(),
                    slope_multipliers: slope_multipliers.clone(),
                    length_calibration: *length_calibration,
                }
            }
        });
        out.axial = cfg.axial.as_ref().map(|a| match a {
            AxialPlasticityConfig::Piecewise { breakpoint_strains, slope_multipliers, length_calibration } => {
                validate_piecewise("axial", breakpoint_strains, slope_multipliers);
                AxialPlasticity::Piecewise {
                    breakpoint_strains: breakpoint_strains.clone(),
                    slope_multipliers: slope_multipliers.clone(),
                    length_calibration: *length_calibration,
                }
            }
        });
        out
    }
}

fn validate_piecewise(channel: &str, breaks: &[f64], slopes: &[f64]) {
    assert_eq!(
        breaks.len(), slopes.len(),
        "[bonds.plasticity.{channel}] Piecewise: `breakpoint_strains` and \
         `slope_multipliers` must have equal length",
    );
    for w in breaks.windows(2) {
        assert!(
            w[1] > w[0],
            "[bonds.plasticity.{channel}] Piecewise: `breakpoint_strains` must be strictly ascending",
        );
    }
    assert!(
        breaks.is_empty() || breaks[0] >= 0.0,
        "[bonds.plasticity.{channel}] Piecewise: `breakpoint_strains` must be non-negative",
    );
}

// ── Per-bond return mappings ────────────────────────────────────────────────

/// Apply the bending-plasticity return-map.
///
/// `theta_bend` is the kinematic bending-angle vector (⊥ to bond axis).
/// `theta_p_bend` is the current plastic anchor. `theta_max_bend ≥ 0` is the
/// largest `|θ_bend|` ever reached on this bond (used to evaluate the
/// monotonic envelope). `k_e = E_b · I / l_b` is the elastic bending
/// stiffness; `model` selects the envelope; `r_b`, `l_b` are bond geometry.
///
/// Returns `(M_bend, θ_p_bend_new, θ_max_bend_new)`.
pub fn update_bending(
    theta_bend: [f64; 3],
    theta_p_bend: [f64; 3],
    theta_max_bend: f64,
    k_e: f64,
    model: &BendingPlasticity,
    r_b: f64,
    l_b: f64,
) -> ([f64; 3], [f64; 3], f64) {
    let theta_bend_mag2 =
        theta_bend[0] * theta_bend[0] + theta_bend[1] * theta_bend[1] + theta_bend[2] * theta_bend[2];
    let theta_bend_mag = theta_bend_mag2.sqrt();
    let theta_max_new = theta_max_bend.max(theta_bend_mag);

    let theta_e = [
        theta_bend[0] - theta_p_bend[0],
        theta_bend[1] - theta_p_bend[1],
        theta_bend[2] - theta_p_bend[2],
    ];
    let theta_e_mag2 = theta_e[0] * theta_e[0] + theta_e[1] * theta_e[1] + theta_e[2] * theta_e[2];
    if theta_e_mag2 <= f64::MIN_POSITIVE {
        return ([0.0; 3], theta_p_bend, theta_max_new);
    }
    let theta_e_mag = theta_e_mag2.sqrt();

    let m_trial_mag = k_e * theta_e_mag;
    // Monotonic envelope cap: evaluate at the larger of the historical
    // maximum and the current kinematic magnitude.
    let m_env_mag = model.envelope(theta_max_new, k_e, r_b, l_b);

    if m_trial_mag <= m_env_mag {
        let m = [k_e * theta_e[0], k_e * theta_e[1], k_e * theta_e[2]];
        (m, theta_p_bend, theta_max_new)
    } else {
        let scale = m_env_mag / m_trial_mag;
        let m = [k_e * theta_e[0] * scale, k_e * theta_e[1] * scale, k_e * theta_e[2] * scale];
        let elastic_offset = m_env_mag / k_e;
        let dir = [theta_e[0] / theta_e_mag, theta_e[1] / theta_e_mag, theta_e[2] / theta_e_mag];
        let theta_p_new = [
            theta_bend[0] - elastic_offset * dir[0],
            theta_bend[1] - elastic_offset * dir[1],
            theta_bend[2] - elastic_offset * dir[2],
        ];
        (m, theta_p_new, theta_max_new)
    }
}

/// Apply the axial-plasticity return-map.
///
/// `eps_axial = (L − L₀)/L₀` is the signed kinematic axial strain;
/// `eps_p_axial` is the signed plastic anchor; `eps_max_axial ≥ 0` is the
/// largest `|ε_axial|` ever reached on this bond. `k_n = E_b · A / L₀` is
/// the elastic axial stiffness (N/m); `l0 = L₀` is the equilibrium bond
/// length. The conservative axial force is `F = k_n · L₀ · ε_e = E_b·A · ε_e`.
///
/// Returns `(F_n, ε_p_new, ε_max_new)`.
pub fn update_axial(
    eps_axial: f64,
    eps_p_axial: f64,
    eps_max_axial: f64,
    k_n: f64,
    l0: f64,
    model: &AxialPlasticity,
) -> (f64, f64, f64) {
    let eps_max_new = eps_max_axial.max(eps_axial.abs());
    let eps_e = eps_axial - eps_p_axial;
    let eps_e_mag = eps_e.abs();
    if eps_e_mag <= f64::MIN_POSITIVE {
        return (0.0, eps_p_axial, eps_max_new);
    }
    let k_eff = k_n * l0;
    let f_trial_mag = k_eff * eps_e_mag;
    let f_env_mag = model.envelope(eps_max_new, k_eff, l0);

    if f_trial_mag <= f_env_mag {
        (k_eff * eps_e, eps_p_axial, eps_max_new)
    } else {
        let dir = eps_e.signum();
        let f_n = f_env_mag * dir;
        let elastic_offset = f_env_mag / k_eff;
        let eps_p_new = eps_axial - elastic_offset * dir;
        (f_n, eps_p_new, eps_max_new)
    }
}

/// Fully-plastic bending-moment cap `M^p = (4/3) σ_0 r_b³` (Guo Eq. 31).
#[inline]
pub fn fully_plastic_moment(yield_stress: f64, r_b: f64) -> f64 {
    (4.0 / 3.0) * yield_stress * r_b.powi(3)
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (f64, f64, f64, f64, f64, f64) {
        // r_b = 1 mm, l_b = 2 mm, E_b = 1 GPa, σ_0 = 10 MPa.
        let r_b: f64 = 1.0e-3;
        let l_b: f64 = 2.0e-3;
        let e_b: f64 = 1.0e9;
        let sigma_y: f64 = 1.0e7;
        let area = std::f64::consts::PI * r_b * r_b;
        let i_ben = std::f64::consts::PI * r_b.powi(4) / 4.0;
        let k_e = e_b * i_ben / l_b;
        let k_n = e_b * area / l_b;
        (r_b, l_b, sigma_y, k_e, k_n, e_b)
    }

    fn guo(sigma_y: f64) -> BendingPlasticity {
        BendingPlasticity::GuoBending { yield_stress: sigma_y }
    }

    // ── Envelope shape ──────────────────────────────────────────────────────

    #[test]
    fn fully_plastic_formula() {
        let r_b: f64 = 2.0e-3;
        let sigma_y: f64 = 5.0e7;
        let m_p = fully_plastic_moment(sigma_y, r_b);
        assert!((m_p - (4.0 / 3.0) * sigma_y * r_b.powi(3)).abs() / m_p < 1e-14);
    }

    #[test]
    fn guo_bending_envelope_elastic_and_capped() {
        let (r_b, l_b, sigma_y, k_e, _, _) = fixture();
        let model = guo(sigma_y);
        let m_p = fully_plastic_moment(sigma_y, r_b);
        let theta_yield = m_p / k_e;
        let env_under = model.envelope(0.5 * theta_yield, k_e, r_b, l_b);
        assert!((env_under - 0.5 * m_p).abs() / m_p < 1e-12);
        let env_over = model.envelope(5.0 * theta_yield, k_e, r_b, l_b);
        assert!((env_over - m_p).abs() / m_p < 1e-12);
    }

    #[test]
    fn bending_piecewise_three_breakpoints_user_example() {
        // "elastic K_e up to 1% strain, then K_e/2, then K_e/10, then flat".
        let (r_b, l_b, _, k_e, _, _) = fixture();
        let model = BendingPlasticity::Piecewise {
            breakpoint_strains: vec![0.01, 0.02, 0.03],
            slope_multipliers:  vec![0.5,  0.1,  0.0],
            length_calibration: None,
        };
        let scale = l_b / r_b;
        let theta_bp = [0.01 * scale, 0.02 * scale, 0.03 * scale];
        let m_at_1 = k_e * theta_bp[0];
        let m_at_2 = m_at_1 + 0.5 * k_e * (theta_bp[1] - theta_bp[0]);
        let m_at_3 = m_at_2 + 0.1 * k_e * (theta_bp[2] - theta_bp[1]);
        assert!((model.envelope(theta_bp[0], k_e, r_b, l_b) - m_at_1).abs() / m_at_1 < 1e-12);
        assert!((model.envelope(theta_bp[1], k_e, r_b, l_b) - m_at_2).abs() / m_at_2 < 1e-12);
        assert!((model.envelope(theta_bp[2], k_e, r_b, l_b) - m_at_3).abs() / m_at_3 < 1e-12);
        let env_far = model.envelope(10.0 * theta_bp[2], k_e, r_b, l_b);
        assert!((env_far - m_at_3).abs() / m_at_3 < 1e-12);
    }

    #[test]
    fn axial_piecewise_three_breakpoints_user_example() {
        // Same shape as the bending test but on the axial channel.
        let (_r, _l, _, _, k_n, e_b) = fixture();
        let l0 = 2.0e-3;
        let area = std::f64::consts::PI * (1.0e-3_f64).powi(2);
        let k_eff = k_n * l0;        // = E_b · A
        assert!((k_eff - e_b * area).abs() / k_eff < 1e-12);
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01, 0.02, 0.03],
            slope_multipliers:  vec![0.5,  0.1,  0.0],
            length_calibration: None,
        };
        let f_at_1 = k_eff * 0.01;
        let f_at_2 = f_at_1 + 0.5 * k_eff * 0.01;
        let f_at_3 = f_at_2 + 0.1 * k_eff * 0.01;
        assert!((model.envelope(0.01, k_eff, l0) - f_at_1).abs() / f_at_1 < 1e-12);
        assert!((model.envelope(0.02, k_eff, l0) - f_at_2).abs() / f_at_2 < 1e-12);
        assert!((model.envelope(0.03, k_eff, l0) - f_at_3).abs() / f_at_3 < 1e-12);
        let env_far = model.envelope(0.30, k_eff, l0);
        assert!((env_far - f_at_3).abs() / f_at_3 < 1e-12);
    }

    // ── Return-map: bending ─────────────────────────────────────────────────

    #[test]
    fn bending_elastic_below_yield() {
        let (r_b, l_b, sigma_y, k_e, _, _) = fixture();
        let model = guo(sigma_y);
        let m_p = fully_plastic_moment(sigma_y, r_b);
        let theta_yield = m_p / k_e;
        let theta_bend = [0.5 * theta_yield, 0.0, 0.0];
        let (m, theta_p_new, theta_max_new) =
            update_bending(theta_bend, [0.0; 3], 0.0, k_e, &model, r_b, l_b);
        assert!((m[0] - k_e * theta_bend[0]).abs() / (k_e * theta_bend[0]).abs() < 1e-12);
        assert_eq!(theta_p_new, [0.0; 3]);
        assert!((theta_max_new - 0.5 * theta_yield).abs() < 1e-15);
    }

    #[test]
    fn bending_moment_capped_at_m_p() {
        let (r_b, l_b, sigma_y, k_e, _, _) = fixture();
        let model = guo(sigma_y);
        let m_p = fully_plastic_moment(sigma_y, r_b);
        let theta_yield = m_p / k_e;
        let theta_bend = [3.0 * theta_yield, 0.0, 0.0];
        let (m, theta_p_new, theta_max_new) =
            update_bending(theta_bend, [0.0; 3], 0.0, k_e, &model, r_b, l_b);
        let m_mag = (m[0].powi(2) + m[1].powi(2) + m[2].powi(2)).sqrt();
        assert!((m_mag - m_p).abs() / m_p < 1e-12);
        assert!((theta_p_new[0] - (theta_bend[0] - theta_yield)).abs() < 1e-12);
        assert!((theta_max_new - 3.0 * theta_yield).abs() < 1e-15);
    }

    #[test]
    fn bending_unload_elastic_slope_then_reverse_kinematic_hardening() {
        let (r_b, l_b, sigma_y, k_e, _, _) = fixture();
        let model = guo(sigma_y);
        let m_p = fully_plastic_moment(sigma_y, r_b);
        let theta_yield = m_p / k_e;

        let theta_load = [3.0 * theta_yield, 0.0, 0.0];
        let (_, theta_p_load, theta_max_load) =
            update_bending(theta_load, [0.0; 3], 0.0, k_e, &model, r_b, l_b);
        let residual = theta_p_load[0];
        assert!(residual > 0.0);

        // Partial unload.
        let theta_unload = [theta_load[0] - 0.5 * theta_yield, 0.0, 0.0];
        let (m_unload, _, _) =
            update_bending(theta_unload, theta_p_load, theta_max_load, k_e, &model, r_b, l_b);
        let theta_e_after = theta_unload[0] - residual;
        assert!((m_unload[0] - k_e * theta_e_after).abs() / (k_e * theta_e_after).abs() < 1e-12);

        // Reverse past negative yield.
        let theta_reverse = [residual - 2.0 * theta_yield, 0.0, 0.0];
        let (m_rev, theta_p_rev, _) =
            update_bending(theta_reverse, theta_p_load, theta_max_load, k_e, &model, r_b, l_b);
        assert!((m_rev[0].abs() - m_p).abs() / m_p < 1e-12);
        assert!(m_rev[0] < 0.0);
        let expected_anchor = theta_reverse[0] + theta_yield;
        assert!((theta_p_rev[0] - expected_anchor).abs() / expected_anchor.abs() < 1e-12);
    }

    #[test]
    fn guo_trilinear_config_expands_to_piecewise_with_correct_breakpoints() {
        // Round-trip: configure GuoTrilinear → expand → check the resulting
        // Piecewise model's envelope hits the right values at the two
        // breakpoints and the cap.
        let sigma_y = 5.0e6;
        let e_b = 1.0e9;
        let cfg = PlasticityConfig {
            bending: Some(BendingPlasticityConfig::GuoTrilinear { yield_stress: sigma_y }),
            axial: None,
        };
        let model = BondPlasticityModel::from_config(Some(&cfg), Some(e_b));
        let bending = model.bending.expect("trilinear should produce a bending model");

        // It should expand to a Piecewise with the right breakpoints.
        match bending {
            BendingPlasticity::Piecewise { ref breakpoint_strains, ref slope_multipliers, .. } => {
                assert_eq!(slope_multipliers, &vec![0.5, 0.0],
                    "trilinear → slopes [K_e/2, 0]");
                let pi = std::f64::consts::PI;
                let eps_e_expected = sigma_y / e_b;
                let eps_p_expected = eps_e_expected * (32.0 - 3.0 * pi) / (3.0 * pi);
                assert!((breakpoint_strains[0] - eps_e_expected).abs() / eps_e_expected < 1e-12);
                assert!((breakpoint_strains[1] - eps_p_expected).abs() / eps_p_expected < 1e-12);
            }
            _ => panic!("expected Piecewise variant"),
        }

        // And the envelope at the breakpoints should hit M^e and M^p exactly.
        let (r_b, l_b, _sy, k_e, _, _) = fixture();
        let _ = _sy;
        // Override σ_y for this test:
        let sigma_y = 5.0e6;
        let m_p = fully_plastic_moment(sigma_y, r_b);
        let i_ben = std::f64::consts::PI * r_b.powi(4) / 4.0;
        let m_e = sigma_y * i_ben / r_b;
        let theta_e = sigma_y * l_b / (e_b * r_b);
        let theta_p = theta_e * (32.0 - 3.0 * std::f64::consts::PI) / (3.0 * std::f64::consts::PI);
        // Rebuild model with the right E_b.
        let model = BondPlasticityModel::from_config(Some(&cfg), Some(e_b));
        let bending = model.bending.unwrap();
        assert!((bending.envelope(theta_e, k_e, r_b, l_b) - m_e).abs() / m_e < 1e-10);
        assert!((bending.envelope(theta_p, k_e, r_b, l_b) - m_p).abs() / m_p < 1e-10);
        assert!((bending.envelope(3.0 * theta_p, k_e, r_b, l_b) - m_p).abs() / m_p < 1e-10);
    }

    #[test]
    fn bending_piecewise_recovers_guo_trilinear_shape() {
        let (r_b, l_b, sigma_y, k_e, _, e_b) = fixture();
        let eps_e = sigma_y / e_b;
        let m_p = fully_plastic_moment(sigma_y, r_b);
        let i_ben = std::f64::consts::PI * r_b.powi(4) / 4.0;
        let m_e = sigma_y * i_ben / r_b;
        let k_ep = k_e / 2.0;
        let theta_e = sigma_y * l_b / (e_b * r_b);
        let theta_p = theta_e + (m_p - m_e) / k_ep;
        let eps_p = theta_p * r_b / l_b;

        let model = BendingPlasticity::Piecewise {
            breakpoint_strains: vec![eps_e, eps_p],
            slope_multipliers:  vec![0.5,   0.0],
            length_calibration: None,
        };
        assert!((model.envelope(theta_e, k_e, r_b, l_b) - m_e).abs() / m_e < 1e-10);
        assert!((model.envelope(theta_p, k_e, r_b, l_b) - m_p).abs() / m_p < 1e-10);
        assert!((model.envelope(2.0 * theta_p, k_e, r_b, l_b) - m_p).abs() / m_p < 1e-10);
    }

    // ── Return-map: axial ───────────────────────────────────────────────────

    #[test]
    fn axial_elastic_below_yield() {
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01],
            slope_multipliers:  vec![0.0],
            length_calibration: None,
        };
        let eps_axial = 0.005;
        let (f, eps_p_new, eps_max_new) = update_axial(eps_axial, 0.0, 0.0, k_n, l0, &model);
        let expected_f = k_n * l0 * eps_axial;
        assert!((f - expected_f).abs() / expected_f.abs() < 1e-12);
        assert_eq!(eps_p_new, 0.0);
        assert!((eps_max_new - 0.005).abs() < 1e-15);
    }

    #[test]
    fn axial_force_capped_at_envelope_when_yielded() {
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01],
            slope_multipliers:  vec![0.0],
            length_calibration: None,
        };
        let k_eff = k_n * l0;
        let f_yield = k_eff * 0.01;
        let (f, eps_p_new, eps_max_new) = update_axial(0.03, 0.0, 0.0, k_n, l0, &model);
        assert!((f - f_yield).abs() / f_yield < 1e-12);
        assert!((eps_p_new - 0.02).abs() / 0.02 < 1e-12);
        assert!((eps_max_new - 0.03).abs() / 0.03 < 1e-12);
    }

    #[test]
    fn axial_compression_capped_symmetrically() {
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01],
            slope_multipliers:  vec![0.0],
            length_calibration: None,
        };
        let k_eff = k_n * l0;
        let f_yield = k_eff * 0.01;
        let (f, eps_p_new, eps_max_new) = update_axial(-0.03, 0.0, 0.0, k_n, l0, &model);
        assert!((f - (-f_yield)).abs() / f_yield < 1e-12);
        assert!((eps_p_new - (-0.02)).abs() / 0.02 < 1e-12);
        assert!((eps_max_new - 0.03).abs() / 0.03 < 1e-12, "max tracks |ε|");
    }

    #[test]
    fn axial_unload_elastic_then_reverse_kinematic_hardening() {
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01],
            slope_multipliers:  vec![0.0],
            length_calibration: None,
        };
        let k_eff = k_n * l0;
        let f_yield = k_eff * 0.01;
        let (_, eps_p_load, eps_max_load) = update_axial(0.03, 0.0, 0.0, k_n, l0, &model);
        // Partial unload: ε = 0.025 → ε_e = 0.005 (elastic), F = 0.5·f_yield.
        let (f_unload, eps_p_after, _) =
            update_axial(0.025, eps_p_load, eps_max_load, k_n, l0, &model);
        assert!((f_unload - 0.5 * f_yield).abs() / f_yield < 1e-12);
        assert!((eps_p_after - eps_p_load).abs() < 1e-15);
        // Continue all the way to ε = 0 → ε_e = −0.02 → F = −f_yield.
        let (f_zero, eps_p_zero, _) =
            update_axial(0.0, eps_p_load, eps_max_load, k_n, l0, &model);
        assert!((f_zero - (-f_yield)).abs() / f_yield < 1e-12);
        assert!((eps_p_zero - 0.01).abs() / 0.01 < 1e-12);
    }

    #[test]
    fn axial_multi_segment_traces_envelope() {
        // Three-segment user example, monotonic axial loading.
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01, 0.02, 0.03],
            slope_multipliers:  vec![0.5,  0.1,  0.0],
            length_calibration: None,
        };
        let k_eff = k_n * l0;
        let f_at = |eps: f64| {
            let f1 = k_eff * 0.01;
            let f2 = f1 + 0.5 * k_eff * 0.01;
            let f3 = f2 + 0.1 * k_eff * 0.01;
            if eps <= 0.01 { k_eff * eps }
            else if eps <= 0.02 { f1 + 0.5 * k_eff * (eps - 0.01) }
            else if eps <= 0.03 { f2 + 0.1 * k_eff * (eps - 0.02) }
            else { f3 }
        };
        let mut eps_p = 0.0;
        let mut eps_max = 0.0;
        for eps in [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.10_f64] {
            let (f, eps_p_new, eps_max_new) =
                update_axial(eps, eps_p, eps_max, k_n, l0, &model);
            let expected = f_at(eps);
            assert!(
                (f - expected).abs() / expected.max(1e-30) < 1e-10,
                "axial envelope mismatch at ε = {}: got {} vs expected {}",
                eps, f, expected,
            );
            eps_p = eps_p_new;
            eps_max = eps_max_new;
        }
    }

    #[test]
    fn axial_cycling_freezes_cap_at_eps_max() {
        // After loading to ε_max = 0.02 (F_cap = 0.015·k_eff), unloading
        // through zero and pushing into the negative side should yield at
        // F = −0.015·k_eff (not at the larger F_env(|ε|) the monotonic curve
        // would say in segment past first yield).
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.01, 0.02, 0.03],
            slope_multipliers:  vec![0.5,  0.1,  0.0],
            length_calibration: None,
        };
        let k_eff = k_n * l0;
        // Monotonic to 0.02 — cap should be 0.015·k_eff.
        let (_, eps_p_load, eps_max_load) =
            update_axial(0.02, 0.0, 0.0, k_n, l0, &model);
        let f_cap = 0.015 * k_eff;
        // Go all the way to ε = -0.02 (symmetric extremum). Trial elastic
        // force is way past the cap, so we yield on the negative side at the
        // same |F| = f_cap (kinematic hardening — the cap doesn't grow until
        // |ε| > 0.02).
        let (f_rev, _eps_p_rev, eps_max_rev) =
            update_axial(-0.02, eps_p_load, eps_max_load, k_n, l0, &model);
        assert!((f_rev.abs() - f_cap).abs() / f_cap < 1e-12);
        assert!(f_rev < 0.0);
        // Now push to ε = -0.025 — past prior |ε_max| of 0.02. Cap grows.
        let (f_grow, _, eps_max_grow) =
            update_axial(-0.025, _eps_p_rev, eps_max_rev, k_n, l0, &model);
        let expected_cap = 0.015 * k_eff + 0.1 * k_eff * 0.005;   // segment 3 slope
        assert!((f_grow.abs() - expected_cap).abs() / expected_cap < 1e-12);
        assert!((eps_max_grow - 0.025).abs() / 0.025 < 1e-12);
    }

    // ── Crack-band length rescaling ─────────────────────────────────────────

    #[test]
    fn axial_piecewise_length_calibration_none_matches_unregularized() {
        // With `length_calibration = None` the envelope must match the
        // unregularized one to bit precision — no rescaling at all.
        let (_r, _l, _, _, k_n, _) = fixture();
        let l0 = 2.0e-3;
        let model_off = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.02, 0.04],
            slope_multipliers:  vec![0.3,  0.0],
            length_calibration: None,
        };
        let model_on_at_calib = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.02, 0.04],
            slope_multipliers:  vec![0.3,  0.0],
            length_calibration: Some(l0),
        };
        let k_eff = k_n * l0;
        // At l_bond = l_calib, length_rescale = 1 ⇒ identical envelopes.
        for eps in [0.005_f64, 0.02, 0.03, 0.04, 0.10] {
            let v_off = model_off.envelope(eps, k_eff, l0);
            let v_on  = model_on_at_calib.envelope(eps, k_eff, l0);
            assert!(
                (v_off - v_on).abs() / v_off.max(1e-30) < 1e-15,
                "rescaling at l_bond = l_calib must be a no-op (eps = {eps}: off = {v_off}, on = {v_on})",
            );
        }
    }

    #[test]
    fn axial_piecewise_post_yield_strain_extent_scales_inversely_with_bond_length() {
        // Calibrated at l_ref = 2 mm with breakpoints [0.02, 0.04] (post-yield
        // strain extent = 0.02). At l_bond = l_ref/2 = 1 mm the extent should
        // double (rescale factor α = 2), so the second breakpoint is at
        // ε = 0.02 + 0.02·2 = 0.06; at l_bond = 2·l_ref = 4 mm the extent halves
        // (α = 0.5) so the second breakpoint is at ε = 0.02 + 0.02·0.5 = 0.03.
        // Yield (ε = 0.02) is preserved at every bond length.
        let (_r, _l, _, _, k_n, _) = fixture();
        let l_ref = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.02, 0.04],
            slope_multipliers:  vec![0.3,  0.0],
            length_calibration: Some(l_ref),
        };

        // Helper: at given l_bond, locate the strain where the envelope hits its
        // plastic plateau (∂F/∂ε = 0 → segment past the second breakpoint).
        let f_at = |eps: f64, l_bond: f64| {
            let k_eff = k_n * l_bond;
            model.envelope(eps, k_eff, l_bond)
        };

        // Elastic yield breakpoint preserved at every bond length.
        for &l_bond in &[0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3] {
            let k_eff = k_n * l_bond;
            // Force at yield ε = 0.02 should be k_eff · 0.02 regardless of l_bond.
            let f_y = f_at(0.02, l_bond);
            assert!(
                (f_y - 0.02 * k_eff).abs() / (0.02 * k_eff) < 1e-12,
                "yield breakpoint must not rescale (l_bond = {l_bond:.3e}: f_y = {f_y})",
            );
        }

        // At l_bond = l_ref/2 (α = 2), the plateau starts at ε = 0.06 (not 0.04)
        // AND the hardening slope is 0.3/α = 0.15·k_eff (post-yield slopes
        // rescale as 1/α to keep σ invariant in plastic-displacement space).
        let l_short = l_ref / 2.0;
        let k_eff_short = k_n * l_short;
        let alpha_short = l_ref / l_short;       // = 2
        let slope_h_short = 0.3 / alpha_short * k_eff_short;
        let f_y         = f_at(0.02,  l_short);
        let f_eps_0_055 = f_at(0.055, l_short);
        let expected    = f_y + slope_h_short * (0.055 - 0.02);
        assert!((f_eps_0_055 - expected).abs() / expected < 1e-12,
            "rescaled hardening slope at l_ref/2: f(0.055) = {f_eps_0_055}, expected {expected}");
        // Just above 0.06: on the plateau (slope = 0), value frozen at f(0.06).
        let f_eps_0_07 = f_at(0.07, l_short);
        let f_eps_0_06 = f_at(0.06, l_short);
        assert!((f_eps_0_07 - f_eps_0_06).abs() / f_eps_0_06 < 1e-12);

        // Crack-band invariance: F at the same plastic-displacement u must be
        // the same regardless of l_bond. This is the cleanest invariant to
        // assert — and it directly implies fracture-energy invariance ∫F du.
        // k_eff = E·A is constant across bond lengths (K_n = E·A/L₀ scales
        // inversely with L₀, so K_n · L₀ = E·A is fixed).
        let k_eff = k_n * 2.0e-3;     // E·A, the material-invariant constant
        for &u in &[1.0e-5_f64, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5, 6.0e-5, 8.0e-5] {
            let f_for = |l_bond: f64| {
                let eps = u / l_bond + 0.02;
                model.envelope(eps, k_eff, l_bond)
            };
            let f_short = f_for(l_ref / 2.0);
            let f_calib = f_for(l_ref);
            let f_long  = f_for(2.0 * l_ref);
            assert!((f_short - f_calib).abs() / f_calib.abs() < 1e-12,
                "F(u={u:.2e}) should be l_bond-invariant: l_ref/2 → {f_short}, l_ref → {f_calib}");
            assert!((f_long - f_calib).abs() / f_calib.abs() < 1e-12,
                "F(u={u:.2e}) should be l_bond-invariant: 2·l_ref → {f_long}, l_ref → {f_calib}");
        }
    }

    #[test]
    fn axial_piecewise_fracture_energy_per_unit_fiber_length_is_invariant() {
        // Crack-band invariance: ∫ F(u) du across a fixed plastic-displacement
        // range u_max is independent of l_bond when length_calibration is set.
        // This is the actual fracture-energy invariant (Hillerborg).
        let (_r, _l, _, _, k_n, _) = fixture();
        let l_ref = 2.0e-3;
        let model = AxialPlasticity::Piecewise {
            breakpoint_strains: vec![0.02, 0.04],
            slope_multipliers:  vec![0.3,  0.0],
            length_calibration: Some(l_ref),
        };
        let k_eff = k_n * 2.0e-3;     // E·A, length-invariant
        let energy = |l_bond: f64| {
            let u_max = 0.02 * l_ref * 2.0;       // = 8e-5 m
            let n = 200_000;
            let mut e = 0.0;
            let du = u_max / n as f64;
            for i in 0..n {
                let u = (i as f64 + 0.5) * du;
                let eps = u / l_bond + 0.02;
                let f = model.envelope(eps, k_eff, l_bond);
                e += f * du;
            }
            e
        };
        let e_short = energy(l_ref / 2.0);
        let e_calib = energy(l_ref);
        let e_long  = energy(2.0 * l_ref);
        let rel = |a: f64, b: f64| (a - b).abs() / a.max(b);
        assert!(rel(e_short, e_calib) < 1e-4,
            "∫F du: l_ref/2 vs l_ref: {e_short} vs {e_calib}, rel diff {}",
            rel(e_short, e_calib));
        assert!(rel(e_long, e_calib) < 1e-4,
            "∫F du: 2·l_ref vs l_ref: {e_long} vs {e_calib}, rel diff {}",
            rel(e_long, e_calib));
    }
}
