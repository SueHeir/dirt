# fiber_bond_breakage — breakage-criterion validation harness

Companion to [`../fiber_bond`](../fiber_bond). Same Rust binary (`fiber_bond`),
same recorder + per-scenario validator dispatch, but every config here
exercises a different breakage variant from the `[bonds.breakage]` menu
in `breakage::BreakageConfig`.

## Quick start

```bash
# Build (one-time)
cargo build --release --example fiber_bond --no-default-features

# Run + validate any scenario
cargo run --release --example fiber_bond --no-default-features -- \
    examples/fiber_bond_breakage/<config>.toml
python3 examples/fiber_bond_breakage/validate.py \
    examples/fiber_bond_breakage/<scenario>/data/fiber_bond.csv
```

## Scenarios

| Config | Criterion variant | Setup | Prediction | Status |
|---|---|---|---|---|
| `axial_stress_constant.toml`       | `AxialStress + Constant`              | touching fiber, axial pull          | `ε = σ_max / E`                    | PASS, **0.0%** |
| `axial_strain_constant.toml`       | `AxialStrain + Constant`              | touching fiber, axial pull          | `ε = ε_max`                        | PASS, **0.0%** |
| `axial_stress_weibull.toml`        | `AxialStress + Weibull` (per-bond)    | touching fiber, axial pull          | weakest-bond Weibull sample        | PASS, 0.4% |
| `combined_stress.toml`             | `CombinedStress` (Guo Eq. 16, P-C)    | spaced fiber, kinematic tip-bend    | `y_tip = σ_max·L_c²/(3·E·r_b)`     | PASS, 30%  |
| `combined_strain.toml`             | `CombinedStrain` (migration Eq. 1.7-1) | spaced fiber, kinematic tip-bend   | `y_tip = ε_max·L_c²/(3·r_b)`       | PASS, 30%  |
| `interaction_linear_stress.toml`   | `InteractionLinearStress` (Clemmer)   | spaced fiber, kinematic tip-bend, bending-only channel | same as CombinedStress              | PASS, 16%  |

## Notes on the prediction quality

The **axial** scenarios are exact to machine precision because uniform
axial pull produces uniform strain along the fiber — every bond sees the
same `ε`, so `σ_anchor = σ_mid = σ_global` and the first bond to break
hits its threshold at exactly the predicted strain.

The **cantilever-bend** scenarios (`combined_*` and
`interaction_linear_stress`) match the small-deformation Euler-Bernoulli
analytical prediction `y_break = σ_max·L_c²/(3·E·r_b)` to within ~30 %.
The discrepancy is structural, not a bug:

* The bonded-sphere chain has slightly higher relative bending angle at
  the pinned end than the continuum cantilever (a discrete-chain
  artefact at the [[pin]] boundary), so the anchor bond reaches its
  yield threshold earlier than EB beam theory predicts.
* At `v_z = -0.5 m/s` the pull is comparable to the natural bending
  period (T_bend ≈ 5.5 ms), so the chain is *not* perfectly
  quasi-static — there's some dynamic stress overshoot.

These are well-known discrete-DEM-vs-continuum-EB mismatches. The
criterion code paths are working correctly — they fire at the right
*qualitative* moment (when the anchor bond reaches its yield), just not
at the exact EB-predicted tip displacement. Tolerance is set to 35 %
for these scenarios.

## InteractionLinearStress — why bending-only?

The InteractionLinear envelope is

```
B = max{0, |F_n|/(A·σ_a,c) + |F_t|/(A·σ_s,c) + r_b·|M_b|/(I·σ_b,c) + r_b·|M_t|/(I_p·σ_t,c)}
break when B ≥ 1
```

A natural multi-channel cantilever-bend test would activate both the
shear and bending channels (axial and twist are ~zero in transverse
pull). But the shear channel has a **dynamic damping trap**: at the
moving tip's bond, the relative tangential velocity `v_t ≈ v_pull`
combined with critical bond-shear damping `γ_t = 2·√(m*·K_t) ≈ 2.6
N·s/m` produces a transient `F_t = γ_t·v_t ≈ 1.3 N → σ_shear ≈ 400 kPa`
at the very first integration step. That trips any moderate `σ_shear,c`
threshold at step 0 — before the chain has bent at all.

To dodge the trap without losing the test, the canonical config sets
only the bending channel's threshold; the other three are left as
`None`, so the envelope reduces to `r_b·|M_b|/(I·σ_b,c) ≥ 1` — same
prediction as `combined_stress.toml`, but routed through the
`InteractionLinearStress` code path (multi-channel sum + threshold dump
+ `Option<ThresholdDistribution>` handling all exercised).

A real multi-channel test that genuinely loads bending **and** shear
quasi-statically would need either a much slower pull (so dynamic damping
is negligible) or much lower per-channel bond damping (`β_shear ≪ 1`).
Flagged for a Phase-3 follow-up.

## Per-bond Weibull thresholds (MPI-stable)

The `axial_stress_weibull` scenario exercises the `per_bond_uniform_samples`
hash-based sampler: each bond's threshold is a pure function of its
canonical tag pair and `[bonds].seed`, so the same fiber decomposed
across any number of MPI ranks would produce **bit-identical** per-bond
thresholds. The recorder dumps every bond's sampled thresholds to
`bond_thresholds.csv` at setup, and the validator reads them directly
to predict the weakest bond's break — no need to re-implement the
SplitMix64 / SmallRng / Weibull sampler in Python.

## Carryover

* Multi-channel interaction with genuine bend + shear coupling — needs
  a quasi-static pull (`v_pull · T_bend / L_c ≪ 1`) or reduced bond-
  channel damping.
* Statistical Weibull-CDF validation across many seeded runs.
* Coupled plastic + breakage on the same fiber.
