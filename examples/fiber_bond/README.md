# fiber_bond — bonded-particle fiber validation harness

A single binary, multiple config files. Each config exercises one
deformation mode of the BPM bond model and writes a CSV that the
companion `validate.py` checks against the analytical predictions in
Guo et al. 2018 (*Chem. Eng. Sci.* **175**, 118–129).

## Quick start

```bash
# Build
cargo build --release --example fiber_bond --no-default-features

# Run + validate a scenario
cargo run --release --example fiber_bond --no-default-features -- \
    examples/fiber_bond/axial_elastic.toml
python3 examples/fiber_bond/validate.py \
    examples/fiber_bond/axial_elastic/data/fiber_bond.csv

cargo run --release --example fiber_bond --no-default-features -- \
    examples/fiber_bond/cantilever_bending.toml
python3 examples/fiber_bond/validate.py \
    examples/fiber_bond/cantilever_bending/data/fiber_bond.csv
```

## Scenarios

| Config | Mode | Closed-form prediction | Status |
|---|---|---|---|
| `axial_elastic.toml`           | Axial tension       | `σ = E·ε` (Guo Eq. 1)                       | PASS, 0.005% error |
| `cantilever_bending.toml`      | Cantilever bending  | `y(x) = F·x²·(3L−x)/(6·E_b·I)` (Guo Sec. 2.1) | PASS, ~0.6% error |
| `bending_vibration.toml`       | Free bending vibration | `T = 1.787·L²·√(ρ_l/EI)` (Guo Eq. 18, discrete-mass form) | PASS, ~4.7% error |
| `axial_plastic_piecewise.toml` | Axial plastic loading | Piecewise-linear hardening envelope (this code's config) | PASS, < 0.1% error |
| `bending_plastic_guo.toml`     | Guo three-step bending plasticity (trilinear) | `|M_bend| ≤ M_p = (4/3)·σ_0·r_b³` (Guo Eq. 31); envelope follows Eqs. 27, 29, 32, 33, 35 | PASS, M traces all three trilinear regimes |

## Breakage scenarios — see [`../fiber_bond_breakage`](../fiber_bond_breakage)

All breakage scenarios (criterion variants, Weibull thresholds, etc.)
live in the sibling `fiber_bond_breakage/` directory. They share this
example's binary (`fiber_bond`) — just point it at the config there.

Both scenarios use a 1 GPa, ν = 0.25, ρ = 2500 kg/m³ "BPM" material with
critical per-channel bond damping (β = 1.0) plus a global viscous fix to
damp the lowest bending mode.

## Particle geometry — touching vs. spaced

The two scenarios use **different fiber CSVs** for a reason:

* `axial_elastic.toml` → `fiber_11.csv` — 11 spheres of radius 1 mm at
  2 mm centre-to-centre spacing (spheres touching). Bond AR = 1.
* `cantilever_bending.toml` → `fiber_11_spaced.csv` — 11 spheres of
  radius 1 mm at **4 mm centre-to-centre spacing** (2 mm gap between
  sphere surfaces). Bond AR = 2.

Touching-sphere chains have **bond aspect ratio** `L_bond / (2·r_bond) = 1`.
At AR = 1 the bonded-particle model carries transverse load through a
coupled shear / bending pathway that adds extra stiffness vs. the
Euler-Bernoulli continuum limit. Empirically, a touching-sphere chain at
these parameters lands ~5× stiffer than EB beam theory.

| Variant | L_bond | r_bond | AR | measured / EB |
|---|---|---|---|---|
| Touching, r_bond = 1 mm   | 2 mm | 1 mm  | 1.0 | 0.22 |
| Thin bond, r_bond = 0.5 mm | 2 mm | 0.5 mm | 1.0 *(length unchanged)* | 0.21 |
| Spaced, L_bond = 4 mm     | 4 mm | 1 mm  | 2.0 | **0.99** |

Conclusion: **bond LENGTH (not bond radius)** is what controls
EB-compliance. Shrinking the bond cross-section scales all four channel
stiffnesses together and leaves the relative bending/shear coupling
unchanged. Spacing the spheres genuinely lengthens the bond cylinder,
weakens the shear lever-arm coupling, and recovers EB to ~1 %.

This is the bonded-sphere analog of the design choice in
Guo 2018, which uses **sphero-cylinder** elements (capsule shape) so
that the bond cylinder is naturally slender, independent of the
hemispherical end-cap radius.

Mass is still lumped at sphere centres (massless bonds), so a spaced
chain has a different mass distribution than a continuum beam of the
same length. For static deflection this is irrelevant; for vibration
tests it will shift natural frequencies by the usual lumped-mass
correction.

## What gets recorded

Per sampled step (`record_every = 200` steps), the CSV captures:

* Endpoint positions and velocities (leftmost / rightmost atom by
  initial x — picked at setup).
* Middle-bond kinematics: rest length, current length, axial strain,
  shear displacement magnitude (|Δs|), bending angle magnitude
  (|Δθ_bend|), twist angle (Δθ · n̂).
* Middle-bond plastic state from `BondHistoryEntry`: `θ_p_bend`,
  `ε_p_axial`, `θ_max_bend`, `ε_max_axial` (all zero in elastic-only
  runs).
* Bond count and cumulative `bonds_broken`.

The recorder intentionally does **not** capture force/moment directly —
they are reconstructed from the bond stiffness and the recorded
kinematics inside `validate.py`. That keeps the recorder cheap and
keeps it working unchanged for future runs with axial/bending plastic
configs.

## A note on the vibration prediction

The textbook Guo Eq. 18 form `T_bend = 1.787·L²·√(ρ_l/EI)` uses
`ρ_l = ρ·A` — the continuum mass per unit length of a fully-filled beam
cross-section. Our spaced-sphere chain has **massless bonds and
discrete mass** lumped at the sphere centres, so its actual `ρ_l` is
`M_chain / L`, where `M_chain = Σ atom masses` (≈ 37 % of the continuum
value at our 4 mm spacing). The validator uses the discrete form
(`ρ_l_discrete = M_chain / L`) for the PASS/FAIL prediction and also
prints the continuum value as a reference. Both predictions are shown
in the per-run summary.

## A note on `bending_plastic_guo`

The scenario uses the **trilinear** Guo 2018 envelope (Eq. 32):

* elastic up to `M^e = σ_0·I/r_b`  (slope `K_e = E_b·I/l_b`)
* elasto-plastic with slope `K_ep = K_e/2` up to `M^p = (4/3)·σ_0·r_b³`
* perfectly plastic at `M^p`

This is configured directly via `kind = "guo_trilinear"`; at config-build
time `BondPlasticityModel::from_config` derives the two breakpoints in
extreme-fibre-strain space (`ε_e = σ_0/E_b`,
`ε_p = ε_e · (32 − 3π) / (3π) ≈ 2.395·σ_0/E_b`) and stores the result as
the generic `Piecewise` runtime variant. The same Piecewise return-map
runs at every step.

The `M(θ_bend)` hysteresis plot shows the DEM trajectory tracing the
full trilinear envelope on the first monotonic loading: elastic slope,
the elasto-plastic K_ep segment, then the perfectly-plastic plateau at
M_p. On unloading the DEM returns along the elastic slope K_e —
kinematic-hardening behaviour matching Guo Fig. 9 Path V.

Under repeated loading at the *same* peak force, plastic flow at any
given bond happens only on cycle 1 — kinematic hardening leaves the
anchor `θ_p` positioned so subsequent cycles reach the cap exactly with
slope K_e and never push past. The tip *does* continue to accumulate
permanent deflection across all three cycles (Fig.
`bending_plastic_timeline.png`) — a *geometric* softening effect from
large-deformation reconfiguring the chain between cycles, not new
plastic flow at the mid-bond.

Time-gated load schedule lives in `main.rs::three_step_force_at` and is
activated when the output directory contains `bending_plastic_guo`;
every other config gets a no-op.

## Carryover

The original carryover list is fully cleared:

* ~~MPI-stable Weibull threshold sampling~~ — landed via
  `per_bond_uniform_samples` (SplitMix64 → SmallRng → Weibull). Same
  bond → same threshold under any MPI decomposition; six unit tests
  pin the contract.
* ~~`fiber_bond_breakage` sibling example~~ — folded into this directory
  as the three `axial_*_constant` and `axial_stress_weibull` scenarios;
  reuses the same Rust binary and the same per-scenario validator
  dispatch in `validate.py`.

Possible Phase-3 follow-ups (not on the active list):

* Multi-criterion breakage scenarios (`CombinedStrain`, `InteractionLinear*`).
* Statistical Weibull-CDF validation across many seeded runs (Phase-2
  parameter sweep over `[bonds].seed`).
* Coupled bending-plastic + breakage on the same fiber.
