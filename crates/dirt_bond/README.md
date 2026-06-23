# dirt_bond

Bonded Particle Model (BPM) forces for DIRT — elastic bonds, breakage, and plasticity.

## What it does

`dirt_bond` adds the `DemBondPlugin`, which treats each bond between two
particles as a **solid cylindrical beam** resisting four independent
deformation channels: axial stretch/compression, transverse shear, twist
(torsion about the bond axis), and bending (rotation perpendicular to the
axis). Stiffness can be derived from material properties (Young's / shear
modulus + bond radius) or supplied directly; damping is set as a
critical-damping ratio; bonds can break under a configurable failure
criterion and yield under configurable plasticity models.

Bonds are created at setup either by **auto-bonding** initially touching
pairs or by **loading a LAMMPS `Bonds` section**. Periodic minimum-image
distances are used so closed-loop / wrapping bonds get the correct rest
length. Across MPI ranks, `ghost_cutoff` is automatically extended to keep
bonded partners (and their 1-3 shared neighbours) visible as ghosts.
Granular contact between bonded pairs is suppressed via `soil_core`'s
`BondStore`, so the bond is the sole interaction until it breaks.

## Bond geometry

Per bond the cross-section is a solid cylinder of radius
`r_b = bond_radius_ratio · min(R_i, R_j)`, giving area `A = π r_b²`, polar
second moment `J = ½ π r_b⁴` (torsion), bending second moment
`I = ¼ π r_b⁴ = ½ J`, and length `L = r₀`.

## Four-channel force model

| Channel         | Stiffness (material mode) |
|-----------------|----------------------------|
| Normal          | `K_n   = E · A / L`        |
| Shear           | `K_t   = G · A / L`        |
| Twist (torsion) | `K_tor = G · J / L`        |
| Bending         | `K_bend = E · I / L`       |

The shear force acts at the bond mid-point, producing a lever-arm torque
`τ = (L/2) n̂ × F_t` on both particles. Per-channel damping uses a
critical-damping ratio `β ∈ [0, 1]` (`γ = 2 β √(m*·K)` for forces,
`2 β √(I*·K)` for moments) with reduced mass / reduced MOI of the pair;
each channel accepts an optional raw-`γ` override.

## Breakage

A bond breaks when the active criterion (`breakage::BreakageCriterion`)
trips. The menu (`breakage::BreakageConfig`) covers three families ×
force/stress/strain measures plus an `Unbreakable` default:

- `axial_force` / `axial_stress` / `axial_strain` — independent tensile +
  shear branches
- `combined_stress` / `combined_strain` — extreme-fibre beam criterion
  (Potyondy–Cundall / Guo), tensile and shear branches
- `interaction_linear_force` / `_stress` / `_strain` — single linear
  damage envelope `Σ |X_i| / X_{i,c} ≥ 1` over the four channels (LAMMPS
  `bpm/rotational`)

Each threshold is a `ThresholdDistribution` (`kind = …`), one of three:

- `constant` — the same `value` for every bond.
- `weibull` — length-scaled 2-parameter Weibull (`mean`, `m`, `l_calib`,
  optional `l_min`) giving the weakest-link size effect.
- `crack_band` — deterministic, length-rescaled threshold (Bažant 1976 /
  Hillerborg–Modéer–Petersson 1976): the part of the threshold above
  `eps_yield` scales as `l_ref / max(l_bond, l_min)` so the per-bond plastic +
  brittle energy budget × bond length stays invariant under mesh refinement.
  Fields: `value_ref`, `l_ref`, optional `eps_yield` (use `0.0` for
  force/stress criteria, `> 0` for strain criteria), optional `l_min`. This is
  the correct regularization for the strain criteria (`axial_strain`,
  `combined_strain`, `interaction_linear_strain`).

Per-bond thresholds are drawn once at creation from the bond's tag pair plus
`seed`, so the breakage pattern is **MPI-decomposition independent and
bit-reproducible**.

## Plasticity

The bending and axial channels are independently optional
(`plasticity::PlasticityConfig`); omit a channel to keep it purely elastic.

- bending: `guo_bending` (elastic–perfectly-plastic cap at
  `M^p = (4/3) σ_0 r_b³`), `guo_trilinear` (Guo 2018 Eq. 32: elastic →
  `K_ep = K_e/2` elasto-plastic → perfectly-plastic trilinear envelope;
  expands at setup to a two-breakpoint `piecewise` and needs
  `[bonds].youngs_modulus`), or `piecewise`
- axial: `piecewise` (breakpoints in axial strain, slope multipliers)

The `piecewise` variants (and `guo_trilinear`, via its expansion) accept an
optional `length_calibration` (m): when set, the post-yield strain breakpoints
rescale at runtime by `length_calibration / l_bond` (the elastic yield
breakpoint is preserved), the crack-band regularization that keeps per-bond
plastic dissipation × bond length invariant under mesh refinement. Omit it to
recover the unregularized envelope.

Plastic anchors and max-strain history are carried per bond in
`BondHistoryStore` and survive MPI communication and atom reordering.

## Key types

| Item | Role |
|------|------|
| `DemBondPlugin` | registers resources, setup (auto-bond / file load / ghost extend / init), and per-step force system |
| `BondConfig` | deserialized `[bonds]` TOML section |
| `BondHistoryStore` | per-atom `BondHistoryEntry` list (shear Δs, rotation Δθ, sampled thresholds, plastic state); implements `AtomData` |
| `BondBreakage` / `BondPlasticity` | active criterion / plasticity model, built at setup |
| `BondMetrics` | per-step strain average + cumulative bonds-broken, published to thermo as `bond_strain`, `bonds_broken`, `bond_missing` |
| `breakage::*`, `plasticity::*` | criterion and plasticity menus |

## Configuration

```toml
[bonds]
auto_bond = true
bond_tolerance = 1.001
bond_radius_ratio = 1.0
ghost_cutoff_multiplier = 2.5   # MPI: ghost skin to cover bond + 1-3 reach

# Material mode (paper-standard beam theory):
youngs_modulus = 1.0e9     # E (Pa) → K_n, K_bend
shear_modulus  = 4.0e8     # G (Pa) → K_t, K_tor
# Direct overrides used when the matching modulus is absent:
# normal_stiffness / shear_stiffness / twist_stiffness / bending_stiffness

# Critical-damping ratios (raw *_damping overrides also available):
beta_normal = 0.05
beta_shear  = 0.05
beta_twist  = 0.05
beta_bending = 0.05

seed = 0   # per-bond threshold RNG seed

[bonds.breakage]
kind    = "combined_stress"
tensile = { kind = "weibull", mean = 5.0e7, m = 8.0, l_calib = 0.020 }
shear   = { kind = "constant", value = 3.0e7 }
# Threshold kinds: "constant", "weibull", or "crack_band" (length-rescaled).
# Crack-band example (deterministic, mesh-invariant; use with a strain criterion):
# tensile = { kind = "crack_band", value_ref = 0.05, l_ref = 0.020, eps_yield = 0.0, l_min = 0.0 }

[bonds.plasticity.bending]
kind         = "guo_bending"
yield_stress = 1.23e8
# Or the Guo 2018 trilinear envelope (needs [bonds].youngs_modulus):
# kind         = "guo_trilinear"
# yield_stress = 1.23e8

[bonds.plasticity.axial]
kind               = "piecewise"
breakpoint_strains = [0.01, 0.02, 0.03]
slope_multipliers  = [0.5, 0.1, 0.0]
# length_calibration = 0.020   # optional crack-band length regularization (m)

# Load explicit bonds instead of auto-bonding:
# file = "bonds.lammps"
# format = "lammps_data"
```

## Usage

This crate exports `DemBondPlugin` (and `BondConfig`, the breakage/plasticity
menus, the history stores). It depends only on `dirt_atom`, `soil_core`, and
the GRASS app/scheduler — **not** on `dirt_core`. A full DEM application is
assembled through the `dirt_core` umbrella crate, which re-exports
`DemBondPlugin` from its prelude alongside the core/granular plugin groups:

```rust
use dirt_core::prelude::*;          // re-exports DemBondPlugin, CorePlugins, GranularDefaultPlugins

let mut app = App::new();
app.add_plugins(CorePlugins)        // atom data, neighbor list, dump/thermo I/O
    .add_plugins(GranularDefaultPlugins) // Hertz-Mindlin contact + Verlet integration
    .add_plugins(DemBondPlugin);    // bond forces (suppresses contact on bonded pairs)
app.start();
```

If you depend on `dirt_bond` directly (without the umbrella), import the plugin
from this crate instead — `use dirt_bond::DemBondPlugin;` — and supply your own
core/granular plugins. `DemBondPlugin` needs the granular contact plugin present
so non-bonded pairs still interact, and contact on bonded pairs is suppressed
via `soil_core`'s `BondStore`.

## Thermo output

`BondMetrics` publishes three keys each thermo step:

| Key | Meaning |
|-----|---------|
| `bond_strain` | Mean axial strain `δ/r₀` over all live bonds this step (0 if none) |
| `bonds_broken` | Cumulative count of bonds broken since the run started |
| `bond_missing` | Bonds skipped this step because the partner atom was not visible (ghost cutoff too small / MPI split through a bond) |

## Concepts

**Bond vs. contact.** A bonded pair is registered in `soil_core`'s `BondStore`,
and the granular contact plugin skips any pair found there — so a bond is the
*sole* interaction between its two particles until it breaks. On breakage the
pair is removed from `BondStore` and normal Hertz/Hooke contact resumes
automatically (relevant for overlapping bonded spheres, which would otherwise
double-count force).

**Sign convention.** Each bond is owned by the **lower-tag** atom. Forces are
applied as `+F` on atom *i* (lower tag) and `−F` on atom *j* (higher tag) —
Newton's third law by construction. Moments follow the same ownership.

**Per-bond history lifecycle.** A bond is created at setup (auto-bond of
touching pairs, or a loaded `Bonds` section). At creation its breakage
thresholds are drawn once from `(tag_i, tag_j, seed)`, so the failure pattern is
decomposition-independent and bit-reproducible. Each step the shear `Δs` and
twist/bend `Δθ` accumulators (and any plastic anchors / max-strain history) are
carried in `BondHistoryStore`, which survives MPI communication and atom
reordering. When the active criterion trips, the bond is removed.

> `bond_type` is parsed and stored on each bond entry but is **not currently
> consumed** by the force model — it is reserved for future per-type bond
> parameter sets.

## Examples

In the repo's `examples/`: `bond_fiber_tensile` (recovers the input
Young's modulus from a pulled fiber), `bond_fiber_tensile_overlap`
(contact-suppression check on overlapping bonded spheres),
`bond_cantilever`, and `bond_mpi_drift`.

```bash
cargo run --release --example bond_fiber_tensile --no-default-features -- \
    examples/bond_fiber_tensile/config.toml
```

## License

MIT OR Apache-2.0
