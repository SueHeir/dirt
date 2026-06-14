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

Each threshold is a `ThresholdDistribution`: `constant` or a length-scaled
`weibull`. Per-bond thresholds are drawn once at creation from the bond's
tag pair plus `seed`, so the breakage pattern is **MPI-decomposition
independent and bit-reproducible**.

## Plasticity

The bending and axial channels are independently optional
(`plasticity::PlasticityConfig`):

- bending: `guo_bending` (elastic–perfectly-plastic cap at
  `M^p = (4/3) σ_0 r_b³`), `guo_trilinear` (Guo Eq. 32), or `piecewise`
- axial: `piecewise` (breakpoints in axial strain, slope multipliers)

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

[bonds.plasticity.bending]
kind         = "guo_bending"
yield_stress = 1.23e8

[bonds.plasticity.axial]
kind               = "piecewise"
breakpoint_strains = [0.01, 0.02, 0.03]
slope_multipliers  = [0.5, 0.1, 0.0]

# Load explicit bonds instead of auto-bonding:
# file = "bonds.lammps"
# format = "lammps_data"
```

## Usage

```rust
use dirt_core::prelude::*;

let mut app = App::new();
app.add_plugins(CorePlugins)
    .add_plugins(GranularDefaultPlugins)
    .add_plugins(DemBondPlugin);
app.start();
```

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
