# Fiber Crossover Coulomb Friction (`bench_fiber_crossover`)

Validates **inter-fiber contact + friction** at a single crossover — the bridge
from one fiber toward yarn / weave mechanics. Two bonded-sphere fibers cross
perpendicularly and touch at one contact. This exercises the *contact* between
the two fibers (Hertz normal + Mindlin tangential + Coulomb cap), which is
distinct from the intra-fiber *bonds* validated by the `fiber_bond` /
`bond_fiber_tensile` examples.

**Test chosen: (B) Coulomb crossover friction.** Option (B) is the cleaner
quantitative check because the sliding limit `F_slide = μ·N` is exact — no
contact-stiffness calibration is needed, and the validation reduces to a single
dimensionless ratio `F_slide / N == μ` that is independent of the precise Hertz
overlap.

## Physics

The crossover is a sphere–sphere contact with normal `n̂ = +ẑ`. The lower fiber
is frozen; the upper fiber is dragged tangentially (`+x`) at constant velocity
with its height held fixed, so the crossover overlap `δ` — hence the normal
load — is a known constant.

- **Hertz normal:** `N = (4/3) E* √(R*) δ^{3/2}`, with `R* = R/2`,
  `1/E* = (1−ν₁²)/E₁ + (1−ν₂²)/E₂`.
- **Mindlin tangential:** an incremental spring `F_t = k_t · s` accumulates with
  tangential displacement `s` (static rise), capped by the **Coulomb slider**
  `|F_t| ≤ μ |N|`. Once `k_t · s` reaches `μN`, the contact slides and
  `F_slide = μ N` (plateau).

**Isolating the crossover force without a per-contact force API.** The recorder
runs in the `Force` schedule phase, *after* `hertz_mindlin_contact` and
`dem_bond_force` but *before* any `PostForce` fix. There, `atoms.force` holds
only contact + bond contributions, so summing `atoms.force` over **all**
upper-fiber spheres makes every intra-fiber bond cancel (Newton's third law),
leaving exactly the crossover contact force on the upper fiber:
`Σ F_z = N` (the measured normal reaction) and `Σ F_x = −F_t` (the measured
tangential reaction). No core crate is modified.

## Material Properties

| Property | Value |
|----------|-------|
| Sphere radius `R` | 1 mm |
| Spheres per fiber | 7 (6 bonds each) |
| Young's modulus `E` | 1.0 × 10⁷ Pa (soft polymer) |
| Poisson ratio `ν` | 0.3 |
| Density `ρ` | 1200 kg/m³ |
| Sliding friction `μ` | 0.4 |
| Drag speed | 1 mm/s |
| Time step `dt` | 2 × 10⁻⁶ s |

Bonds are loaded from an explicit list (`bonds.data`, intra-fiber only).
`auto_bond` is **disabled**: the two crossover spheres touch within the
auto-bond tolerance, so auto-bonding would weld the fibers together and destroy
the contact-friction test.

## Parameter Sweep

The imposed crossover overlap is swept — `δ ∈ {8, 12, 16, 20, 26, 32, 40} µm`
(set by the upper fiber's height `z = 2R − δ`) — producing normal loads `N`
spanning ~3.7 × 10⁻³ to ~4.1 × 10⁻² N (an 11× range). For each case the
tangential force is dragged into its sliding plateau and `F_slide` is averaged
over the plateau window.

## Validation Criteria

- **Per-case ratio:** `|F_slide/N − μ| ≤ 0.05` for every case.
- **Slope:** the least-squares slope of `F_slide` vs `N` (through the origin)
  must satisfy `|m − μ| ≤ 0.06`.

`graph` prints a PASS/FAIL table and exits non-zero on FAIL.

## How to Run

```bash
# Full pipeline (generate configs, run all cases, validate + plot):
python3 examples/bench_fiber_crossover/sweep.py

# Or stage by stage:
python3 examples/bench_fiber_crossover/sweep.py generate
python3 examples/bench_fiber_crossover/sweep.py start
python3 examples/bench_fiber_crossover/sweep.py graph

# Single representative case (standalone):
cargo run --release --example bench_fiber_crossover --no-default-features -- \
    examples/bench_fiber_crossover/config.toml
```

## Expected Plots

- `plots/fslide_vs_N.png` — `F_slide` vs `N` across the sweep: a straight line
  through the origin whose slope equals `μ`, overlaid with the `μN` Coulomb
  limit and the fitted slope.
- `plots/ft_vs_displacement.png` — tangential force vs tangential displacement
  for the representative case: the linear static rise followed by the `μN`
  plateau (the curve tracks the instantaneous `μN` line through sliding).

## Assumptions

- **Quasi-static drag.** The drag speed (1 mm/s) is far below the contact
  signal speed, so inertial / damping contributions to the tangential force are
  negligible and the plateau reflects the Coulomb limit.
- **Single, persistent crossover contact.** The drag distance (60 µm) keeps the
  same sphere pair engaged (it separates only past ~280 µm) while the Mindlin
  spring saturates within ~7 µm.
- **Geometric drift.** As the upper sphere slides off-centre, the contact normal
  tilts slightly, so the measured `Σ F_x` is the `x`-projection of the
  tangential force and reads a few percent below `μN`. The validation uses the
  *measured* `N` (`Σ F_z`), so `F_slide/N` stays close to `μ` regardless.
- Rolling and twisting friction are off; only sliding friction is exercised.

## References

- Mindlin, R. D. & Deresiewicz, H. "Elastic spheres in contact under varying
  oblique forces." *J. Appl. Mech.* **20** (1953) 327–344.
- Johnson, K. L. *Contact Mechanics.* Cambridge University Press (1985).
- Coulomb friction limit: `F_t ≤ μ F_n`.
