# DMT / SJKR Cohesion Benchmark

Validates two attractive contact models that are **physically distinct from JKR**
and exercises DIRT's adhesion-model selection:

- **DMT** (Derjaguin–Muller–Toporov) adhesive **pull-off**, `F = 2·π·w·R*`.
- **SJKR** (simplified-JKR) **area-proportional cohesion**, `F_coh(δ) = c·π·R*·δ`.

`bench_jkr_adhesion` already covers the JKR pull-off constant `(3/2)·π·w·R*`; this
benchmark covers the *other* two branches of `dirt_granular::contact` and shows
that selecting the model changes the measured physics (the DMT/JKR pull-off ratio
is exactly `4/3`).

## Physics

DIRT's `dirt_granular::contact` selects the normal attractive term from the
material columns and the `[dem] adhesion_model` key:

| Model | Selected by | Attractive normal force | Range |
|-------|-------------|-------------------------|-------|
| JKR   | `surface_energy = w`, `adhesion_model = "jkr"` (default) | constant `(3/2)πwR*` | extended (gap regime) |
| **DMT** | `surface_energy = w`, `adhesion_model = "dmt"` | constant `2πwR*` | **overlap only** (no gap regime) |
| **SJKR** | `cohesion_energy = c` (no `surface_energy`) | `c·π·R*·δ`, ∝ contact area | overlap only |

### DMT pull-off — `F = 2·π·w·R*`

DMT adds a **constant** attractive force but, unlike JKR, does **not** extend the
interaction range (`delta_pulloff = 0`), so there is no gap/adhesion-only regime:
the pull-off is realized *inside* geometric overlap. With `restitution = 1` (no
velocity damping) the net normal force during a quasi-static contact is

```
F_n(δ) = (4/3)·E*·√R*·δ^{3/2}  −  2·π·w·R*
```

which is most tensile as the overlap `δ → 0⁺`, tending to exactly `−2πwR*`. The
recorder captures that peak tension as the measured pull-off force. Because the
approach is slow, the last overlap step sits at `δ ≈ 10⁻¹¹ m`, so the residual
Hertz term is negligible and the raw peak equals `2πwR*` to ~0.01%.

### SJKR cohesion — `F_coh(δ) = c·π·R*·δ`

SJKR ("simplified JKR") makes the cohesive force proportional to the circular
contact area `A = π·R*·δ`. It is **linear in overlap** and **vanishes at
separation** — there is no constant pull-off plateau, which is what makes it
qualitatively different from JKR/DMT. Its net normal force (with `restitution = 1`) is

```
F_n(δ) = (4/3)·E*·√R*·δ^{3/2}  −  c·π·R*·δ
```

To isolate the cohesion term without assuming the Hertz coefficient, the sweep
**differences** an SJKR run against a pure-Hertz **baseline** run (same geometry,
`c = 0`) at matched overlap. With `restitution = 1` the normal force is
conservative and a pure function of `δ`, so the shared Hertz term cancels exactly:

```
F_coh(δ) = F_n^{Hertz}(δ) − F_n^{SJKR}(δ) = c·π·R*·δ
```

The benchmark checks that this difference is linear in `δ` with slope `c·π·R*`,
and that the fitted slope is in turn linear in `c` with slope-of-slopes `π·R*`.

> Note: DIRT's DMT/SJKR are *simplified* explicit force models. They reproduce
> the analytical pull-off magnitude (DMT) and area-law cohesion (SJKR) exactly,
> but not the full Maugis contact-area / force–overlap law. It is those closed-form
> magnitudes, not the detailed contact compliance, that are validated here.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 70 GPa | Pa |
| Poisson's ratio ν | 0.22 | — |
| Restitution e | 1.0 | — (no normal damping → conservative force) |
| Density ρ | 2500 | kg/m³ |
| Radius R | 5 | mm |
| Effective radius R\* = R/2 | 2.5 | mm |
| Effective modulus E\* | 36.78 | GPa |
| Approach velocity | 2 | mm/s |

Both spheres are the same material (`R* = R/2`); friction and gravity are zero
(pure normal interaction).

## Parameter Sweep

- **DMT arm** — work of adhesion `w` (= `surface_energy`): 0.1, 0.2, 0.5, 1.0, 2.0, 5.0 J/m²,
  plus one **JKR reference** case at `w = 1.0` to demonstrate model selection.
- **SJKR arm** — cohesion energy density `c` (= `cohesion_energy`): 1e6, 2e6, 5e6,
  1e7, 2e7 Pa (J/m³), differenced against one shared Hertz baseline (`c = 0`).

## Validation Criteria

| Check | Tolerance |
|-------|-----------|
| DMT per-case `F_pulloff` vs `2πwR*` | ≤ 2% relative error |
| DMT linear in `w` (through origin) | R² ≥ 0.999, slope vs `2πR*` ≤ 2% |
| DMT/JKR pull-off ratio vs `4/3` | ≤ 2% (adhesion-model selection) |
| SJKR per-case slope `F_coh(δ)/δ` vs `cπR*` | ≤ 5%, R² ≥ 0.995 |
| SJKR slope-of-slopes vs `πR*` (linear in `c`) | ≤ 5%, R² ≥ 0.999 |

In practice every DMT case matches to **< 0.25%** (linear fit R² = 1.000000), the
DMT/JKR ratio comes out **1.3332** vs `4/3`, and every SJKR slope matches to
**< 0.01%** with R² = 1.000000 — because both attractive terms are exact
closed-form expressions.

## How to Run

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/bench_dmt_sjkr_cohesion/sweep.py

# Or one stage at a time:
python3 examples/bench_dmt_sjkr_cohesion/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_dmt_sjkr_cohesion/sweep.py start       # build, run all cases -> data/*.csv
python3 examples/bench_dmt_sjkr_cohesion/sweep.py graph       # validate + write plots/
```

### LAMMPS comparison

None. DIRT's DMT/SJKR are simplified explicit force models with no exact LAMMPS
counterpart (LAMMPS' `granular ... jkr`/`sjkr` use different force–overlap laws),
so a code-to-code overlay would compare different physics. Validation is against
the closed-form pull-off and area-law expressions only.

### Single case (default config, DMT arm)

```bash
cargo run --release --example bench_dmt_sjkr_cohesion --no-default-features \
  --features precision-double -- examples/bench_dmt_sjkr_cohesion/config.toml
```

## Expected Plots

### DMT pull-off force vs work of adhesion
![DMT pull-off vs w](plots/dmt_pulloff_vs_w.png)

Measured DMT pull-off (markers) lands on the `F = 2πwR*` line (solid); the JKR
line `1.5πwR*` (dashed) is shown for contrast — selecting `adhesion_model = "dmt"`
moves the response onto the upper line (ratio `4/3`).

### SJKR cohesion vs overlap
![SJKR cohesion vs overlap](plots/sjkr_cohesion_vs_overlap.png)

The isolated cohesion force (Hertz − SJKR difference, markers) is linear in the
overlap `δ` and lands on the `c·π·R*·δ` area-law line (solid) for every `c`.

## Assumptions

- **3D simulation**, two equal spheres (one frozen, one free), so `R* = R/2`.
- **No friction, no gravity** — pure normal interaction.
- `restitution = 1` so the normal contact force is **conservative** (no velocity
  damping); this makes the DMT `δ→0` intercept and the SJKR Hertz-cancellation
  exact.
- DIRT's **simplified constant-force DMT** and **area-law SJKR** models: the
  pull-off / cohesion *magnitudes* are exact, not the full contact-area compliance.
- The spheres start ≥ 1.1 diameters apart (insertion overlap check) and approach
  slowly so the contact is finely resolved.

## References

1. B.V. Derjaguin, V.M. Muller, Yu.P. Toporov, "Effect of contact deformations on
   the adhesion of particles", *J. Colloid Interface Sci.* 53(2):314–326, 1975. (DMT)
2. K.L. Johnson, K. Kendall, A.D. Roberts, "Surface energy and the contact of
   elastic solids", *Proc. R. Soc. Lond. A* 324:301–313, 1971. (JKR reference)
3. K.L. Johnson, *Contact Mechanics*, Cambridge University Press, 1985.
4. For the "simplified JKR" (SJKR) area-proportional cohesion model as used in DEM
   codes, see the LIGGGHTS `cohesion model sjkr` documentation.
