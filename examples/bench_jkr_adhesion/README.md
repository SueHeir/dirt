# JKR Adhesion Pull-off Benchmark

Validates DIRT's adhesive contact model by bringing two identical spheres into
adhesive contact and slowly separating them, then measuring the **pull-off
force** — the peak tension the contact sustains before it snaps apart — against
the analytical JKR result. The work of adhesion is swept and the measured
pull-off force is checked for linearity and the correct slope.

## Physics

When two elastic spheres adhere, separating them requires overcoming an
attractive contact force that rises to a maximum, the pull-off force, before the
contact breaks. For a work of adhesion `w` (surface energy) and effective radius
`R*`, the analytical pull-off forces are:

- **JKR** (Johnson–Kendall–Roberts):
  ```
  F_pulloff = (3/2) · π · w · R*
  ```
- **DMT** (Derjaguin–Muller–Toporov):
  ```
  F_pulloff = 2 · π · w · R*
  ```

with `R* = R/2` for two equal spheres (and `R* = R` for a sphere on a flat
wall).

**Which model does DIRT use?** Both. `dirt_granular::contact` implements adhesion
as a *constant attractive force* selected by the `adhesion_model` key:
`F_adh = (3/2)·π·w·R*` for `"jkr"` (the default) and `F_dmt = 2·π·w·R*` for
`"dmt"`, where `w` is the material `surface_energy`. Because the adhesive force
is constant, the pull-off force is exactly that constant. In the **gap regime**
(surfaces separated, no geometric overlap) the net normal force equals this
plateau with no Hertz spring or velocity damping mixed in — that is the clean
pull-off force the benchmark measures. This example uses **JKR** and compares to
`(3/2)·π·w·R*`.

> Note: DIRT's JKR/DMT is a *simplified* explicit constant-force model. It
> reproduces the analytical pull-off magnitude exactly, but it is not the full
> JKR contact-area / force–overlap law (no hysteretic neck, no adhesive
> contribution to the contact stiffness). It is the pull-off magnitude, not the
> detailed contact compliance, that is validated here.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 70 GPa | Pa |
| Poisson's ratio ν | 0.22 | — |
| Restitution e | 0.5 | — |
| Density ρ | 2500 | kg/m³ |
| Radius R | 5 | mm |
| Effective radius R\* = R/2 | 2.5 | mm |
| Approach velocity | 2 | mm/s |

Both spheres are the same material, so `R* = R/2`. Friction is zero (pure normal
pull-off).

## Parameter Sweep

- **Work of adhesion** `w` (= `surface_energy`): 0.1, 0.2, 0.5, 1.0, 2.0, 5.0 J/m²

For each `w` the left sphere is frozen and the right sphere is launched slowly
(2 mm/s) inward; it makes adhesive contact and then separates. The approach is
slow so the few-nm-wide adhesion window is sampled by many steps, and the
constant gap-regime tension is captured cleanly as the measured pull-off force.

## Validation Criteria

| Check | Tolerance | Notes |
|-------|-----------|-------|
| Per-case `F_pulloff` vs `(3/2)πwR*` | ≤ 2% relative error | |
| Linear in `w` (fit through origin) | R² ≥ 0.999 | F_pulloff must scale linearly with w |
| Fitted slope vs `(3/2)πR*` | ≤ 2% relative error | correct JKR coefficient |
| All 6 cases complete | 6/6 | |

In practice every case matches the analytical pull-off to **< 0.001%** and the
linear fit is exact (R² = 1.000000), because the model's adhesion force is an
exact closed-form constant.

## How to Run

Everything is driven by `sweep.py`, which takes one of three commands. With no
argument it runs all three in order.

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/bench_jkr_adhesion/sweep.py

# Or one stage at a time:
python3 examples/bench_jkr_adhesion/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_jkr_adhesion/sweep.py start      # build, run all 6 cases -> data/*.csv
python3 examples/bench_jkr_adhesion/sweep.py graph      # validate against JKR theory + write plots/
```

`graph` reads the existing `data/sweep_results.csv`, so you can re-validate and
re-plot without re-running the simulations.

### LAMMPS comparison

None. DIRT's JKR/DMT is a simplified constant-force model with no exact LAMMPS
counterpart (LAMMPS' `pair_style granular ... jkr` uses the full Maugis
contact-area model with a different force–overlap law), so a code-to-code
overlay would compare different physics. Validation is against the closed-form
pull-off line only.

### Single case (default config)

```bash
cargo run --release --example bench_jkr_adhesion --no-default-features -- examples/bench_jkr_adhesion/config.toml
```

## Expected Plots

### Pull-off force vs work of adhesion
![Pull-off vs surface energy](plots/pulloff_vs_surface_energy.png)

Measured pull-off force (markers) lands on the JKR theory line
`F = (3/2)πwR*` (solid) across the whole sweep — a straight line through the
origin whose slope is the JKR coefficient.

### Force vs separation
![Force vs separation](plots/force_separation.png)

The per-case contact normal force vs surface separation. In the overlap regime
(separation < 0) the force is the Hertz repulsion minus adhesion; in the gap
regime (separation > 0) it is the flat tensile plateau `−F_adh`, whose
magnitude is the pull-off force and grows linearly with `w`.

## Assumptions

- **3D simulation** with two spherical particles; one is frozen, the other free.
- **Two equal spheres**, so `R* = R/2`.
- **No friction** — pure normal pull-off.
- **No gravity** — the spheres interact only through the contact.
- DIRT's **simplified constant-force JKR** model (default `adhesion_model = "jkr"`):
  the pull-off *magnitude* is exact, but the detailed contact compliance is not
  the full JKR area model (see Physics note).
- The spheres start ≥ 1.1 diameters apart (insertion overlap check) and approach
  slowly so the nm-scale adhesion window is well sampled.

## References

1. K.L. Johnson, K. Kendall, A.D. Roberts, "Surface energy and the contact of
   elastic solids", *Proc. R. Soc. Lond. A* 324:301–313, 1971.
2. B.V. Derjaguin, V.M. Muller, Yu.P. Toporov, "Effect of contact deformations on
   the adhesion of particles", *J. Colloid Interface Sci.* 53(2):314–326, 1975.
3. K.L. Johnson, *Contact Mechanics*, Cambridge University Press, 1985.
