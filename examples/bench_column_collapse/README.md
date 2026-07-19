# Granular Column-Collapse Benchmark

Releases a quasi-2D rectangular column of grains on a flat floor and measures the
final runout `L_f` as a function of the initial aspect ratio `a = H / L0`, to test
the experimental column-collapse scaling laws of Lube et al. (2004) and Lajeunesse
et al. (2004). The column is held against a removable vertical **gate** wall while
it settles, then the gate is removed at runtime (`Walls::deactivate_by_name`) and
the column collapses and spreads. The floor is a frictional `dirt_wall` plane,
which is what arrests the spreading deposit and sets the runout.

If a LAMMPS binary is on `PATH`, an intended counterpart is also run and overlaid.
The retained comparison exposed an initial-condition mismatch, described below;
it is not currently a valid code-to-code validation.

## Physics

A column of grains of initial width `L0` and height `H` is released on a flat
floor and spreads to a final runout `L_f`. Experiments collapse onto two regimes
in the aspect ratio `a = H / L0`:

```
(L_f - L0)/L0 ≈ 1.2 · a          for a ≲ 2–3   (low-aspect, linear)
(L_f - L0)/L0 ≈ 1.6 · a^(2/3)    for a ≳ 3      (high-aspect, 2/3 power)
```

The prefactors are experimental and material-dependent; the benchmark validates
the **scaling exponents and the regime change**, not the exact constants. The
exponents fitted per regime should approach **1** (linear) and **2/3** (power).

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 7 × 10⁷ | Pa (softened, keeps `dt` reasonable for a bed) |
| Poisson's ratio ν | 0.245 | — |
| Density ρ | 2500 | kg/m³ |
| Radius R | 1.5 | mm (d = 3 mm) |
| Restitution e | 0.926 | — |
| Friction μ | 0.16 | — (particle–particle **and** particle–wall) |
| Column width L0 | 24 | mm (8 diameters) |
| Slab width W | 9 | mm (3 diameters, quasi-2D) |
| Timestep dt | 4 × 10⁻⁶ | s |

## Parameter Sweep

- **Aspect ratio** `a = H/L0 ∈ {0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5}` —
  an **11-point** sweep (7 in the linear regime, 5 in the power regime) so each
  least-squares exponent is fit from many points, not a coarse handful.
- **Seed averaging.** Each aspect ratio is run at **3 insertion seeds** and the
  runout is averaged; the `seed` field on `[[particles.insert]]` sets the RNG.
  Single-seed packing scatter in this quasi-2D (3-diameter-deep) column is ~±20-25%
  of the runout, so averaging is needed for a stable exponent.
- `a` is varied by the **particle count** (settled column height H) at fixed L0,
  using a packing fraction of 0.60 to size N. Counts run from ~110 (a = 0.5) to
  ~1100 (a = 5), kept modest so each case finishes in a few minutes.
- Each case runs two stages: `settle` (80 000 steps — pack the loosely-inserted
  column against the gate) then `collapse` (200 000 steps — gate removed on the
  first step, column spreads).

## Validation Criteria

| Check | Tolerance |
|-------|-----------|
| Actual particle-count aspect vs requested aspect | within 2% |
| Linear-regime exponent (a ≤ 3) vs 1 | within ±0.25 |
| Power-regime exponent (a ≥ 3) vs 2/3 | within ±0.25 |

`graph` fits the runout exponent in each regime by least squares on log–log axes
and exits non-zero if insertion completeness or either exponent gate fails.

## Status — WITHHELD after initial-condition audit

The old comparison was not an identical DIRT–LAMMPS experiment. DIRT's random
non-overlap insertion exhausted its attempt budget: retained logs show `79/110`
particles at requested `a=0.5`, `140/220` at `a=1`, and `675/1100` at `a=5`.
The final case therefore had actual `a = H/L0 ≈ 3.07`, but the old figure placed
it at the requested label `a=5`. LAMMPS used a taller loose-fill region and placed
the full requested count. It also used a single packing seed while DIRT was
described as a three-seed average.

That mismatch explains the systematic solver separation and corrupts the fitted
DIRT exponent. It does **not** establish a DIRT model defect, cross-code agreement,
or a finite-size law. The harness now computes the aspect ratio from the actual
particle count, plots that value, and rejects any case more than 2% from its
requested aspect.

The next valid campaign must generate a complete initial particle coordinate set
once and feed that identical set to both DIRT and LAMMPS for each declared seed.
Until those matched reruns exist, the benchmark and its exponent claims are
withheld from scientific validation.

## How to Run

Everything is driven by `sweep.py`; with no argument it runs all three stages.

```bash
# Everything: generate configs → build & run → extract runout & plot
python3 examples/bench_column_collapse/sweep.py

# Or one stage at a time:
python3 examples/bench_column_collapse/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_column_collapse/sweep.py start      # build + run DIRT (+ LAMMPS if on PATH) -> data/*.csv
python3 examples/bench_column_collapse/sweep.py graph      # fit exponents + write plots/
```

### Single case (default config)

```bash
cargo run --release --example bench_column_collapse --no-default-features -- examples/bench_column_collapse/config.toml
```

The binary is a thin recorder: it removes the gate on the first `collapse` step
and, after the run, dumps every particle's `(x, y, z, radius)` at rest to
`<output_dir>/data/column_collapse_results.csv`. All runout extraction, regime
fitting, and plotting live in `sweep.py`.

## Expected Plots

### Runout scaling
![Runout scaling](plots/runout_scaling.png)

Normalized runout `(L_f − L0)/L0` vs aspect ratio `a` on log–log axes, with the
two experimental scaling lines (`1.2 a` and `1.6 a^(2/3)`) overlaid. DIRT is shown
as filled circles; if LAMMPS was available, its runout is overlaid as open squares.

### Deposit profile
![Deposit profile](plots/deposit_profile.png)

A side-view (x–z) snapshot of the settled deposit for the representative `a = 2`
case, with the initial column width `L0` marked.

## LAMMPS cross-check

If a LAMMPS binary (`lmp_serial`, `lmp`, `lmp_mpi`, or `lammps`) is on `PATH`, the
`start` stage runs a counterpart sweep and `graph` overlays it. At present this is
a diagnostic only because the initial particle sets are not matched.

The LAMMPS model is the equivalent of DIRT's Hertz–Mindlin granular contact, same
material, geometry, and protocol:

| DIRT | LAMMPS |
|------|--------|
| Hertz normal, `youngs_mod` / `poisson_ratio` / `restitution` | `pair_style granular hertz/material E e nu` |
| Mindlin tangential, `k_t = 8 G* √(R* δ)`, Coulomb `μ` | `tangential mindlin NULL <damp> μ` (NULL → same `k_t`) |
| viscoelastic damping from `e` | `damping tsuji` |
| `[gravity] gz = −9.81` | `fix gravity 9.81 vector 0 0 -1` |
| frictional `dirt_wall` floor / back / side planes | `fix wall/gran granular … zplane/xplane/yplane` (with friction) |
| removable gate (`deactivate_by_name` at collapse) | `fix wall/gran … xplane NULL L0`, then `unfix gate` before the collapse run |

LAMMPS's final deposit is dumped as `(id, x, y, z, radius)`, converted to the same
`x,y,z,radius` CSV, and measured with the same `measure_column()`. That common
measurement is useful, but it does not rescue unmatched inputs: LAMMPS placed the
full requested count while DIRT often placed only 60–72%. The old claim that the
codes missed the experimental exponent "the same way" is withdrawn. A new
cross-code conclusion requires identical file-based initial coordinates and the
same seed ensemble.

## Assumptions

- **Quasi-2D** slab geometry (thin in y, confined by frictionless side walls).
- **Hertz** normal contact with viscoelastic damping; **Mindlin** tangential
  spring + Coulomb friction on both particle–particle and particle–wall contacts.
- Gate release is an instantaneous support removal (no gate-drag artifact), the
  standard idealization for this benchmark.
- The scaling laws are for cohesionless, dry grains; no cohesion/adhesion is used.

## Floor friction

Basal friction is essential here — it is what arrests the spreading deposit and
sets the runout. The floor is a frictional `dirt_wall` plane: `dirt_wall` applies a
**Mindlin tangential (sliding) spring with a `μ|F_n|` Coulomb cap** on plane walls
(using the material's `friction` via `friction_ij`), the wall analogue of the
particle–particle tangential path in `dirt_granular`.

> Historical note: an earlier version of `dirt_wall` resolved only the normal force,
> so a released column slid into a one-grain-thick sheet that ran to the domain
> boundary. Adding particle–wall sliding friction to the core crate (it benefits every
> wall-bounded granular example) fixed that failure mode — the deposit now arrests as a
> finite pile. It did **not** by itself bring the fitted exponents into tolerance: the
> linear-regime exponent is still outside the ±0.25 band, so the bench remains a
> documented FAIL (see *Status*), limited by the noisy few-point fit rather than by wall
> friction.

## References

1. G. Lube, H.E. Huppert, R.S.J. Sparks, M.A. Hallworth, "Axisymmetric collapses
   of granular columns", *J. Fluid Mech.* 508 (2004) 175–199.
2. E. Lajeunesse, A. Mangeney-Castelnau, J.P. Vilotte, "Spreading of a granular
   mass on a horizontal plane", *Phys. Fluids* 16 (2004) 2371–2381.
3. N.J. Balmforth, R.R. Kerswell, "Granular collapse in two dimensions",
   *J. Fluid Mech.* 538 (2005) 399–428.
