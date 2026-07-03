# Granular Column-Collapse Benchmark

Releases a quasi-2D rectangular column of grains on a flat floor and measures the
final runout `L_f` as a function of the initial aspect ratio `a = H / L0`, to test
the experimental column-collapse scaling laws of Lube et al. (2004) and Lajeunesse
et al. (2004). The column is held against a removable vertical **gate** wall while
it settles, then the gate is removed at runtime (`Walls::deactivate_by_name`) and
the column collapses and spreads. The floor is a frictional `dirt_wall` plane,
which is what arrests the spreading deposit and sets the runout.

If a LAMMPS binary is on `PATH`, the same sweep is **also** run in LAMMPS with the
equivalent granular model and overlaid on the runout-vs-aspect-ratio plot as a
code-to-code cross-check (see *LAMMPS cross-check* below). LAMMPS is optional — the
example runs and validates against the experimental laws with no LAMMPS present.

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
| Poisson's ratio ν | 0.25 | — |
| Density ρ | 2500 | kg/m³ |
| Radius R | 1.5 | mm (d = 3 mm) |
| Restitution e | 0.5 | — |
| Friction μ | 0.5 | — (particle–particle **and** particle–wall) |
| Column width L0 | 24 | mm (8 diameters) |
| Slab width W | 9 | mm (3 diameters, quasi-2D) |
| Timestep dt | 4 × 10⁻⁶ | s |

## Parameter Sweep

- **Aspect ratio** `a = H/L0 ∈ {0.5, 1, 2, 3, 4, 5}`, spanning both regimes.
- `a` is varied by the **particle count** (settled column height H) at fixed L0,
  using a packing fraction of 0.60 to size N. Counts run from ~110 (a = 0.5) to
  ~1100 (a = 5), kept modest so each case finishes in a few minutes.
- Each case runs two stages: `settle` (80 000 steps — pack the loosely-inserted
  column against the gate) then `collapse` (200 000 steps — gate removed on the
  first step, column spreads).

## Validation Criteria

| Check | Tolerance |
|-------|-----------|
| Linear-regime exponent (a ≤ 3) vs 1 | within ±0.25 |
| Power-regime exponent (a ≥ 3) vs 2/3 | within ±0.25 |

`graph` fits the runout exponent in each regime by least squares on log–log axes
and exits non-zero if either fit is outside the band. **It currently FAILs** (a
known limitation — see *Status* below): the fitted exponents are **1.57** (linear
regime, target 1.0 — **outside** ±0.25) and **0.58** (power regime, target 2/3 —
inside ±0.25), so `graph` exits 1.

## Status — known FAIL

This benchmark does **not** currently validate to tolerance and its gate exits
non-zero; the documentation here and in `examples/VALIDATION.md` reflect that FAIL
rather than reporting it green. Adding particle–wall sliding friction to `dirt_wall`
(see *Floor friction*) fixed the earlier failure mode where the released column slid
into a one-grain-thick sheet that ran to the domain boundary — the deposit now
arrests as a finite pile — but the fitted **linear-regime exponent (1.57) still lands
outside the ±0.25 band**. The miss is not a wall-friction gap (friction is present and
applied on the floor); it is the noisy few-point log–log fit at these deliberately
modest particle counts (~110–1100), a single seed, and a coarse 6-point aspect sweep
with diameter-scale runout quantization (the `a = 3` and `a = 4` cases fall in the same
runout bin). Recovering clean exponents would need finer/seed-averaged runs. **No
tolerance is loosened to force a pass**; the bench is kept as an honest, visible FAIL.

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
`start` stage runs the **same** sweep in LAMMPS and the `graph` stage overlays it.
This is an optional cross-code check; **LAMMPS never gates the PASS/FAIL** — only
the DIRT-vs-theory exponents do.

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
`x,y,z,radius` CSV the DIRT recorder writes, and the runout `L_f` is extracted with
the **same** `measure_column()` — so the two codes are compared on equal footing.

Both codes track the increasing runout and bracket the reference lines, but neither
lands both fitted exponents inside the ±0.25 band. Fitted exponents (this run): linear
regime (a ≤ 3) — DIRT **1.57**, LAMMPS **1.12** (target 1.0); power regime (a ≥ 3) —
DIRT **0.58**, LAMMPS **1.03** (target 2/3). Only the DIRT-vs-theory exponents gate the
result, and DIRT's linear-regime exponent is outside the band, so the bench FAILs (see
*Status*). The large per-code spread is expected and is exactly why the fit is unreliable
here: with only modest particle counts and a single seed, the regime split at a = 3 is
soft and the few-point log–log fits are noisy; both codes reproduce the increasing,
near-power-law runout and sit on the experimental band, but the exponents do not coincide
to a decimal or reliably fall inside the validation tolerance.

> Caveat: LAMMPS `create_atoms random` rejects overlapping placements, so the loose
> insert column is sized (and the minimum separation relaxed) to seat all `N` grains;
> the placed count matches DIRT to within a few percent. The two codes' loose-fill
> microstructures still differ, which contributes to the point-by-point runout
> differences.

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
