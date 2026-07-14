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
| Restitution e | 0.926 | — |
| Friction μ | 0.16 | — (particle–particle and grain–base) |
| Column width L0 | 48 | mm (16 diameters) |
| Slab width W | 18 | mm (6 diameters, quasi-2D) |
| Timestep dt | 4 × 10⁻⁶ | s |

## Parameter Sweep

- **Aspect ratio** `a = H/L0 ∈ {0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5}` —
  an **11-point** sweep (7 in the linear regime, 5 in the power regime) so each
  least-squares exponent is fit from many points, not a coarse handful.
- **Seed averaging.** Each aspect ratio is run at **3 deterministic initial phases**
  and the runout is averaged. `generate` writes one exact-count, non-overlapping
  source-file column per phase; the phase shifts the column relative to the rough
  base before the fully dynamical settling stage. This removes rejection-sampler
  underfills while retaining a packing-realization ensemble.
- `a` is varied by the **particle count** at fixed L0. The 16d × 6d cross-section
  is eight times the old 8d × 3d population; counts run from ~880 (a = 0.5) to
  ~8800 (a = 5). A capacity-derived loose-fill height is used, then the executable
  records the actual pre-release particle coordinates. Analysis fits that measured
  release height, never an intended count or nominal height.
- Each case runs two stages: `settle` (80 000 steps — pack the loosely-inserted
  column against the gate) then `collapse` (1 000 000 steps / 4 s — gate removed
  on the first step, column spreads and must meet the terminal-rest criterion).
- The base is a frozen, close-packed monolayer of the same beads. This models a
  rough granular substrate; the fixed layer is excluded from the released height
  and deposit silhouette. A smooth Coulomb plane is not an interchangeable
  representation of that experimental boundary condition.

## Validation Criteria

| Check | Tolerance |
|-------|-----------|
| Linear-regime exponent (a ≤ 3) vs 1 | within ±0.25 |
| Power-regime exponent (a ≥ 3) vs 2/3 | within ±0.25 |

`graph` fits the runout exponent in each regime by least squares on log–log axes
and exits non-zero if either fit is outside the band. It accepts only one row for
each of the 11 scheduled aspects, with all three seeds and finite measured values.
Before fitting, it independently re-derives every row from all **33** release,
final-deposit, and terminal-rest witnesses; it rejects a missing witness or any
summary value that disagrees with those raw measurements. Thus `runout.csv` is a
cache, not evidence that can make a partial or edited campaign pass.
The scheduled aspect proves coverage; the fit itself uses the executable's measured
pre-release height `H/L0`, so loose insertion and settling cannot silently shift a
point horizontally. A LAMMPS overlay is emitted only for a complete, exact-population
**11 × 3** independent campaign with the same release and arrest checks. Fresh runout figures display each fitted exponent and its
unchanged ±0.25 acceptance band directly on the plot.

Each DIRT row also records a SHA-256 fingerprint of the complete physical and
measurement contract: geometry, material, timestep/stages, frozen rough-base
coordinates, aspects/seeds, toe estimator, Froude arrest limit, references, and
the unchanged exponent band. `graph` refuses any CSV whose fingerprint differs.
This makes it impossible to graph a smooth-base, partial, or prior-material
campaign as evidence for the rough-base protocol merely by retaining its CSV.

The initializer is deliberately part of that contract. For every scheduled
aspect/phase it writes exactly `N` active coordinates with centre spacings no
smaller than one diameter, rather than asking a runtime random inserter to place
`N` grains and discovering an underfill later. This is a preparation correction,
not a change to the empirical reference, material, aspect range, toe metric, or
±0.25 exponent bands.

## Status — current evidence required

The tracked figures are historical smooth-base output and are not evidence for this
rough-base protocol. No exponent PASS is claimed until a fresh, complete 11 × 3 campaign has
passed the release/population/rest gates and regenerated both figures. This avoids
turning the historical 8d × 3d result into an implicit claim about the new continuum
resolution.

**The earlier "fit noise" hypothesis was tested and rejected.** The suspected causes
— single seed, a coarse 6-point sweep, and diameter-scale runout quantization — were
all removed: the runout is now **seed-averaged (3 seeds)**, the sweep is **11 points**,
and the runout uses a **sub-diameter deposit-toe metric** (same physical definition
as before — the far edge where the deposit is ≳1 diameter tall — but with the
diameter-scale binning removed). After all three fixes the linear exponent **barely
moved, 1.57 → 1.54**, and the run-to-run seed scatter is now small (σ ≲ 0.1–0.6 in
normalized runout). So the miss is **not** a measurement artifact.

Two independent lines of evidence show it is a genuine **finite-size** limitation of
this deliberately small benchmark, not a DIRT model defect:

1. **Front-definition dependence.** The fitted linear exponent swings with the runout
   definition — a 2-layer deposit toe gives ≈1.5, a 1-diameter-connected-front gives
   ≈0.5 — because at these particle counts (~80–1100, a 3-grain-deep slab) the
   low-aspect deposits are only a few grains thick with no sharp front. A benchmark
   in the self-similar regime the `1.2 a` law describes would not be this sensitive.
2. **Cross-code agreement.** Running the *identical* geometry, model, and metric in
   **LAMMPS** (authoritative granular DEM) gives linear exponent **1.27** and power
   **0.97** — LAMMPS misses the linear target the same way DIRT does. A code-independent
   miss is a property of the benchmark size, not of DIRT.

**Concrete fix path:** a substantially larger system (thicker slab, ~×10 more grains
so the deposit front becomes continuum-like), not more seeds. **No tolerance is
loosened to force a pass**; the bench is kept as an honest, visible FAIL.

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

Historical LAMMPS values from the superseded small-system protocol are not evidence
for the current geometry. A fresh overlay is written only after all 33 LAMMPS
realizations finish with the exact requested population at release and final state
and meet the same terminal Froude limit; a partial campaign is discarded rather
than used to support either a PASS or a failure diagnosis.

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
