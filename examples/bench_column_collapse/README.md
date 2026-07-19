# Granular Column-Collapse Benchmark

Releases a quasi-2D rectangular column of grains on a flat floor and measures the
final runout `L_f` as a function of the initial aspect ratio `a = H / L0`, to test
the planar experimental column-collapse scaling law of Lajeunesse et al. (2004).
The column is held against a removable vertical **gate** wall while
it settles, then the gate is removed at runtime (`Walls::deactivate_by_name`) and
the column collapses and spreads. The floor is a frictional `dirt_wall` plane,
which is what arrests the spreading deposit and sets the runout.

If a LAMMPS binary is on `PATH`, the same sweep is **also** run in LAMMPS with the
equivalent granular model and overlaid on the runout-vs-aspect-ratio plot as a
code-to-code cross-check (see *LAMMPS cross-check* below). Both solvers receive the
same generated active-grain coordinates for every seed, so this isolates solver
behavior rather than conflating it with different random packings. LAMMPS is
optional — the example runs and validates against the experimental laws with no
LAMMPS present.

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

### Reference-geometry boundary

The acceptance reference is specifically the planar experiment of Lajeunesse,
Mangeney-Castelnau, and Vilotte, *Physics of Fluids* **16** (2004), 2371–2381,
doi:[10.1063/1.1736611](https://doi.org/10.1063/1.1736611). Its horizontal-plane,
gate-release geometry is the applicable external reference for this quasi-2D
benchmark. Lube et al.'s axisymmetric experiment is useful background, but its
high-aspect scaling is not interchangeable with the planar `2/3` law and is not
used as a second target. `graph` prints the cited external reference immediately
before evaluating DIRT's two fitted exponents; no DIRT result is used to establish
those target exponents.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 7 × 10⁷ | Pa (softened, keeps `dt` reasonable for a bed) |
| Poisson's ratio ν | 0.245 | — |
| Density ρ | 2500 | kg/m³ |
| Radius R | 1.5 | mm (d = 3 mm) |
| Restitution e | 0.926 | — |
| Friction μ | 0.16 | — (particle–particle and grain–base) |
| Column width L0 | 96 | mm (32 diameters) |
| Slab width W | 30 | mm (10 diameters) |
| Timestep dt | 4 × 10⁻⁶ | s |

## Parameter Sweep

- **Aspect ratio** `a = H/L0 ∈ {0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5}` —
  an **11-point** sweep (7 in the linear regime, 5 in the power regime) so each
  least-squares exponent is fit from many points, not a coarse handful.
- **Seed averaging.** Each aspect ratio is run at **3 deterministic initial
  fabrics** and the runout is averaged. `generate` writes one exact-count,
  non-overlapping, slightly dilated close-packed source column per seed. After a
  common supported first layer, each seed selects a reproducible sequence of
  non-repeating ABC registries and bounded in-plane perturbations. Thus the
  realizations are neither translated copies nor perfect all-contact crystals.
  These controlled realizations do not substitute for a measured
  packing-preparation sensitivity study; a completed dynamic campaign is still
  required before reporting
  statistical agreement with the experiment.
- **Fabric provenance.** `generate` writes a local protocol manifest containing
  a SHA-256 digest of the frozen base and of each of the 33 canonical active
  coordinate sets. `start` recomputes those digests immediately before launch,
  and `graph` recomputes them before accepting raw witnesses. Thus a copied,
  stale, or hand-edited source file cannot be relabelled as a distinct seed
  realization merely because its count and outer envelope remain plausible.
  The manifest is preparation provenance, not dynamic evidence: release,
  terminal-population, and sustained-arrest witnesses are still required.
- **Release-width witness.** The source grid covers the declared 32d column
  width (rather than leaving a gate-side void), and `graph` rejects any run
  whose recorded pre-release active-grain envelopes span less than 95% of
  `L0`. Population alone cannot establish that the controlled initial width
  used to normalize runout was physically present.
- `a` is varied by the **particle count** at fixed L0. The 32d × 10d section
  uses roughly 5,000 active grains at `a = 0.5` and 50,000 at `a = 5`.
  Crucially, that population is not computed from an infinite-packing volume
  fraction: finite edge losses and phase-shifted layers made the earlier source
  prepare a nominal `a = 0.5` column at `H/L0 = 0.610`. The driver now searches
  the exact coordinate generator and uses the largest common seed population
  that fits the scheduled source height (within one fcc layer). The executable
  still records the post-settlement coordinates and fits that measured height;
  the source qualification prevents the *input* from silently relabelling the
  11-point schedule. This is not an adjustment to a fit or acceptance criterion.
- Each case runs two stages: `settle` (80 000 steps — pack the loosely-inserted
  column against the gate) then `collapse` (1 000 000 steps / 4 s — gate removed
  on the first step, column spreads and must meet the terminal-rest criterion).
- The base is a frozen **hexagonally close-packed** monolayer of the same beads,
  extending from the back wall through the full 0.60 m downstream domain. This
  models a rough granular substrate throughout the possible runout; the plane
  below is only a containment support, not a mid-runout boundary change. The
  fixed layer is excluded from the released height and deposit silhouette. A
  smooth Coulomb plane (or a square bead grid) is not interchangeable with it.

## Validation Criteria

| Check | Tolerance |
|-------|-----------|
| Linear-regime exponent (a ≤ 3) vs 1 | within ±0.25 |
| Power-regime exponent (a ≥ 3) vs 2/3 | within ±0.25 |

`graph` fits the runout exponent in each regime by least squares on log–log axes
and exits non-zero if either fit is outside the band. It accepts only one row for
each of the 11 scheduled aspects, with all three seeds and finite measured values.
Before fitting, it independently re-derives every row from all **33** release,
final-deposit, and sustained-rest witnesses; it rejects a missing witness or any
summary value that disagrees with those raw measurements. Thus `runout.csv` is a
cache, not evidence that can make a partial or edited campaign pass.
Each completed realization also carries a content receipt: SHA-256 digests of
the generated configuration, active source, rough-base source, recorder source,
and all four raw witnesses. Reuse and graphing recompute the receipt and reject
a stale or mixed case before it can enter a seed average. This is an ordinary
reproducibility guard, not a cryptographic signature and not physical evidence
in place of a completed campaign.

The optional LAMMPS overlay has a separate per-case receipt. It binds each
external-code release/deposit/arrest witness to the rendered LAMMPS input and
the SHA-256 identity of the executable that produced it. Graphing withholds an
overlay if either changes. This records solver provenance; it does not make a
LAMMPS result an experimental target or replace the required complete campaign.
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
aspect/seed it writes exactly `N` active coordinates in a deterministic,
stacking-disordered, slightly dilated fabric with centre spacings no smaller
than one diameter, rather than asking a runtime random inserter to place `N`
grains and discovering an underfill later. Its triangular layers span the full
release width; a bounded seed-specific in-plane perturbation removes the
special all-contact crystal network without changing the walls or column
envelope.
The independently recorded release witness must confirm coverage before a fit
is allowed. This is a preparation correction, not a change to the empirical
reference, material, aspect range, toe metric, or ±0.25 exponent bands.

### Independent observer

`graph` also invokes `independent_observer.py`. This standard-library program
does **not** import `sweep.py`, use `runout.csv`, or use its gridded toe
measurement. It independently reads all raw release, final, terminal, and
arrest witnesses; checks population, release width, and sustained Froude rest;
it separately counts the generated active-column and rough-base source files
and requires the release witness to contain those exact populations and the
same frozen support. Equal release/final snapshots alone would otherwise let a
runtime underfill evade the independent check. It then identifies the immutable
rough-base coordinates from the release witness and excludes them from the final
deposit before estimating the same stated
two-layer toe via a continuous, footprint-anchored particle silhouette with a
one-diameter physical gap. It fits the same externally cited
planar targets and unchanged ±0.25 bands. A missing witness or disagreement
fails closed, so it is corroboration rather than a criterion engineered from the
driver's result. The implementation was AI-assisted (2026-07-19); it is not a
new experiment and cannot remove the need for a completed 11 × 3 campaign.

## Status — current evidence required

The tracked figures are historical output and are not evidence for this full-length,
continuum-resolution rough-base protocol. No exponent PASS is claimed until a fresh, complete 11 × 3
campaign has passed the release/population/rest gates and regenerated both figures.
This avoids turning historical 8d × 3d or mixed-boundary output into an implicit
claim about the new continuum-resolution rough-base protocol.

Arrest is not inferred from a convenient final frame. DIRT records the maximum
grain speed every 25,000 collapse steps, and the last four samples must each
satisfy `Fr = v_max / sqrt(g d) <= 0.05`; LAMMPS writes and is checked against
the same window. A candidate with a low final speed but a recent velocity burst
is rejected before either a fit or a cross-code overlay can be made. This
strengthens the existing terminal Froude condition and does not alter the
runout estimator, aspect schedule, empirical targets, or ±0.25 bands.

Independent reproduction (2026-07-18): `generate` produced all 33 exact-count,
seed-distinct, non-overlapping sources. The geometry-qualified source spans
`H/L0 = 0.4997` for every `a = 0.5` seed, `2.9798` for `a = 3`, and `4.9915`
for `a = 5`; its minimum active spacing in those boundary cases is at least
1.0128 diameters. An installed LAMMPS 22-Jul-2025 binary independently consumed
the exact `a=0.5`, seed-0 source and its frozen base for a 20-settle/5-collapse
smoke: release and final populations were both 8,062, across gate removal.
This is an initialization and population-integrity check, not a DIRT/LAMMPS
agreement result or an exponent PASS.

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
2. **Cross-code agreement (superseded boundary).** The earlier LAMMPS comparison
   used the same small, mixed-boundary protocol and also missed the target. It
   informs that old protocol, but is not evidence for this full-length rough-base
   case. A fresh matched LAMMPS campaign is required before any code-to-code claim.

**Current correction:** source populations are derived from the actual bounded
coordinate generator, rather than from an incompatible infinite-packing formula.
The original 11 × 3 schedule, experimental reference, two-layer toe rule, and
±0.25 exponent bands are unchanged. **No tolerance is loosened to force a pass**;
a fresh campaign remains required.

## How to Run

Everything is driven by `sweep.py`; with no argument it runs all three stages.

```bash
# Everything: generate configs → build & run → extract runout & plot
python3 examples/bench_column_collapse/sweep.py

# Or one stage at a time:
python3 examples/bench_column_collapse/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_column_collapse/sweep.py start      # build + run DIRT (+ LAMMPS if on PATH) -> data/*.csv
python3 examples/bench_column_collapse/sweep.py graph      # fit exponents + write plots/

# Independent aspect/seed witnesses can run concurrently. This changes only
# scheduling; graphing still requires the complete 11 × 3 physical ensemble.
python3 examples/bench_column_collapse/sweep.py start --jobs 4

# A batch worker may run exactly one declared witness.  This records raw
# release/final/arrest evidence only; it neither fits nor claims a result.
python3 examples/bench_column_collapse/sweep.py start --case 2,0

# `start` resumes an interrupted campaign, but reuses a case only after it
# independently rechecks its release/final populations, release geometry,
# deposit readability, and four-sample Froude arrest witness.  Invalid or
# partial case evidence is deleted and rerun; `graph` still rebuilds every
# summary row from all 33 raw witnesses.  Use --rerun only to intentionally
# regenerate all DIRT witnesses under the same frozen protocol.
python3 examples/bench_column_collapse/sweep.py start --jobs 4 --rerun
```

For a distributed campaign, run `generate` once, dispatch each of the 33
declared `--case ASPECT,SEED` pairs, then run `start` without `--case` to
independently re-admit every completed witness and produce `runout.csv`.
Finally, `graph` is the only command that fits the exponents and can return
PASS.  A successful individual worker is therefore not a partial validation
claim.

### Representative current-protocol case

`config.toml` is the `a=2`, seed-0 member of the current 32d × 10d,
rough-base protocol.  Its two coordinate sources are generated artifacts, so
generate them first; this is deliberately not the old self-contained 8d × 3d
illustration, which had different material, boundaries, and duration.

```bash
python3 examples/bench_column_collapse/sweep.py generate
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
material, geometry, protocol, **and exact generated initial coordinates**:

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
`write_lammps_input` rejects a missing, truncated, malformed, or non-finite DIRT
source file before launching LAMMPS; it no longer synthesizes a separate random
packing.

Historical LAMMPS values from the superseded small-system protocol are not evidence
for the current geometry. A fresh overlay is written only after all 33 LAMMPS
realizations finish with the exact requested population at release and final state
and meet the same terminal Froude limit; a partial campaign is discarded rather
than used to support either a PASS or a failure diagnosis.

`graph` independently reconstructs that optional overlay from all 33 raw LAMMPS
release/deposit/arrest witnesses and requires exact agreement with the cached
summary. If any witness is missing, non-arrested, population-inconsistent, or the
summary disagrees, the overlay is withheld. This is a provenance check only: it
does not make LAMMPS a pass criterion and it does not establish experimental
validity of an unrun DIRT campaign.

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

1. E. Lajeunesse, A. Mangeney-Castelnau, J.P. Vilotte, "Spreading of a granular
   mass on a horizontal plane", *Phys. Fluids* 16 (2004) 2371–2381,
   doi:10.1063/1.1736611. This is the acceptance reference.
2. G. Lube, H.E. Huppert, R.S.J. Sparks, M.A. Hallworth, "Axisymmetric collapses
   of granular columns", *J. Fluid Mech.* 508 (2004) 175–199. Background only;
   its axisymmetric regime must not be substituted for the planar target.
3. N.J. Balmforth, R.R. Kerswell, "Granular collapse in two dimensions",
   *J. Fluid Mech.* 538 (2005) 399–428.

## Authorship and validation boundary

This revision was AI-assisted. The raw DIRT and LAMMPS campaigns have not been
executed as part of this source-only change, so it makes no new PASS, solver-parity,
or experimental-replication claim. The LAMMPS guard was checked by Python syntax
compilation and fail-closed inspection of absent evidence; a complete fresh 11 × 3
campaign in each solver remains necessary to validate the physical result.
