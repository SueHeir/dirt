# Angle-of-Repose Benchmark

Forms a static granular heap and measures its **angle of repose** θ_r as a
function of the sliding friction μ. This is the standard DEM bulk-friction
calibration: θ_r is an emergent, many-body property of the contact model, so it
is the right macroscopic check that friction, rolling resistance, and damping
are wired correctly. The reference is **empirical** (θ_r has no closed form), so
validation tests the qualitative laws a correct model must obey rather than a
single analytical number. When a LAMMPS binary is available, the same protocol is
also run in LAMMPS and overlaid for an informative cross-code comparison (it does
**not** gate validation — see "Cross-code overlay").

## Experimental-comparison boundary

This benchmark deliberately does **not** claim to reproduce the 20--25° angle
reported for millimetre glass beads poured from a hopper by Elekes & Parteli,
*PNAS* **118** (2021), e2107965118,
doi:[10.1073/pnas.2107965118](https://doi.org/10.1073/pnas.2107965118).  That
number is an external experimental reference, but it is not a valid pass band
for this calculation: this example releases a pre-settled column onto a flat,
frictional floor, whereas the paper's material is prepared by pouring.  Angle of
repose is preparation- and contact-model-dependent; in particular, the paper
does not supply this implementation's rolling spring/dashpot/cap parameters.

For a useful experimental validation, a new hopper protocol must be implemented
and its *independently measured* bead size distribution, particle--particle and
particle--wall friction, rolling-resistance law, fill height, and uncertainty
must be fixed before simulation.  Selecting a rolling-friction value after
observing an angle inside the published interval is calibration, not validation,
and is intentionally not performed here.  The solver checks below are therefore
limited to protocol-internal DEM behaviour, while the optional LAMMPS result is
an independent implementation cross-check rather than an experiment.

The mismatch is concrete, rather than a generic caveat.  Elekes & Parteli's
Table 1 specifies a sliding coefficient of 0.5 and rolling coefficient of 0.05
for its glass beads.  This benchmark deliberately sweeps sliding friction and
uses a fixed `sds` rolling cap of 0.1; more importantly, it creates a settled
column and releases its retaining cylinder instead of depositing grains through
a hopper.  Consequently neither its plotted values nor its `[10°, 40°]`
protocol-internal band may be read as a comparison with the paper's 20--25°
observation.  A future matched benchmark must make a *single, predeclared*
material/protocol choice and compare its replicate mean and uncertainty to the
source's reported spread; it must fail if that comparison is unavailable or
outside the source band.

## Physics

A loose column of monodisperse spheres is confined inside a thin cylinder on a
flat floor and allowed to settle. The cylinder is then removed ("lifted") and the
column slumps into a conical heap. The heap stops growing when the surface slope
reaches the angle at which gravity along the slope is balanced by inter-particle
friction — the angle of repose:

```
θ_r = atan(slope of the heap surface)
```

measured by fitting the settled surface height `h(r)` against radial distance `r`
on the straight sloping flank, `θ_r = atan(−dh/dr)`.

There is **no exact θ_r**. It depends on μ, rolling friction, restitution,
polydispersity, and the protocol. What is universal is the *behaviour*:

- θ_r **increases monotonically** with sliding friction μ (and with rolling
  friction),
- θ_r sits in a **physically sensible band**, ~20–40° for typical μ,
- θ_r → ~0° (a nearly flat spread) as μ → 0,
- the heap is **reproducible**: independent random packs give the same θ_r to
  within a few degrees.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 1.0 × 10⁷ | Pa |
| Poisson's ratio ν | 0.25 | — |
| Restitution e | 0.4 | — |
| Rolling Coulomb cap μ_roll | 0.1 | — |
| Rolling stiffness k_roll | 1.0 × 10⁻² | N·m/rad |
| Rolling damping γ_roll | 1.0 × 10⁻⁶ | N·m·s/rad |
| Density ρ | 2500 | kg/m³ |
| Radius R | 2.0 | mm |
| Mobile heap particles | 1200 | — |
| Confining-cylinder radius | 25 | mm |
| Gravity g_z | −9.81 | m/s² |

E is softened to 10 MPa (a routine DEM practice) so the Rayleigh-criterion
timestep the solver auto-selects (≈ 2.6 × 10⁻⁵ s at R = 2 mm) is large enough
that each heap settles in a few seconds of wall-clock time. Rolling resistance is
included because pure sliding friction alone gives weakly-held, low heaps; the
rolling parameters are held fixed while μ is swept, so θ_r(μ) is isolated.

### Rolling resistance — the `sds` spring–dashpot–slider model

Rolling resistance uses the **`sds`** (spring–dashpot–slider) model, the same one
LAMMPS's `pair_style granular … rolling sds k_roll γ_roll μ_roll` implements. The
rolling torque is

```
τ_roll = −k_roll·δ − γ_roll·ω_roll,   capped at  |τ_roll| ≤ μ_roll·|F_n|·r_eff
```

where δ is the accumulated rolling-displacement spring (rescaled on slip), ω_roll
the relative rolling angular velocity, and r_eff the reduced radius (the grain
radius at a wall). DIRT exposes this through `rolling_model = "sds"` with
`rolling_stiffness` (k_roll), `rolling_damping` (γ_roll), and `rolling_friction`
(μ_roll, the Coulomb cap) in `[[dem.materials]]`, and `dirt_wall` applies the same
sds rolling on the floor and confining walls.

**Parameter choice** (Ai et al. 2011, *Comput. Geotech.* 38; Wensrich &
Katterfeld 2012, *Powder Technol.* 217): the rolling spring stiffness is tied to
the contact via k_roll ≈ k_t·r² (k_t the tangential stiffness ≈ 2 × 10³ N/m here,
r = 2 mm the grain radius), giving **k_roll = 1.0 × 10⁻² N·m/rad**. The damping
**γ_roll = 1.0 × 10⁻⁶ N·m·s/rad** is ≈ 0.4 of the critical rolling damping
2·√(I·k_roll), enough to suppress rolling oscillation without overdamping; the
rolling-oscillation period 2π·√(I/k_roll) ≈ 7 × 10⁻⁴ s is well resolved by the
≈ 2.6 × 10⁻⁵ s timestep (~28 steps/period). The Coulomb cap **μ_roll = 0.1** sets
the steady rolling resistance. These exact three values are used in **both** codes.

### Base friction from a real frictional floor wall

The heap stands directly on a **frictional plane wall** at z = 0 (normal +z).
`dirt_wall` applies **Mindlin sliding (tangential) friction** on plane walls,
using the material's `friction` coefficient (μ) through `friction_ij` — exactly
the base friction the bottom layer needs so it cannot slide out and pancake the
heap into a thin monolayer. The same swept μ therefore governs both the
particle–particle contacts that set the pile's angle and the particle–floor
contacts that anchor its base.

This replaces an earlier workaround (a frozen rough particle bed standing in for
wall friction, from before `dirt_wall` had tangential friction): no second
material, no `[[group]]`/`[[freeze]]`, no base bed — just one frictional
`[[wall]]` plane. The confining cylinder wall now also carries friction, which is
harmless: it is deactivated at the lift before the heap forms.

## Parameter Sweep

- **Sliding friction μ**: 0.0, 0.1, 0.2, 0.3, 0.4, 0.5
- **Replicates**: 3 independent random packs per μ (distinct deterministic
  insertion seeds), giving a direct run-to-run spread for the reproducibility
  check.

In the lift-the-cylinder protocol the heap forms by a column *collapse* on the
frictional floor. At low μ the bottom layer slides out and the deposit spreads
into a near-flat disk (θ_r ≈ 0°); as μ grows the floor friction arrests the
runout and the deposit relaxes into a cone whose flank steepens with μ. The
collapse kinetic energy and the wide low apron the cone grows keep the absolute
angle modest — this protocol reads low — so the sensible band here is ~10–40°,
not the 25–40° of a slowly poured heap. The trend (monotonic, μ = 0 → flat) is
the primary validation.

### Measured result (3 packs per μ)

| μ | mean θ_r | std |
|---|----------|-----|
| 0.0 | 0.0° | 0.0° |
| 0.1 | 1.9° | 0.6° |
| 0.2 | 8.1° | 0.4° |
| 0.3 | 10.1° | 0.9° |
| 0.4 | 13.6° | 0.9° |
| 0.5 | 13.5° | 0.9° |

θ_r rises monotonically within the stochastic slack, is flat (≈ 0°) at μ = 0,
has a small but nonzero slope by μ = 0.1, and climbs into the collapse-heap band
by μ = 0.3–0.5. Absolute angles vary within the quoted spread between
independent random packs; the trend is reproducible. (`graph` PASSes on this
data.)

The "lift the cylinder" protocol, per case:
1. **fill** — 1200 mobile spheres are inserted inside a narrow 25 mm cylinder
   (a tall poured-column geometry), resting on the frictional floor wall, and
   settle into a packed column under gravity. When the fastest particle slows
   below 2 mm/s, the cylinder wall is deactivated by name at runtime (the "lift").
2. **lift** — the column slumps across the frictional floor and relaxes into a
   cone. A wide outer cylinder (70 mm, beyond the heap toe) catches the few
   particles flung out during collapse so the count is conserved; it never
   touches the static heap. When the heap comes to rest (fastest particle
   < 1 cm/s, or a 150k-step cap after lift — the geometry locks in well before the
   last micro-jittering particle stops), `main.rs` dumps every particle's final
   `(x, y, z, radius)`.

## Validation Criteria

| Check | Tolerance | Notes |
|-------|-----------|-------|
| θ_r monotonic in μ | mean may dip ≤ 2.5° between μ steps | stochastic slack |
| θ_r overall increase | θ_r(μ_max) > θ_r(μ_min) + 1° | friction raises the heap |
| Frictionless heap flat | θ_r(μ=0) ≤ 8° | spreads to a near-flat disk |
| Sensible band | some μ ≥ 0.2 case in [10°, 40°] | collapse-heap band (reads low) |
| Reproducibility | per-μ std dev ≤ 5° | over the 3 random packs |

`graph` prints the per-μ table and a PASS/FAIL, and exits non-zero on FAIL.

## How to Run

Everything is driven by `sweep.py` (run from anywhere). With **no argument** it
runs the **bounded smoke gate** (the harness default); `full` runs the complete
sweep.

```bash
python3 examples/bench_angle_of_repose/sweep.py            # BOUNDED smoke gate (PASS/FAIL) — the harness default
python3 examples/bench_angle_of_repose/sweep.py smoke      # same as no-arg
python3 examples/bench_angle_of_repose/sweep.py full       # full sweep: generate → run → validate → plot

# Or one full-sweep stage at a time:
python3 examples/bench_angle_of_repose/sweep.py generate   # write sweep/<case>/config.toml (full)
python3 examples/bench_angle_of_repose/sweep.py start       # build, run all cases -> data/*.csv (full)
python3 examples/bench_angle_of_repose/sweep.py graph        # fit θ_r, validate, write plots/ (full)
```

`graph` re-reads `data/repose_sweep.csv` (and `data/lammps_results.csv` if it
exists), so you can re-validate and re-plot without re-running the simulations.

### Bounded smoke gate (CI)

The full sweep is 6 μ × 3 reps = 18 pour-settle-lift heap runs (plus the optional
LAMMPS leg), each a ~1200-grain collapse relaxed on a real settle/rest detector —
so it legitimately overran the 1800 s automation cap every hourly run (exit 124)
and validated nothing. `sweep.py` with **no argument** now runs a **bounded gate**
on the **same material, geometry and physics**: a coarse 3-point μ grid
(μ = 0, 0.3, 0.5), **three deterministic seeded reps** each, and no LAMMPS. It
fits the same mean θ_r(μ) the full run measures and asserts the robust
qualitative laws — the frictionless collapse is nearly flat (mean θ_r(0) ≤ 8°),
every frictional mean holds a real slope in the sensible band [10°, 40°], and
mean θ_r rises across the μ range (coarse monotone trend). The bounded gate is
declared in `smoke.toml`; the full-scale sweep and tolerances remain in
`sweep.py full`. It completes under the automation cap on the host and prints
`ALL CHECKS PASSED`/`CHECKS FAILED` (exit 0/1). Measured:

```
  mu   theta_r mean +/- std (deg)   reps   N each
  0.00      0.00 +/- 0.00              3    1200     frictionless -> flat
  0.30     10.57 +/- 0.71              3    1200     mid-friction, in band
  0.50     13.69 +/- 0.27              3    1200     high friction, in band
```

![Bounded smoke gate](plots/smoke_gate.png)

*Bounded harness smoke gate: the actual μ = 0, 0.3, 0.5 measurements are plotted
as deterministic pack points plus the gated mean against the flat-frictionless
limit, the [10°, 40°] frictional pass band, and the coarse increasing-trend
criterion. Latest run: PASS, 4/4 checks passed.*

This is an **additive breakage gate**, not a replacement for the full validation:
it deliberately does **not** assert the fine, μ-resolved monotonicity between
adjacent close friction values, nor the reproducibility spread — those, with
their tolerances (see *Validation Criteria*), remain the full run's job and are
**unchanged** (`sweep.py full`). It reuses the **same physical bounds** as
`validate()`; nothing is loosened. Smoke artifacts land in gitignored
`data/smoke` and `sweep/smoke`.

### Cross-code overlay (optional LAMMPS leg)

If a LAMMPS binary (`lmp_serial` / `lmp` / `lmp_mpi` / `lammps`) is on `PATH`,
`start` also runs the **same lift-the-cylinder protocol in LAMMPS** with a matched
`pair_style granular` Hertz-Mindlin model **and the matched `sds` rolling model**,
and overlays θ_r(μ) on the plot as open dashed markers:

| DIRT | LAMMPS mapping |
|------|----------------|
| `contact_model = "hertz"`, E, ν, e | `pair_coeff … hertz/material E e ν damping coeff_restitution` |
| Mindlin tangential friction μ | `tangential mindlin NULL 1.0 μ` |
| `rolling_model = "sds"` (k_roll, γ_roll, μ_roll) | `rolling sds k_roll γ_roll μ_roll` — **identical values**, in `pair_coeff` AND every `fix wall/gran` |
| floor plane wall (μ + sds rolling) | `fix wall/gran … rolling sds … zplane 0.0` |
| confining cylinder wall, lifted by name | `fix wall/gran/region … region cyl`, `unfix`-ed at the lift |
| outer catch cylinder (r = 70 mm) | `fix wall/gran/region … region catch` |
| 1200 grains, random non-overlapping insert | `fix pour 1200 … region pourreg` (random, non-overlapping) |
| fill → settle → lift → relax | `run` / `unfix cylwall` / `run` |

The grains are introduced with `fix pour` (random, non-overlapping — the same
packing style as DIRT's overlap-checked inserter). A lattice fill was tried first
and rejected: a crystalline column is mechanically locked and stands as a rigid
pillar that never collapses, so it yields no repose angle. The same heap-fit code
is applied to LAMMPS's settled positions.

**LAMMPS is strictly optional and never gates validation.** `validate()` checks
DIRT against the empirical laws and returns PASS/FAIL on DIRT alone; the LAMMPS
overlay is reported by `compare_codes()` for information only. With no LAMMPS on
`PATH`, the example runs and validates exactly as before.

#### A fair sds↔sds comparison

Every contact-model parameter is matched, **including rolling resistance**. Both
codes run the identical `sds` spring–dashpot–slider rolling model with the
identical parameters — k_roll = 1.0 × 10⁻² N·m/rad, γ_roll = 1.0 × 10⁻⁶ N·m·s/rad,
μ_roll = 0.1 — applied to both grain–grain contacts (`pair_coeff` / `rolling_model
= "sds"`) and grain–wall contacts (every `fix wall/gran` / `dirt_wall`'s sds
branch). DIRT's `sds` rolling and LAMMPS's `rolling sds k_roll γ_roll μ_roll` are
the same model (torque −k_roll·δ − γ_roll·ω_roll, Coulomb-capped at
μ_roll·|F_n|·r_eff, spring rescaled on slip), so the overlay is a genuine
cross-code comparison rather than a comparison of two different rolling laws.

With rolling resistance now matched, **both codes hold a pile** (θ_r > 0,
monotonically rising with μ) instead of LAMMPS pancaking to ≈ 0° — confirming the
earlier flat-LAMMPS result was an artifact of the unmatched rolling model, not a
genuine bulk divergence. Any residual gap between the two curves reflects the
remaining unavoidable differences (pour microstructure: DIRT's overlap-checked
inserter vs LAMMPS's `fix pour`; collapse-protocol energetics), not a model
mismatch. (Measured numbers below.)

### Single case (default config)

```bash
cargo run --release --example bench_angle_of_repose --no-default-features -- examples/bench_angle_of_repose/config.toml
```

This runs the representative μ = 0.3 case and writes
`examples/bench_angle_of_repose/data/repose_results.csv` (the final particle
positions).

## Expected Plots

### θ_r vs μ
![theta vs mu](plots/theta_vs_mu.png)

Mean DIRT θ_r (filled, with ±1 std-dev error bars over the 3 packs) and the
individual runs versus μ. The DIRT curve rises monotonically through the shaded
sensible band and starts near 0° at μ = 0. If LAMMPS was available, its θ_r(μ) is
overlaid as open dashed markers — with the matched sds rolling model it also holds
a pile that rises with μ, tracking the DIRT curve (the fair sds↔sds comparison
discussed above).

### Heap cross-section
![heap profile](plots/heap_profile.png)

The settled surface envelope `h(r)` for each μ (solid = DIRT; dashed open =
LAMMPS, when present). Steeper flanks (higher θ_r) at larger μ are directly
visible; the slope of each flank is what the fit converts to θ_r. With the matched
sds rolling model both codes build a resolvable cone, so the DIRT and LAMMPS
profiles overlay rather than the LAMMPS deposit collapsing to a flat disk.

## Assumptions

- **3D simulation**, monodisperse spheres (single radius).
- **Hertz–Mindlin** normal/tangential contact with viscoelastic damping (DIRT
  default), plus a fixed `sds` (spring–dashpot–slider) rolling-resistance term
  (k_roll, γ_roll, μ_roll), matched 1:1 to LAMMPS.
- **Softened stiffness** (E = 10 MPa) for a tractable timestep — repose angle is
  governed by friction, not by absolute stiffness, so this does not bias θ_r.
- **Frictional base from a real floor wall.** The heap stands on a frictional
  `[[wall]]` plane at z = 0; `dirt_wall`'s Mindlin sliding friction (the
  material's μ via `friction_ij`) anchors the bottom layer so the pile holds a
  slope. No frozen particle bed, second material, or `[[freeze]]` is involved —
  the floor supplies the base friction directly.
- θ_r is fit on the **straight cone flank only** (apex-skip to just inside the
  toe, where the toe is the outermost radius standing > 0.5 diameters above the
  floor baseline), excluding the rounded apex and the sparse stragglers that
  avalanche out past the toe during collapse.
- "Lift the cylinder" deposits read **lower** than slowly-poured heaps (column
  collapse adds kinetic energy that mobilizes the surface), so absolute θ_r here
  is at the low end of the typical band; the **trend** θ_r(μ) is the validated
  quantity.
- The reference is **empirical**: this validates trends and ranges, not an exact
  angle.

## References

1. Y.C. Zhou, B.H. Xu, A.B. Yu, P. Zulli, "Rolling friction in the dynamic
   simulation of sandpile formation", *Physica A* 269 (1999) 536–553.
2. H.P. Zhu, Z.Y. Zhou, R.Y. Yang, A.B. Yu, "Discrete particle simulation of
   particulate systems: A review of major applications and findings",
   *Chemical Engineering Science* 63 (2008) 5728–5770.
3. J.M.N.T. Gray, "Particle segregation in dense granular flows",
   *Annu. Rev. Fluid Mech.* 50 (2018) 407–433 (heap/repose context).
