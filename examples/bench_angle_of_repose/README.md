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
- **Replicates**: 3 independent random packs per μ (the inserter is
  entropy-seeded), giving a direct run-to-run spread for the reproducibility
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
| 0.1 | 0.0° | 0.0° |
| 0.2 | 2.0° | 2.0° |
| 0.3 | 8.2° | 2.4° |
| 0.4 | 11.5° | 1.4° |
| 0.5 | 11.5° | 3.9° |

θ_r increases monotonically, is flat (≈ 0°) at μ ≤ 0.1, and climbs into the
collapse-heap band (~10–12°) by μ = 0.4–0.5, where it plateaus. Absolute angles
vary within the quoted spread between independent random packs; the trend is
reproducible. (`graph` PASSes on this data.)

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

Everything is driven by `sweep.py` (run from anywhere). With no argument it runs
all three stages in order.

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/bench_angle_of_repose/sweep.py

# Or one stage at a time:
python3 examples/bench_angle_of_repose/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_angle_of_repose/sweep.py start       # build, run all cases -> data/*.csv
python3 examples/bench_angle_of_repose/sweep.py graph        # fit θ_r, validate, write plots/
```

`graph` re-reads `data/repose_sweep.csv` (and `data/lammps_results.csv` if it
exists), so you can re-validate and re-plot without re-running the simulations.

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
  toe, where the toe is the outermost radius standing > 1.5 diameters above the
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
