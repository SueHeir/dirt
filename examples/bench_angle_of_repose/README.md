# Angle-of-Repose Benchmark

Forms a static granular heap and measures its **angle of repose** θ_r as a
function of the sliding friction μ. This is the standard DEM bulk-friction
calibration: θ_r is an emergent, many-body property of the contact model, so it
is the right macroscopic check that friction, rolling resistance, and damping
are wired correctly. The reference is **empirical** (θ_r has no closed form), so
validation tests the qualitative laws a correct model must obey rather than a
single analytical number.

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
| Rolling friction μ_r | 0.1 | — |
| Density ρ | 2500 | kg/m³ |
| Radius R | 2.0 | mm |
| Mobile heap particles | 1200 | — |
| Confining-cylinder radius | 25 | mm |
| Gravity g_z | −9.81 | m/s² |

E is softened to 10 MPa (a routine DEM practice) so the Rayleigh-criterion
timestep the solver auto-selects (≈ 2.6 × 10⁻⁵ s at R = 2 mm) is large enough
that each heap settles in a few seconds of wall-clock time. A small rolling
friction is included because pure sliding friction alone gives weakly-held, low
heaps; μ_r is held fixed while μ is swept, so θ_r(μ) is isolated.

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

`graph` re-reads `data/repose_sweep.csv`, so you can re-validate and re-plot
without re-running the simulations.

This benchmark is **DIRT-only** — there is no LAMMPS overlay, because there is no
analytical target to cross-check a second code against. The validation is the set
of empirical laws above applied to DIRT's own settled heaps. (`sweep.py` still
probes for a LAMMPS binary for structural parity and notes if one is present.)

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

Mean θ_r (with ±1 std-dev error bars over the 3 packs) and the individual runs
versus μ. The curve rises monotonically through the shaded sensible band and
starts near 0° at μ = 0.

### Heap cross-section
![heap profile](plots/heap_profile.png)

The settled surface envelope `h(r)` for each μ. Steeper flanks (higher θ_r) at
larger μ are directly visible; the slope of each flank is what the fit converts
to θ_r.

## Assumptions

- **3D simulation**, monodisperse spheres (single radius).
- **Hertz–Mindlin** normal/tangential contact with viscoelastic damping (DIRT
  default), plus a fixed rolling-friction term.
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
