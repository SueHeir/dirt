# Rolling-Resistance Decay Benchmark

Validates DIRT's **constant-torque rolling-resistance model** by launching a
single sphere in *pure rolling* (v = ωR) on a flat frictional floor wall under
gravity and measuring how the rolling-resistance couple decelerates it. The
measured deceleration is checked against the **exact** analytical rate for that
model. Where a LAMMPS binary is available, the identical setup is run with
LAMMPS's granular `rolling sds` contact on a flat `wall/gran` floor and overlaid
on the plots as a code-to-code cross-check.

## Physics

A sphere of radius R, mass m, moment of inertia I = (2/5) m R², rolls without
slipping on a flat floor under gravity (normal force F_n = mg). DIRT's
*constant* rolling-resistance model applies a couple that opposes the rolling
spin:

```
τ_r = μ_r · F_n · r_eff           (a pure torque, no associated force)
```

Mindlin static (sliding) friction at the contact enforces the rolling
constraint (v = ωR). Writing Newton's equations for the translation and the
spin and eliminating the (unknown) static friction force gives a **constant**
deceleration:

```
m·a   = f                         (translation; f = static friction)
I·ω̇  = R·f − τ_r                 (spin; lever R·f spins up, couple τ_r down)
v = ωR  ⇒  a = R·ω̇               (rolling constraint)

⇒  a · (I/R + mR) = τ_r
⇒  a = τ_r / ((7/5) m R) = (5/7) · μ_r · g · (r_eff / R)
```

For a **flat wall** the effective contact radius is simply `r_eff = R` (no
curvature correction), so the rate collapses to the **exact**

```
a = (5/7) · μ_r · g
```

The factor **5/7** is the signature of the rolling constraint: friction supplies
*both* the translational deceleration and the spin-down torque, so the inertia
enters as `I/R + mR = (7/5) m R`. Both v(t) and ω(t) decay **linearly** at this
rate while remaining in pure rolling, until the sphere stops.

**The floor is a real `[[wall]]` plane.** `dirt_wall` now carries the full
friction trio on every wall type — normal force, Mindlin **sliding** friction
(material `friction`/`friction_ij`), and **rolling resistance** (material
`rolling_friction`/`rolling_friction_ij`, with `constant` and `sds` models per
`rolling_model`; see `crates/dirt_wall/src/lib.rs::wall_rolling_torque`). A
sphere rolling on a wall plane therefore feels both the static friction that
holds the no-slip constraint and the rolling-resistance couple that decelerates
it. Earlier this benchmark faked the floor with a giant frozen sphere (R_f = 5 m)
because particle–wall contact lacked friction; that workaround carried an
`r_eff/R ≈ 0.999` curvature factor and a ~2 % systematic. With the real flat
wall, `r_eff = R` exactly and that systematic is gone.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 1×10⁸ | Pa |
| Poisson's ratio ν | 0.30 | — |
| Density ρ | 2500 | kg/m³ |
| Restitution e | 0.5 | — |
| Sliding friction μ | 0.5 | — |
| Sphere radius R | 5 | mm |
| Gravity g | 9.81 | m/s² |
| Initial speed v₀ | 0.03 | m/s |

E is deliberately soft: the rolling-decay rate is **independent of E** (it
depends only on μ_r and g), and a soft contact permits a larger stable timestep
(dt = 1×10⁻⁵ s), keeping each case to a few seconds. μ is healthy so the static
friction never saturates and the sphere stays in pure rolling (no sliding).

The sphere has no initial-spin TOML knob, so `main.rs` sets ω = v₀/R once at the
first step (it also seats the sphere on the wall, z = R, with a hair of overlap
so normal contact is live from step 0).

## Parameter Sweep

- **Rolling-friction coefficient μ_r**: 0.02, 0.05, 0.10

Each case runs for up to 40 000 steps (0.4 s) — long enough for the slowest
(smallest μ_r) sphere to roll to a stop. `sweep.py` linear-fits the v(t) decay
over the moving window (vₓ between 5 % and 95 % of v₀, before the sphere first
reaches rest) and compares the fitted deceleration to `a = (5/7) μ_r g`.

## Validation Criteria

| Check | Tolerance | Notes |
|-------|-----------|-------|
| Fitted deceleration vs theory | ≤ 2 % relative error | per μ_r (theory is exact) |
| Pure rolling maintained | \|vₓ − ωR\| ≤ 1 % of v₀ | over the fit window |
| All 3 cases complete | 3/3 | |

**Result (no LAMMPS required):** PASS. Measured vs predicted deceleration:

| μ_r | a_fit (m/s²) | a_pred (m/s²) | rel. err | max slip |
|-----|--------------|---------------|----------|----------|
| 0.02 | 0.1401 | 0.1401 | 0.0 % | 7×10⁻⁷ |
| 0.05 | 0.3504 | 0.3504 | 0.0 % | 5×10⁻⁵ |
| 0.10 | 0.7007 | 0.7007 | 0.0 % | 5×10⁻⁴ |

The fit matches the exact theory to better than 0.1 % for every case — with a
flat wall there is no curvature systematic, so the only residual is the
linear-fit/discretisation floor. The slip stays far below 1 % of v₀, so the
sphere is in pure rolling throughout: it is the **rolling-resistance couple**
(not sliding friction) that decelerates it.

## How to Run

Everything is driven by `sweep.py`, which takes one of three commands. With no
argument it runs all three in order.

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/bench_rolling_decay/sweep.py

# Or one stage at a time:
python3 examples/bench_rolling_decay/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_rolling_decay/sweep.py start       # build, run all 3 cases (DIRT [+ LAMMPS]) -> data/*.csv
python3 examples/bench_rolling_decay/sweep.py graph       # fit + validate vs theory + write plots/
```

`graph` reads the existing `data/sweep.csv` (and `data/sweep_lammps.csv` if
present), so you can re-validate and re-plot without re-running the simulations.

### LAMMPS comparison

If a LAMMPS binary (`lmp_serial`, `lmp`, `lmp_mpi`, or `lammps`) is on `PATH`,
`start` also runs each case in LAMMPS and overlays it on the figures —
**DIRT as filled markers, LAMMPS as open markers**. LAMMPS is optional; without a
binary the benchmark runs DIRT only and still validates against theory.

The LAMMPS model mirrors DIRT's: a single sphere rolling on a flat granular
floor (`fix wall/gran ... zplane 0.0`) under gravity, launched in pure rolling.
LAMMPS has no constant-torque rolling model, so it uses `granular` with
`rolling sds` — a stiff spring–dashpot–slider whose torque saturates at the same
cap `μ_r·F_n·r_eff` (with `r_eff = R` for the flat wall). In steady rolling the
SDS torque sits at that cap, so the two codes reproduce the **same**
deceleration: they agree to within `|Δa| ≤ 4×10⁻⁴ m/s²` across all μ_r. (The SDS
spring can make vₓ oscillate slightly about zero *after* the sphere stops, so the
fit uses only the initial monotone decay.)

### Single case (default config)

```bash
cargo run --release --example bench_rolling_decay --no-default-features -- examples/bench_rolling_decay/config.toml
```

This runs the μ_r = 0.05 case and writes `data/rolling_decay_results.csv`
(columns `t, x, vx, omega`).

## Expected Plots

### Velocity decay
![Velocity decay](plots/velocity_decay.png)

Measured vₓ(t) for each μ_r (solid) overlaid on the predicted line
`v(t) = v₀ − a·t` (dashed). The curves are straight (constant deceleration) and
sit exactly on theory.

### Deceleration vs μ_r
![Deceleration vs mu_r](plots/deceleration_vs_mu_r.png)

Fitted deceleration vs rolling friction — DIRT (filled) and LAMMPS (open) markers
on the `a = (5/7) μ_r g` theory line. The relation is linear in μ_r, as the
model predicts.

## Assumptions

- **3D simulation**, single sphere rolling on a flat `[[wall]]` floor (z = 0).
- **Pure rolling start**: ω = v₀/R is imposed at t = 0 (no sliding transient).
- **Healthy sliding friction** (μ = 0.5) so the rolling constraint is enforced by
  static friction without saturating — the decelerator is rolling resistance,
  not sliding.
- **Flat wall**: `r_eff = R` exactly, so the theory `a = (5/7) μ_r g` carries no
  curvature correction.
- **Constant rolling-resistance model** (`rolling_model = "constant"`), the model
  under test. (DIRT also offers a `"sds"` model, not exercised here.)

## References

1. J. Ai, J.-F. Chen, J.M. Rotter, J.Y. Ooi, "Assessment of rolling resistance
   models in discrete element simulations", *Powder Technology* 206 (2011)
   269–282. (Model A, "constant directional torque".)
2. Y.C. Zhou, B.D. Wright, R.Y. Yang, B.H. Xu, A.B. Yu, "Rolling friction in the
   dynamic simulation of sandpile formation", *Physica A* 269 (1999) 536–553.
