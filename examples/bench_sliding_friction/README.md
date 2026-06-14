# Coulomb Sliding-Friction Benchmark

Validates DIRT's tangential (Coulomb) friction in the **gross-sliding** regime
using the classic "ball thrown onto a rough floor" problem. A single sphere is
launched horizontally with translational speed `v0` and **zero initial spin**;
kinetic friction decelerates its center and spins it up until the contact point
stops sliding, after which it rolls without slipping. The primary reference is the
exact rigid-body result; where a LAMMPS binary is available, the identical setup is
run with LAMMPS's granular `wall/gran` floor and overlaid as a code-to-code
cross-check. LAMMPS is optional — without it, the benchmark still fully validates
against theory.

## Physics

A sphere (radius `R`, mass `m`, moment of inertia `I = (2/5) m R²`) lands on a
flat floor with center velocity `v0` (along `+x`) and `ω = 0`. While the contact
point slides (`v > ωR`), the kinetic-friction force is `F = μ N = μ m g`, giving:

- **Constant deceleration of the center:**
  ```
  a = μ g            (v(t) = v0 − μ g t)
  ```
- **Constant angular spin-up:**
  ```
  α = μ g / ((2/5) R) = (5 μ g) / (2 R)        (ω(t) = α t,  Rω(t) = R α t)
  ```
- **Slip → roll transition** (when `v = ωR`):
  ```
  t* = 2 v0 / (7 μ g)
  ```
- **Final rolling speed** (no-slip, `v = ωR`, conserved thereafter):
  ```
  v_final = (5/7) v0        (independent of μ)
  ```

The kink at `t*` (where the surface speed `Rω` rises to meet the falling center
speed `v`) and the flat `(5/7) v0` plateau are the signatures being checked.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 70 | GPa |
| Poisson's ratio ν | 0.22 | — |
| Restitution e | 0.3 | — (damps the vertical settling mode fast) |
| Density ρ | 2500 | kg/m³ (glass) |
| Projectile radius R | 5 | mm |
| Gravity g | 9.81 | m/s² |
| Floor | `dirt_wall` z-plane at z = 0, normal +z | (perfectly flat) |

## Parameter Sweep

- **Friction μ**: 0.2, 0.3, 0.5, 0.7 at fixed `v0 = 1.0 m/s` — checks `a = μ g`
  and the `t*` scaling.
- **Launch speed v0**: 0.5, 1.0, 1.5 m/s at fixed `μ = 0.5` — checks
  `v_final = (5/7) v0` is linear in `v0` and independent of μ.

Each case runs only as long as needed: `~1.6 t*` plus a short rolling plateau.

## Validation Criteria

| Check | Tolerance | Notes |
|-------|-----------|-------|
| Sliding deceleration `a` vs `μ g` | ≤ 8% relative | OLS fit over the central 70% of the slip phase |
| Final rolling speed `v_final` vs `(5/7) v0` | ≤ 3% relative | mean over the rolling plateau |
| Transition time `t*` vs `2 v0 / (7 μ g)` | ≤ 10% relative | first sample where slip closes |
| All cases produce a clean slip phase | every case | fit must succeed |

`graph` exits non-zero if any case fails.

## How to Run

Everything is driven by `sweep.py`; with no argument it runs all three stages.

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/bench_sliding_friction/sweep.py

# Or one stage at a time:
python3 examples/bench_sliding_friction/sweep.py generate   # write sweep/<case>/config.toml
python3 examples/bench_sliding_friction/sweep.py start      # build + run all cases (DIRT [+ LAMMPS]) -> data/*.csv
python3 examples/bench_sliding_friction/sweep.py graph      # validate vs theory + write plots/
```

### LAMMPS comparison

If a LAMMPS binary (`lmp_serial`, `lmp`, `lmp_mpi`, or `lammps`) is on `PATH`,
`start` also runs each case in LAMMPS and overlays it on the figures —
**DIRT as filled markers, LAMMPS as open markers**. LAMMPS is optional; without a
binary the benchmark runs DIRT only and still validates against theory.

The LAMMPS model mirrors DIRT's exactly: one sphere of the same radius/density
launched horizontally at `v0` with **zero spin** (`velocity ... set v0 0 0`) onto
a flat frictional granular floor (`fix wall/gran granular hertz/material <E> <e>
<nu> tangential mindlin NULL <xt> <mu> ... zplane 0.0 NULL`) under gravity, with
the same material (E, ν, restitution, friction μ). The sliding-phase deceleration
`a` and rolling plateau `v_final` are fit the **same way** as DIRT. Because
`a = μg` and `v_final = (5/7)v0` are identical in both codes, the cross-check is
near-exact: across all cases `|Δa| ≤ 4×10⁻³ m/s²` and `|Δv_final| ≤ 1×10⁻³ m/s`.

### Single case (default config)

```bash
cargo run --release --example bench_sliding_friction --no-default-features -- examples/bench_sliding_friction/config.toml
```

## Expected Plots

- **`slip_to_roll.png`** — left: `v_x(t)` (solid) and `Rω(t)` (dashed) for the μ
  sweep, showing the falling center speed meeting the rising surface speed at the
  per-μ transition (vertical lines) and settling onto the common `(5/7) v0`
  plateau, with the LAMMPS `v_x` overlaid as open circles and `Rω` as dotted
  lines on the same color; right: a single-μ zoom overlaying the theory lines
  `v0 − μg t` and `Rα t` and the predicted `t*`.
- **`decel_vs_mu.png`** — fitted sliding deceleration vs μ, on the line `a = g μ`;
  DIRT filled circles, LAMMPS open squares.
- **`vfinal_vs_v0.png`** — measured rolling-plateau speed vs `v0`, on the line
  `(5/7) v0`, confirming independence from μ; DIRT filled squares, LAMMPS open
  diamonds.

## Assumptions

- **The floor is a real frictional `dirt_wall` z-plane.** `dirt_wall` now applies
  **Mindlin tangential (sliding) friction** with a Coulomb cap on plane walls,
  using the material's `friction` coefficient (`friction_ij`) — mirroring the
  particle–particle Hertz–Mindlin model. The flat wall therefore decelerates the
  sliding sphere at `a = μg` and applies the spin-up torque, with no giant-sphere
  hack. Because the floor is **perfectly flat**, there is no curvature systematic:
  the normal load is exactly `m g` and there is no gravity-along-surface leak, so
  `a = μg` and `v_final = (5/7) v0` are exact.
- **Single particle, 3D, monodisperse projectile.**
- **Gross sliding throughout the slip phase** (no micro-stick): the Mindlin spring
  saturates at the Coulomb cap `μ|F_n|` immediately because `v0 ≫ 0` and `ω0 = 0`.
- **Vertical settling is decoupled.** The projectile is placed just touching the
  pole; any brief normal-mode ring-down (period ~0.25 ms, `e = 0.3`) damps in a few
  timesteps, far shorter than `t*` (40–150 ms). The deceleration fit trims the
  first 15% of the slip window to exclude it.

## References

1. K. L. Johnson, *Contact Mechanics*, Cambridge University Press, 1985 (§ rolling
   and sliding of spheres).
2. D. Kleppner and R. Kolenkow, *An Introduction to Mechanics*, 2nd ed.,
   Cambridge University Press, 2014 (the "billiard ball / bowling ball" slip-to-roll
   problem).
