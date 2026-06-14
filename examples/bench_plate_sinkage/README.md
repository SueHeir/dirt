# bench_plate_sinkage — terramechanics pressure–sinkage (Bekker form)

Validates that DIRT reproduces the **empirical pressure–sinkage law** of
terramechanics: a flat plate pressed into a granular bed develops a pressure that
grows as a **power law in sinkage depth**. This is the bearing relation behind
landing-pad and wheel/track sinkage on regolith and soil — directly relevant to
HLS landing-pad bearing capacity.

Because the Bekker relation is *empirical* (its constants are soil-fit, not
derived from a contact law), this benchmark validates the **form and qualitative
behavior**, not specific constants: p(z) is monotone, well fit by p ∝ zⁿ with a
physically sensible exponent, and wider/deeper trends are sane and repeatable.

## Physics

A flat plate of footprint width `b` is driven vertically downward into a settled
granular bed at a constant slow velocity. The vertical reaction force `F` on the
plate is recorded versus its sinkage depth `z`; the plate pressure is `p = F/A`
with `A = b·L_y`. The Bekker pressure–sinkage relation is

```
p = (k_c / b + k_φ) · zⁿ        (Bekker 1956; Wong, Theory of Ground Vehicles)
```

so for a fixed plate `p ∝ zⁿ`, monotonically increasing, with exponent `n`
typically `≈ 0.5–1.3` for granular soils. The `k_c/b` term is the cohesive/edge
contribution (vanishes for an edgeless plate); `k_φ` is the frictional bulk term.

### How the reaction force is measured

The plate is a **downward-facing plane wall** (`normal_z = -1`) clipped to a
finite footprint with `bound_x_low/high`, driven at a constant `velocity`. The
`dirt_wall` plugin already accumulates the scalar contact force on each plane wall
in `WallPlane::force_accumulator` (it exists to drive the servo controller); the
recorder reads that field directly as the vertical load on the plate. **No core
change is required** — the moving-wall reaction force is public API. The recorder
latches the sinkage datum `z = 0` at first bed contact (reaction above a small
threshold), so the plate's free descent through air before contact is irrelevant.

## Material Properties

| Property | Value | Notes |
|---|---|---|
| Particle radius | ~2 mm (uniform 1.8–2.2 mm) | slight polydispersity for a stable bed |
| Density | 2500 kg/m³ | |
| Young's modulus | 5×10⁶ Pa | **softened sphere** — standard for quasi-static DEM terramechanics; keeps the stable timestep large |
| Poisson ratio | 0.3 | |
| Restitution | 0.3 | heavy damping → bed settles quickly |
| Inter-particle friction | 0.5 (0.8 in one case) | bed shear strength |
| Gravity | 5 g (49.05 m/s²) | enhanced for faster settling (quasi-static test) |
| Plate speed | 0.04 m/s | slow enough to stay quasi-static |
| Bed | ~1600 particles, ~4 cm deep | container floor + side walls, periodic in y |

## Parameter Sweep

Four cases, varying the plate footprint width `b` and the bed friction `μ`:

| case | b (mm) | μ |
|---|---|---|
| `b020_mu05` | 20 | 0.5 |
| `b040_mu05` | 40 | 0.5 (representative `config.toml`) |
| `b060_mu05` | 60 | 0.5 |
| `b040_mu08` | 40 | 0.8 |

The `b` triple tests the width trend; the high-friction case tests bed shear
strength. Each case fits `p ∝ zⁿ` over `0 < z ≤ 30 mm`.

## Validation Criteria

For **every** case (PASS requires all):

- power-law fit exponent in the sensible band **`0.4 ≤ n ≤ 1.6`**,
- power-law fit quality **`R² ≥ 0.90`**,
- pressure **monotonically increasing** with sinkage,
- pressure rises from first to last bin (no flat/decreasing curve).

Plus a cross-case sanity check: at fixed `μ`, the **total load** `F = p·A` does not
decrease as the plate width `b` grows (a wider plate supports more load).

## How to Run

```bash
python3 examples/bench_plate_sinkage/sweep.py generate   # write per-case configs
python3 examples/bench_plate_sinkage/sweep.py start      # build + run all sims -> CSV
python3 examples/bench_plate_sinkage/sweep.py graph      # validate + plot
python3 examples/bench_plate_sinkage/sweep.py            # all three, in order
```

A single representative case can be run directly:

```bash
cargo run --release --example bench_plate_sinkage --no-default-features -- \
    examples/bench_plate_sinkage/config.toml
```

If a LAMMPS binary (`lmp_serial` / `lmp` / `lmp_mpi` / `lammps`) is on `PATH`,
each width case is **also** run in LAMMPS as a cross-code overlay (see below).
LAMMPS is **optional** — the example fully runs and validates against the Bekker
power-law form with no LAMMPS present; PASS/FAIL is decided on DIRT alone.

### LAMMPS cross-code overlay

The LAMMPS leg reproduces the *same* setup with the *same* material in the
GRANULAR package: `pair_style granular` with `hertz/material E e nu tangential
mindlin ... mu damping tsuji` (matching DIRT's `E / nu / e / μ`), the same 5 g
gravity, the same loose-insert-then-settle bed (a coarse cubic lattice cloud,
spacing > grain diameter so there is no initial overlap, settled under gravity),
and a flat plate of the same footprint width `b` pressed straight down at the same
constant velocity. The plate is a **one-grain-thick raft of frozen grains**
(`fix move linear`) spanning `x ∈ [-b/2, b/2]` across the full y slice; its
vertical reaction force is the plate↔bed contact force read with
`compute group/group plate bed pair yes` (z component). `p = F/A` vs sinkage `z`
is processed and fit to `p ∝ zⁿ` the **same way** as DIRT (first-contact datum,
binned, log–log least squares). On the pressure–sinkage figure LAMMPS appears as
**open markers with dashed fits** in the matching per-width color.

Because the DEM contact damping differs (LAMMPS `damping tsuji` has no exact
`e → COR` closed form, just as DIRT's β-damping does not), the codes are not
expected to agree pointwise; the overlay checks that both produce the same
qualitative Bekker power law with sensible exponents over the loading branch.

### Outputs

| path | contents | tracked |
|---|---|---|
| `sweep/<case>/` | per-case DIRT configs + run logs | no (gitignored) |
| `data/plate_sinkage_results.csv` | per-step `(time, sinkage, force)` for the standalone run | no |
| `data/curve_<case>.csv` | binned DIRT `p(z)` per case | no |
| `data/sweep.csv` | fitted DIRT `n`, `R²`, `p_max`, flags per case | no |
| `sweep/<case>/in.lammps`, `lammps_trace.txt` | LAMMPS input + `(plate_z, F_z)` trace (if LAMMPS present) | no |
| `data/lammps_curve_<case>.csv` | binned LAMMPS `p(z)` per case | no |
| `data/lammps_results.csv` | fitted LAMMPS `n`, `R²`, `p_max` per case | no |
| `plots/pressure_sinkage.png` | `p` vs `z` log-log, DIRT fit per case (+ LAMMPS overlay) | **yes** |
| `plots/pressure_sinkage_linear.png` | `p` vs `z` linear (monotone view) | **yes** |

## Assumptions

- **Softened spheres** (`E = 5 MPa`) — the Bekker law is a *macroscopic* bearing
  relation; absolute `k_c/k_φ` depend on grain stiffness/shape and are not
  validated. The validated claim is the power-law *form* and a sensible exponent.
- **Enhanced gravity** (5 g) only accelerates settling; the pressing phase is slow
  enough to remain quasi-static.
- **Edgeless-in-y, finite-in-x plate**: the footprint is bounded in x via
  `bound_x_*` and periodic in y, so it is a 2-D-slice plate strip. Width effects
  are exercised through `b` (the x-footprint).
- The footprint is centered; the bed is wider than the widest plate so the plate
  always presses into bulk, not against the side walls.

## References

- M. G. Bekker, *Theory of Land Locomotion*, Univ. of Michigan Press, 1956.
- M. G. Bekker, *Introduction to Terrain–Vehicle Systems*, Univ. of Michigan
  Press, 1969.
- J. Y. Wong, *Theory of Ground Vehicles*, 4th ed., Wiley, 2008 (Ch. 2,
  pressure–sinkage).

## License

MIT OR Apache-2.0
