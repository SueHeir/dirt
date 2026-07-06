# Column-Collapse Runout — MACRO VALIDATION (SPH glass-sphere calibration)

This is the **macro validation gate** of the SPH glass-sphere calibration: the
end-to-end check that the calibrated material (and, later, the SPH closure) must
reproduce. It releases a quasi-2D rectangular column of the canonical glass grains
on a flat floor and measures the final runout `L_f` as a function of the initial
aspect ratio `a = H / L0`, testing the experimental column-collapse scaling laws of
Lube et al. (2004) and Lajeunesse et al. (2004). It is **not** a fitted closure
parameter — it is the macroscopic flow behavior the upstream micro-calibrations
must predict.

The column is held against a removable vertical **gate** wall while it settles,
then the gate is removed at runtime (`Walls::deactivate_by_name`) and the column
collapses and spreads. The floor is a frictional `dirt_wall` plane, which (together
with rolling resistance) arrests the spreading deposit and sets the runout.

If a LAMMPS binary is on `PATH`, the same sweep is **also** run in LAMMPS with the
equivalent granular model (including rolling resistance) and overlaid on the
runout-vs-aspect-ratio plot as a code-to-code cross-check. LAMMPS is optional — the
example runs and validates against the experimental laws with no LAMMPS present.

## Rolling-friction dependency from 03_angle_of_repose

This gate was meant to consume the **calibrated rolling friction `μ_r`** produced
by `03_angle_of_repose`. That upstream deliverable currently reports **no
transferable `μ_r` closure**: at the measured glass sliding friction (`μ_p =
0.16`) its lift-the-cylinder heap is sliding-limited and no swept `μ_r` reaches
the measured 22-26 degree glass repose band.

This example therefore no longer carries a placeholder rolling friction. It runs
the campaign canonical material value recorded in `../calibration.yaml`
(`rolling_friction = 0.10`) and keeps the existing Lube/Lajeunesse exponent gates
strict. If those gates fail, the failure is reported as a remaining macro
model/protocol limitation rather than hidden behind "pending 03" wording.

The current limitation is consistent with the upstream repose result: rolling
resistance alone is too small to recover the missing macro friction when the
protocol is sliding-limited. A future passing 03 closure should replace the
canonical `0.10` here and re-run this harness without changing the exponent
tolerances.

## Physics

A column of grains of initial width `L0` and height `H` is released on a flat
floor and spreads to a final runout `L_f`. Experiments collapse onto two regimes
in the aspect ratio `a = H / L0`:

```
(L_f - L0)/L0 ≈ 1.2 · a          for a ≲ 2–3   (low-aspect, linear)
(L_f - L0)/L0 ≈ 1.6 · a^(2/3)    for a ≳ 3      (high-aspect, 2/3 power)
```

The prefactors are experimental and material-dependent; the gate validates the
**scaling exponents and the regime change**, not the exact constants. The
exponents fitted per regime should approach **1** (linear) and **2/3** (power).

## Material Properties

Canonical glass-bead (ballotini) material, shared across the SPH calibration.

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 7 × 10⁷ | Pa (softened from ~65 GPa real glass; keeps `dt` reasonable) |
| Poisson's ratio ν | 0.245 | — |
| Density ρ | 2500 | kg/m³ |
| Restitution e | 0.926 | — (measured glass–glass COR) |
| Sliding friction μ | 0.16 | — (measured glass–glass) |
| Rolling friction μ_r | 0.10 | — (campaign canonical; 03_angle_of_repose currently has no transferable closure) |
| Radius R | 1.5 | mm (d = 3 mm) |
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
and exits non-zero if either fit is outside the band.

> This gate is intentionally allowed to fail. The current 03_angle_of_repose
> deliverable has no passing `μ_r` closure to transfer, so a non-zero exit here is
> evidence of a remaining macro limitation, not a placeholder setup.

## How to Run

Everything is driven by `sweep.py`; with no argument it runs all three stages.

```bash
# Everything: generate configs → build & run → extract runout & plot
python3 examples/SPH_glass_sphere_calibration/07_column_collapse/sweep.py

# Or one stage at a time:
python3 examples/SPH_glass_sphere_calibration/07_column_collapse/sweep.py generate
python3 examples/SPH_glass_sphere_calibration/07_column_collapse/sweep.py start
python3 examples/SPH_glass_sphere_calibration/07_column_collapse/sweep.py graph
```

### Single case (default config)

```bash
cargo run --release --example sphcal_column_collapse --no-default-features -- examples/SPH_glass_sphere_calibration/07_column_collapse/config.toml
```

The binary is a thin recorder: it removes the gate on the first `collapse` step
and, after the run, dumps every particle's `(x, y, z, radius)` at rest to
`<output_dir>/data/column_collapse_results.csv`. All runout extraction, regime
fitting, and plotting live in `sweep.py`.

## Expected Plots

### Runout scaling
![Runout scaling](plots/runout_scaling.png)

Normalized runout `(L_f − L0)/L0` vs aspect ratio `a` on log–log axes, with the
two experimental scaling lines (`1.2 a` and `1.6 a^(2/3)`) overlaid. The right
panel shows the actual exponent gates (target ±0.25) and PASS/FAIL verdicts. DIRT
is shown as filled circles; if LAMMPS was available, its runout is overlaid as open
squares.

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
| `rolling_friction` (rolling Coulomb cap `μ_r`) | `rolling sds <k_r> <damp> μ_r` |
| viscoelastic damping from `e` | `damping tsuji` |
| `[gravity] gz = −9.81` | `fix gravity 9.81 vector 0 0 -1` |
| frictional `dirt_wall` floor / back / side planes | `fix wall/gran granular … zplane/xplane/yplane` (with friction + rolling) |
| removable gate (`deactivate_by_name` at collapse) | `fix wall/gran … xplane NULL L0`, then `unfix gate` before the collapse run |

LAMMPS's final deposit is dumped as `(id, x, y, z, radius)`, converted to the same
`x,y,z,radius` CSV the DIRT recorder writes, and the runout `L_f` is extracted with
the **same** `measure_column()` — so the two codes are compared on equal footing.

> Caveat: LAMMPS `create_atoms random` rejects overlapping placements, so the loose
> insert column is sized (and the minimum separation relaxed) to seat all `N` grains;
> the placed count matches DIRT to within a few percent. The two codes' loose-fill
> microstructures still differ, which contributes to the point-by-point runout
> differences.

## Assumptions

- **Quasi-2D** slab geometry (thin in y, confined by side walls).
- **Hertz** normal contact with viscoelastic damping; **Mindlin** tangential
  spring + Coulomb friction on both particle–particle and particle–wall contacts;
  **rolling resistance** with a `μ_r` Coulomb cap.
- Gate release is an instantaneous support removal (no gate-drag artifact), the
  standard idealization for this benchmark.
- The scaling laws are for cohesionless, dry grains; no cohesion/adhesion is used.
- **`μ_r` is not a placeholder**. It is the campaign canonical value used because
  `03_angle_of_repose` currently reports no transferable glass-band closure.

## References

1. G. Lube, H.E. Huppert, R.S.J. Sparks, M.A. Hallworth, "Axisymmetric collapses
   of granular columns", *J. Fluid Mech.* 508 (2004) 175–199.
2. E. Lajeunesse, A. Mangeney-Castelnau, J.P. Vilotte, "Spreading of a granular
   mass on a horizontal plane", *Phys. Fluids* 16 (2004) 2371–2381.
3. N.J. Balmforth, R.R. Kerswell, "Granular collapse in two dimensions",
   *J. Fluid Mech.* 538 (2005) 399–428.
