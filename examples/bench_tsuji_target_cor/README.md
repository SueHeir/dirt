# Physical target-COR calibration for Hertz--Tsuji contact

`[[dem.materials]].restitution` remains DIRT's **legacy raw Tsuji input**: it is
also the restitution argument of LAMMPS `hertz/material ... damping tsuji`.
It is not renamed, and existing input decks keep their behaviour.  In particular,
the raw value `0.5` realizes roughly `0.614` for the no-tensile-cutoff convention.

For a requested physical normal COR `e_target`, use DIRT's maintained conversion:

```bash
cargo run --release --example bench_tsuji_target_cor -- calibrate 0.70
```

It prints `target, e_raw, beta`; `e_raw` is placed in the declarative material
block. The campaign invokes this executable, then checks the returned mapping
against a separately implemented dimensionless Hertz contact ODE:

```toml
[[dem.materials]]
# requested physical target COR = 0.70
restitution = 0.601268762 # generated raw Tsuji input; do not substitute 0.70
```

The force mapping is the one already used by DIRT and LAMMPS:

`beta = alpha(e_raw)/sqrt(5)`, with `alpha` the Tsuji--Tanaka--Ishida (1992)
polynomial. DIRT bisection-inverts its documented dimensionless collision
equation; the campaign separately re-integrates the ODE as an independent
oracle and checks the conversion residual. Its
dimensionless form means the mapping is independent of size, density, stiffness,
and velocity; the sweep below tests that claim in the discretized solver.

Run the maintained campaign (DIRT and LAMMPS when available):

```bash
$BENCH_PYTHON examples/bench_tsuji_target_cor/sweep.py
```

It sweeps physical targets 0.50/0.70/0.90, impact velocity 0.25/1.0 m/s, radius
2.5/5 mm, density 1000/2500 kg/m³, and two Rayleigh-timestep fractions.  Each
case gates DIRT rebound against the declared physical target (absolute error
<= 0.015), verifies the independent ODE inversion residual (<= 0.0005), and
requires LAMMPS parity under the same mapped parameter (absolute COR difference
<= 0.015).
The tracked `data/results.csv` is the aggregate receipt from the committed
48-case run. The benchmark environment supplies matplotlib for labelled and
ticked result figures.

![Target COR and measured DIRT/LAMMPS COR](plots/target_cor.png)

*Measured rebound COR versus the declared physical target; the shaded band is the
±0.015 DIRT pass criterion. PASS: all maintained sweep cases are in band.*

![COR calibration error](plots/calibration_error.png)

*DIRT and LAMMPS target-COR errors over the size, density, velocity, and timestep
sweep; dashed lines are the ±0.015 gate. PASS: no plotted point crosses a gate.*

## Scope and limitation

This calibration applies to a single normal Hertz contact using
`limit_damping = false`, the same convention as the LAMMPS comparison. A different
damping cutoff or a mixed-material pair needs its own contact-level calibration;
per-material raw values are geometrically mixed before DIRT computes `beta_ij`.

Reference: Tsuji, Tanaka & Ishida, *Powder Technology* 71 (1992), 239--250;
LAMMPS `doc/src/pair_granular.rst` (`damping tsuji`).
