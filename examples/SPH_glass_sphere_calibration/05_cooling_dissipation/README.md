# Granular-Temperature Dissipation Closure — Glass Beads

Calibrates the **inelastic dissipation closure** (the cooling rate) for the
de-fluidization energy balance of a granular gas. A periodic box of glass beads
is given a Maxwellian velocity field and left to cool through inelastic
collisions; the granular-temperature decay is fit to **Haff's law**, the cooling
time `tc` is extracted, and the dissipation coefficient `omega0` is reported and
compared to kinetic theory and (when available) to the same gas run in LAMMPS.

## Physics

For inelastic spheres with a **velocity-independent** coefficient of restitution
(constant `e`), the granular temperature obeys the de-fluidization energy balance

```
dT/dt = -omega0 * sqrt(T) * T
```

whose solution is Haff's law (Haff 1983):

```
T(t) = T0 / (1 + t/tc)^2        →   late-time log-log slope = -2
tc   = 2 / (omega0 * sqrt(T0))
```

so the cooling time `tc` *is* the dissipation closure: it encodes the inelastic
cooling rate `omega0`. DIRT's contact produces a constant `e` (the rebound
benchmark confirms COR is velocity-independent), so the gas follows this `t^-2`
law — **not** the `t^-5/3` viscoelastic law, which applies only to a
constant-coefficient dashpot.

The kinetic-theory dissipation coefficient (Carnahan–Starling pair correlation
`g0`) is

```
omega0 = (4/3) n d^2 g0 sqrt(pi) (1 - e^2)
```

The fitted `omega0` (from `tc` and `T0`) is the calibrated DEM dissipation rate;
its ratio to `omega0_theory` quantifies how the rough-sphere DEM departs from the
smooth kinetic-theory estimate (friction diverts energy into rotation).

## Material Properties

| Property | Value |
|----------|-------|
| Particles | 800 glass spheres, r = 1.1 mm, ρ = 2500 kg/m³ |
| Box | 40 mm cube, fully periodic, φ ≈ 0.07 |
| Young's modulus | 7.0e7 Pa (softened glass-bead modulus) |
| Poisson ratio | 0.245 |
| Restitution | 0.926 (measured glass–glass COR) |
| Sliding friction | 0.16 (measured glass–glass μ_p) |
| Rolling friction | 0 (keeps contact to normal + tangential → exact LAMMPS match) |
| Initial field | Gaussian, σ = 0.5 m/s per component (T₀ ≈ 0.25 m²/s²) |
| Contact | Hertz normal + Mindlin tangential, viscoelastic (tsuji) damping |

The spheres are **rough** (sliding friction → rotational temperature builds up),
but carry **no rolling resistance**, so the DIRT and LAMMPS contact laws are
identical.

## Parameter Sweep

This calibration is a **single representative case** (the canonical glass-bead
gas above); the "sweep" is the time series of granular temperature that the
recorder writes, from which the Haff law is fit. The dissipation closure is the
fitted `tc` / `omega0`.

## Validation Criteria

| Check | Pass condition |
|-------|----------------|
| Finite, non-negative temperatures | all `T` finite and ≥ 0 |
| Cooling | `T_final < T_initial` |
| No energy growth | `max(T) < 1.5 T₀` |
| Haff's law holds | `1/√T` linear in `t`, R² > 0.99 |
| Late-time decay | log-log slope ∈ [−2.3, −1.7] |

The fitted dissipation rate `omega0 = 2/(tc √T0)` and the kinetic-theory
`omega0_theory` are reported alongside as the calibration output (theory
under-predicts `tc` somewhat because friction diverts energy into rotation).

## How to Run

Driven by `sweep.py` (no argument runs all three stages):

```bash
python3 examples/SPH_glass_sphere_calibration/05_cooling_dissipation/sweep.py            # generate → start → graph
python3 examples/SPH_glass_sphere_calibration/05_cooling_dissipation/sweep.py generate   # write sweep/in.lammps
python3 examples/SPH_glass_sphere_calibration/05_cooling_dissipation/sweep.py start       # build, run DIRT + LAMMPS
python3 examples/SPH_glass_sphere_calibration/05_cooling_dissipation/sweep.py graph       # validate vs Haff + plot
```

Single case directly:

```bash
cargo build --release --no-default-features --example sphcal_cooling_dissipation
./target/release/examples/sphcal_cooling_dissipation \
    examples/SPH_glass_sphere_calibration/05_cooling_dissipation/config.toml
```

`graph` re-reads the CSVs, so you can re-validate/re-plot without re-running.

### LAMMPS cross-check

If a LAMMPS binary (`lmp_serial`/`lmp`/`lmp_mpi`/`lammps`) is on `PATH`, `start`
also runs the same gas in LAMMPS (`pair_style granular`, matched Hertz + Mindlin +
tsuji damping, same timestep) and overlays its cooling curve. An 800-particle gas
is chaotic, so the codes are compared by **cooling law** (the `t^-2` slope, the
fitted `tc`, and `omega0`), not trajectory-by-trajectory. LAMMPS is optional —
without it the example runs and validates DIRT alone.

## Expected Plots

![Haff cooling](plots/haff_cooling.png)

*Left:* `T_total/T₀` vs time (log-log) for DIRT and LAMMPS, with the Haff fit
(annotated with the fitted `tc` and dissipation rate `omega0`) and the −2
reference slope. *Right:* DIRT's energy partition into translational and
rotational temperature.

## Assumptions

- Constant (velocity-independent) restitution, so the cooling law is `t^-2`
  (Haff), not the viscoelastic `t^-5/3`.
- Dilute-to-moderate density (φ ≈ 0.07): the Carnahan–Starling `g0` is used as
  the theory reference; the fitted closure is what DIRT actually produces.
- Homogeneous cooling state assumed for the fit window; clustering/inelastic
  collapse at very late times is excluded by the `T/T0 > 1e-3` fit window.
- Rolling resistance is zero so the LAMMPS contact law matches exactly.

## References

1. P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", *J. Fluid Mech.* 134 (1983) 401–430.
2. N.V. Brilliantov, T. Pöschel, *Kinetic Theory of Granular Gases*, Oxford, 2004.
