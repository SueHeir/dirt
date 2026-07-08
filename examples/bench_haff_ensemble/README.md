# Haff Ensemble Cooling-Time Validation

Runs the sphere, multisphere-clump, and rod Haff cooling examples as seeded
ensembles and gates the fitted Haff cooling time `tc` against a kinetic-theory
estimate, rather than checking only that the cooling curve has the right shape.

The reference is the Enskog/Haff cooling rate for inelastic hard spheres,

```
tc = 2 / (omega0 sqrt(T0)),    omega0 = (4/3) n d^2 g0 sqrt(pi) (1 - e^2)
```

with Carnahan-Starling `g0 = (1 - phi/2)/(1 - phi)^3`. For the clump and rod
cases, `d` is the bounding collision diameter. Because the rate scales as `d^2`,
that proxy is the dominant geometric uncertainty for non-spherical particles; the
gate still uses the same factor-two envelope for every case so a cooling time
that is only order-of-magnitude correct fails.

## How to Run

```bash
python3 examples/bench_haff_ensemble/sweep.py
```

The driver builds each existing Haff example, writes generated per-seed configs
under `examples/bench_haff_ensemble/data/generated/`, runs two deterministic
seeds, fits `1/sqrt(T)` versus time, runs the existing Haff linearity and
late-slope gates for every seed, and overlays an optional LAMMPS comparison for
the first seed when a compatible LAMMPS binary is available.

## Validation Criteria

| Check | Pass condition |
|-------|----------------|
| Finite, non-negative temperatures | all `T` finite and `>= 0` |
| Cooling | `T_final < T_initial` |
| No energy growth | `max(T) < 1.5 T_initial` |
| Haff's law holds | `1/sqrt(T)` linear in `t`, `R^2 > 0.99` |
| Late-time decay | existing per-bench slope gate remains satisfied |
| Kinetic `tc` | ensemble median fitted `tc` lies within `[0.5, 2.0]` of the Enskog/Haff estimate |

## Expected Plot

![Haff ensemble validation](plots/haff_ensemble.png)

*Cooling curves, Haff-fit residuals, and fitted `tc` distributions for the three
Haff examples. The orange band is the kinetic-theory acceptance interval; the
optional LAMMPS curve is shown where the local LAMMPS binary can run the matched
case. Latest committed run: PASS.*

## References

1. P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", *J. Fluid Mech.* 134 (1983) 401-430.
2. N.V. Brilliantov and T. Poschel, *Kinetic Theory of Granular Gases*, Oxford University Press, 2004.
