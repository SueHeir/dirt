# Haff Cooling Benchmark — Single Rough Spheres

Validates free cooling of a granular gas against **Haff's law**. A periodic box of
frictional spheres is given a Maxwellian velocity field and left to cool through
inelastic collisions; the granular temperature decay is compared to theory and,
as an independent cross-check, to the same gas run in LAMMPS.

## Physics

For inelastic spheres with a **velocity-independent** coefficient of restitution
(constant `e`), the granular temperature follows Haff (1983):

```
T(t) = T0 / (1 + t/tc)^2          →   late-time log-log slope = -2
```

DIRT's contact model produces a constant `e` (the rebound benchmark confirms COR
is velocity-independent), so the gas must follow this `t⁻²` law — **not** the
`t⁻⁵ᐟ³` viscoelastic law (Brilliantov–Pöschel), which applies only to a
constant-coefficient dashpot whose effective restitution rises as collisions slow.

The cooling time from kinetic theory (Carnahan–Starling pair correlation `g0`):

```
tc = 2 / (ω0 √T0),   ω0 = (4/3) n d² g0 √π (1 − e²)
```

## Setup

| Property | Value |
|----------|-------|
| Particles | 800 spheres, r = 1.1 mm, ρ = 2500 kg/m³ |
| Box | 40 mm cube, fully periodic, φ ≈ 0.07 |
| Material | E = 50 MPa, ν = 0.3, e = 0.9, μ = 0.3 |
| Rolling friction | 0 (so the LAMMPS cross-check matches exactly) |
| Initial field | Gaussian, σ = 0.5 m/s per component (T₀ ≈ 0.25 m²/s²) |
| Contact | Hertz normal + Mindlin tangential, viscoelastic (tsuji) damping |

The spheres are **rough** (sliding friction → rotational temperature builds up),
but carry **no rolling resistance**, which keeps the contact to normal + tangential
so the DIRT and LAMMPS contact laws are identical.

## How to Run

Driven by `sweep.py` (no argument runs all three stages):

```bash
python3 examples/bench_sphere_haff_cooling/sweep.py            # generate → start → graph
python3 examples/bench_sphere_haff_cooling/sweep.py generate   # write data/in.lammps
python3 examples/bench_sphere_haff_cooling/sweep.py start      # build, run DIRT + LAMMPS
python3 examples/bench_sphere_haff_cooling/sweep.py graph      # validate vs Haff + plot
```

`graph` re-reads the CSVs, so you can re-validate/re-plot without re-running.

### LAMMPS cross-check

If a LAMMPS binary (`lmp_serial`/`lmp`/`lmp_mpi`/`lammps`) is on `PATH`, `start`
also runs the same gas in LAMMPS (`pair_style granular`, matched Hertz + Mindlin +
tsuji damping, same timestep) and overlays its cooling curve. An 800-particle gas
is chaotic, so the codes are compared by **cooling law** (the `t⁻²` slope and Haff
fit), not trajectory-by-trajectory. LAMMPS is optional — without it the benchmark
runs and validates DIRT alone.

## Validation Criteria

| Check | Pass condition |
|-------|----------------|
| Finite, non-negative temperatures | all `T` finite and ≥ 0 |
| Cooling | `T_final < T_initial` |
| No energy growth | `max(T) < 1.5 T₀` |
| Haff's law holds | `1/√T` linear in `t`, R² > 0.99 |
| Late-time decay | log-log slope ∈ [−2.3, −1.7] |

The kinetic-theory `tc` is reported alongside the fitted `tc` as a diagnostic
(translational theory under-predicts `tc` somewhat because friction diverts energy
into rotation).

## Expected Plot

![Haff cooling](plots/haff_cooling.png)

*Left:* `T_total/T₀` vs time (log-log) for DIRT and LAMMPS, with the Haff fit and
the −2 reference slope — both codes follow the same cooling law. *Right:* DIRT's
energy partition into translational and rotational temperature.

## References

1. P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", *J. Fluid Mech.* 134 (1983) 401–430.
2. N.V. Brilliantov, T. Pöschel, *Kinetic Theory of Granular Gases*, Oxford, 2004.
