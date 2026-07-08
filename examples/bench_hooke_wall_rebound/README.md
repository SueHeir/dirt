# bench_hooke_wall_rebound - Hooke wall rebound

Validates DIRT's `dirt_wall` normal force when the material table selects
`contact_model = "hooke"`. A single sphere strikes a real plane wall and
rebounds. The sweep compares measured restitution, contact duration, and peak
overlap against the closed-form damped linear oscillator for a rigid wall.

For a wall contact the reduced mass is the particle mass, and the wall is flat
with infinite radius. DIRT therefore uses the same material-table semantics as
the particle-particle Hooke path, with `kn_ij` for the linear spring and
`beta_ij` for damping:

```text
m_eff x'' + gamma_n x' + kn x = 0
gamma_n = 2 beta sqrt(kn m_eff)
beta = -ln(e) / sqrt(pi^2 + ln(e)^2)
```

The exact references are:

| quantity | reference |
|---|---|
| coefficient of restitution | `COR = e` |
| contact duration | `t_c = pi / (sqrt(kn/m_eff) sqrt(1 - beta^2))` |
| peak overlap | damped-oscillator peak with `x(0)=0`, `x'(0)=v` |

## Running

```bash
python3 examples/bench_hooke_wall_rebound/sweep.py
```

The validation gates are COR within 1%, contact time within 2%, peak overlap
within 2%, and velocity-independence of COR and contact time.

## Result

![Measured vs input COR](plots/cor_validation.png)

*Measured wall rebound restitution compared with the exact `COR = e` reference;
the gray band is the PASS tolerance.*

![Contact duration](plots/contact_duration.png)

*Measured contact duration compared with the exact half-period of the damped
linear oscillator; the gray band is the PASS tolerance.*

![Peak overlap](plots/peak_overlap.png)

*Measured peak overlap compared with the closed-form damped-oscillator peak; the
shaded bands are the PASS tolerance.*

## References

- Y. Tsuji, T. Tanaka, T. Ishida, "Lagrangian numerical simulation of plug flow
  of cohesionless particles in a horizontal pipe", Powder Technology 71:239-250
  (1992).
- H.-G. Schafer, S. Dippel, D.E. Wolf, "Force schemes in simulations of
  granular materials", J. Phys. I France 6:5-20 (1996).
