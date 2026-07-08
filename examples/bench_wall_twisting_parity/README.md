# bench_wall_twisting_parity

Checks that wall twisting friction uses the same local contact law for a plane,
cylinder, sphere, and region wall when their local normal and overlap are
equivalent.

The benchmark places one sphere at a 0.5 mm overlap with each wall, gives it a
pure spin about the local normal, and compares the measured torque with the
plane-wall reference

```
tau = mu_tw * |F_n| * R*
```

where `R*` is the particle radius for wall contacts.

Run it with:

```bash
source ~/projects/.build-env
"$BENCH_PYTHON" examples/bench_wall_twisting_parity/sweep.py
```

Latest result:

```
RESULT: PASS (max_rel_err=0.00e+00)
```

![Wall twisting torque parity](plots/wall_twisting_parity.png)

The figure compares measured torque with the plane-wall reference and shows the
relative-error gate used by the benchmark (`max_rel_err < 1e-12`).
