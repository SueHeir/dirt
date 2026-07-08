# Guo/Curtis Cantilever Bending

This benchmark reproduces the flexible-fiber cantilever bending validation from
Guo, Wassgren, Hancock, Ketterhagen, and Curtis, *Powder Technology* 249 (2013),
386-395, Fig. 3-4. A 10-sphere bonded fiber is fixed at one end with
`[[freeze]]`; a static transverse load is applied to the free-end sphere and the
quasi-static response is compared with Euler-Bernoulli thin-beam theory.

The checked quantities are:

- free-end deflection versus normalized load `F (L-rs)^2 / EI`
- along-fiber deflection distribution at `|y0|/(L-rs) ~= 0.15`
- along-fiber bending-moment distribution at the same load

The benchmark uses its own DIRT simulation output and its own generated plots.
The paper is cited only as the independent validation target.

## Run

```bash
source ~/projects/.build-env
$BENCH_PYTHON examples/bench_curtis_cantilever/sweep.py
```

The sweep writes run-specific TOML files under `runs/`, executes the
`bench_curtis_cantilever` example at five normalized loads, writes summary CSVs
under `data/`, and regenerates the plots below.

## Results

![Tip deflection vs load](plots/tip_deflection_vs_load.png)

*Normalized free-end deflection from DIRT overlaid on the Euler-Bernoulli
small-deflection beam-theory curve. The shaded band is the +/-3% PASS gate.
Latest run: PASS, max relative error 2.988%, with all 9 bonds intact.*

![Moment and deflection profiles](plots/moment_deflection_profiles.png)

*Deflection and bending-moment distributions along the fiber at normalized load
0.45, the analog of the paper's Fig. 4 case. Latest run: PASS, max absolute
profile error 2.828% for deflection and 1.791% for bending moment.*

Current quantitative gate from `data/results.csv` and `data/target_profile.csv`:

| check | latest error | gate |
| --- | ---: | ---: |
| tip deflection curve | 2.988% | 3% relative |
| deflection profile | 2.828% | 3% absolute |
| bending-moment profile | 1.791% | 3% absolute |
| broken bonds | 0 | 0 |

## Reference

Guo, Y., Wassgren, C., Hancock, B., Ketterhagen, W., and Curtis, J. (2013).
Validation and time step determination of discrete element modeling of flexible
fibers. *Powder Technology*, 249, 386-395.
