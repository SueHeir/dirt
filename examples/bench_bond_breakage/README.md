# bench_bond_breakage — bond breakage and plasticity benchmark

A regression benchmark (`sweep.py` driver, quantitative PASS/FAIL gate) for the
`dirt_bond` **breakage** and **plasticity** machinery. It covers deterministic
breakage/plasticity gates and a seeded statistical Weibull weakest-link gate.

The bench reuses the existing **`fiber_bond`** example binary (no new core code,
no new example binary); `sweep.py` generates the geometry + configs, runs them,
and validates the recorded output against **closed-form theory** — no reference
value is back-fitted and no LAMMPS run is needed.

## Run

```bash
# all three phases (generate -> build+run -> validate); exits non-zero on FAIL
python3 examples/bench_bond_breakage/sweep.py
# or individually
python3 examples/bench_bond_breakage/sweep.py generate
python3 examples/bench_bond_breakage/sweep.py start
python3 examples/bench_bond_breakage/sweep.py graph
```

Wall time ≈ 5 s (release build cached). Via the harness:
`~/projects/automation/bin/run-bench.sh examples/bench_bond_breakage`.

## Group A — Guo 2018 trilinear bending, fully-plastic moment cap

A pinned 11-bead fiber is loaded with the built-in three-step transverse tip
schedule (`fiber_bond/main.rs::apply_three_step_load`, activated by the
`bending_plastic_guo` output-dir tag). The bending channel uses the literal
`guo_trilinear` envelope, which caps the bond moment at the fully-plastic moment
(Guo et al. 2018, *Chem. Eng. Sci.* **175**, 118–129, Eq. 31):

```
M^p = (4/3) · sigma_0 · r_b^3
```

Sweeping `sigma_0` (0.5, 1.0, 1.5 MPa) scales `M^p` linearly
(0.67, 1.33, 2.00 mN·m). The peak bond moment is reconstructed from the recorded
kinematics as `M = K_bend · (θ_bend − θ_p_bend)`.

**Gate:** peak reconstructed moment plateaus at `M^p` for every `sigma_0`
(rel. tol 6 %). Measured ratios are 1.000 for all three cases.

| σ₀ (MPa) | M^p (mN·m) | peak (mN·m) | ratio |
|----------|------------|-------------|-------|
| 0.5      | 0.6667     | 0.6667      | 1.000 |
| 1.0      | 1.3333     | 1.3333      | 1.000 |
| 1.5      | 2.0000     | 1.9988      | 0.999 |

The `sigma_0` values are chosen so `M^p` stays below the moment the fixed tip
schedule can drive into the middle bond, so every case reaches its cap.

## Group B — seeded Weibull weakest-link CDF

Sixty independently seeded axial-stress Weibull breakage realizations run the
same 10-bond fiber. The per-bond tensile thresholds are sampled by the solver
from the configured two-parameter Weibull law (`mean = 5 MPa`, `m = 5`) and
written to `bond_thresholds.csv`. Each run keeps the deterministic weakest-link
check: first-break strain must match `min(thr0)/E` within 5%.

The ensemble then checks the measured first-break strains against the analytical
weakest-link CDF for the minimum of `N_bonds` independent Weibull thresholds:

```
F_min(eps) = 1 - exp[-N_bonds · (E·eps/lambda)^m]
lambda = mean / Gamma(1 + 1/m)
```

**Gate:** all 60 per-seed first breaks match their sampled weakest-link
prediction within 5%, and the empirical CDF passes a Kolmogorov-Smirnov gate
`D <= 0.18`. Latest run: max per-seed error 3.8%, KS `D = 0.075`, PASS.

![Weibull CDF and QQ validation](plots/weibull_cdf_qq.png)

*Empirical first-break strain CDF and QQ plot from 60 seeded DIRT realizations
against the analytical weakest-link Weibull CDF. Latest run: PASS, KS `D =
0.075` below the 0.18 gate.*

## Outputs

* `data/guo_trilinear.csv` — one row per Guo case.
* `data/weibull_cdf.csv` — one row per seeded Weibull realization.
* `plots/guo_trilinear_moment.png` — peak moment vs σ₀, M^p line.
* `plots/weibull_cdf_qq.png` — empirical CDF and QQ plot vs weakest-link Weibull theory.

## References

* Guo, Y. et al. (2018), *Chem. Eng. Sci.* **175**, 118–129 — trilinear bending
  envelope and fully-plastic moment cap `M^p = (4/3)σ_0 r_b³` (Eq. 31).
