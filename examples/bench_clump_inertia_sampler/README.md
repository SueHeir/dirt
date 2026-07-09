# Clump Inertia Sampler Determinism

Validates the overlapping-sphere Monte Carlo inertia helper's reproducibility
contract. The example samples a unit sphere with the clump Monte Carlo path so
the exact analytical mass and inertia are known:

```text
m = rho * 4/3*pi*r^3
I = 2/5*m*r^2
```

The left panel checks the actual determinism contract: repeated default-seed
calls and repeated explicit-seed calls must be bitwise identical. The right
panel sweeps explicit seeds and sample counts against the analytical sphere
inertia to show the bounded Monte Carlo spread at the production-sized sample
count used by clump insertion.

```bash
python3 examples/bench_clump_inertia_sampler/sweep.py
python3 examples/bench_clump_inertia_sampler/sweep.py start
python3 examples/bench_clump_inertia_sampler/sweep.py graph
```

![Clump inertia sampler determinism](plots/inertia_sampler_determinism.png)

*Repeatability and seed-spread measurements against the analytical single-sphere
reference. The stdout validation reports the bitwise-repeat and 100 000-sample
inertia-tolerance checks.*
