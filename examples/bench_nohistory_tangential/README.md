# bench_nohistory_tangential — history-free (`linear_nohistory`) tangential model

Validates DIRT's new **history-free** tangential contact model against the
documented LAMMPS `pair_granular` velocity-Coulomb law, and demonstrates that it
is genuinely distinct from the default history-based Mindlin path.

## What it checks

Two identical glass spheres are held at a **fixed** normal overlap `δ` (so the
normal force `F_n = k_n δ` is constant) and driven with a prescribed, **reversing**
relative tangential velocity — a triangle *load → unload → reverse* path
`0 → +V → −V → +V → −V`. This is a loading path along which contact history
matters. At every step the pair contact force is evaluated for **both** tangential
models on independent state, and `(model, step, v_t, F_n, F_t, |ξ|)` is recorded.

`sweep.py` then asserts (8 checks, all must pass):

**`linear_nohistory` reproduces the documented law**
(`doc/src/pair_granular.rst`, *tangential linear_nohistory*):

```
F_t = -min(μ·F_n, η_t·|v_t|) · t̂ ,   t̂ = v_t/|v_t|
```

1. zero accumulated displacement — `|ξ| ≡ 0`;
2. `F_t = 0` at every `v_t = 0` crossing (no elastic memory);
3. the sub-cap branch is **linear through the origin** `F_t = η_t v_t`
   (no elastic offset — the signature of a history-free law);
4. the full `min(μF_n, η_t|v_t|)` shape matches every recorded point;
5. both the sliding (Coulomb-capped) and viscous regimes are exercised
   (anti-trivial guard);
6. in the sliding regime `|F_t| = μ·|F_n|` exactly, with `μ` taken independently
   from the input deck — the documented Coulomb critical force.

![History-free tangential force vs LAMMPS reference](plots/nohistory_tangential_lammps.png)

*DIRT `linear_nohistory` measured tangential force against the documented LAMMPS
velocity-Coulomb law. The dashed lines mark the Coulomb cap and the shaded band /
residual lines show the PASS tolerance used by `sweep.py`; current run: PASS.*

**history (Mindlin) is genuinely distinct** — the point of the goal:

7. it **accumulates** a nonzero tangential displacement `|ξ| > 0`;
8. the **distinguishing behavior**: at `v_t = 0` the history model retains a large
   elastic force (`≈ μF_n`), while the history-free model is exactly zero.

## Reference (independent)

- `μ` is the sliding-friction coefficient from `config.toml` (the LAMMPS input
  deck), not a self-consistent DIRT value.
- The Coulomb cap `|F_t| = μ|F_n|` and the velocity-Coulomb shape
  `min(μF_n, η_t|v_t|)` are the documented LAMMPS `linear_nohistory` expressions.
- `η_t` is identified empirically from the viscous branch (a straight line through
  the origin), then the whole `min()` shape is checked against it.

The tangential damping prefactor is DIRT's Mindlin-consistent
`γ_t = 2β√(5/6)√(k_t m_r)`; `linear_nohistory` keeps this damping and drops **only**
the accumulated spring — exactly the LAMMPS behavior (the tangential force is a
pure function of the instantaneous relative tangential velocity).

## Run

```bash
# build + run + validate (PASS/FAIL, exit 0/1)
python3 examples/bench_nohistory_tangential/sweep.py

# re-validate an existing CSV without rebuilding
python3 examples/bench_nohistory_tangential/sweep.py graph
```

## Config (declarative)

`config.toml` sets the material (`E`, `ν`, restitution, `μ`) and the scenario
(radius, fixed overlap, `dt`, velocity amplitude `V`, and the triangle path legs).
No simulation logic lives in the config — it is pure data.

## Enabling the model in your own runs

```toml
[dem]
contact_model = "hertz"          # or "hooke"
tangential_model = "linear_nohistory"   # default is "history" (Mindlin spring)
```
