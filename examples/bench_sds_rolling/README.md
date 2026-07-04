# bench_sds_rolling — SDS (spring-dashpot-slider) rolling-resistance model

Validates DIRT's **SDS rolling** contact model (`rolling_model = "sds"`, implemented in
`crates/dirt_granular/src/contact.rs`) against **its own analytical behaviour** in both of
its regimes. `bench_rolling_decay` covers only the `constant`-torque rolling model; this
bench is the missing companion for the spring-dashpot-slider variant.

## Physical setup

Two identical spheres (R = 5 mm) are stacked along +z. The **lower** sphere is a frozen
anchor (`[[freeze]]` — it can neither translate nor spin). The **upper** sphere rests on it
under gravity, seated by `main.rs` at the static Hertz overlap so the normal force is
`F_n = m g` from step 0 (no settling ring-down), and is given a pure **rolling** spin
`ω = (ω0, 0, 0)` about a horizontal axis.

`ω` about x̂ is perpendicular to the contact normal `n̂ = +ẑ`, so its projection onto n̂ is
zero: this is *rolling*, not twisting. Sliding friction is turned **off** (`friction = 0`) so
the tangential surface slip a horizontal spin would otherwise drive cannot add a couple —
the **only** torque ⊥ n̂ is the SDS rolling resistance under test. That couple is a pure
torque (no force), so the sphere never translates; `ω_x` simply spins down while `ω_y`,
`ω_z` and the horizontal position stay at zero (checked as `omega_perp`, `drift`).

With `F_n = m g`, equal-sphere `r_eff = R/2`, and sphere inertia `I = (2/5) m R²`.

## What the SDS model does (contact.rs)

The rolling "displacement" integrates the rolling velocity, `δ̇ = ω`, and the couple is

    τ = −k_r·δ − γ_r·ω      capped at   |τ| ≤ τ_max = μ_r·F_n·r_eff .

This is a spring (`k_r`, `rolling_stiffness`) + dashpot (`γ_r`, `rolling_damping`) + slider
(`μ_r`, `rolling_friction`). The benchmark validates the two regimes separately.

### 1. Elastic regime (Coulomb cap disengaged, large μ_r)

Below the cap the spin obeys the **exact damped linear oscillator**

    I·δ̈ + γ_r·δ̇ + k_r·δ = 0 ,   ω = δ̇ ,   δ(0)=0, ω(0)=ω0 .

`sweep.py` compares `ω(t)` point-by-point to the closed-form solution of that ODE:

- **over-damped** (`γ_r² ≫ 4 k_r I`): near-exponential decay set by the larger eigenvalue
  `|s₁| = (γ_r + √(γ_r²−4k_rI)) / 2I`. Match: **0.10 %** of ω0.
- **under-damped**: a decaying oscillation whose spring restoring force *reverses* the spin
  (ω goes negative). Match: **0.56 %** of ω0.

As a discriminating control the script also computes the **springless** (`k_r = 0`)
pure-dashpot curve `ω = ω0 e^{−(γ_r/I)t}`; it is off by **2.4 % / 131 %**, so the rolling
spring is genuinely exercised (the test is not passable with the spring dropped).

### 2. Coulomb-cap regime (slider saturated, small μ_r)

Under a large sustained spin (and `γ_r = 0` with a stiff spring), the spring+dashpot torque
exceeds the cap after a few steps, so the slider holds `τ = τ_max = μ_r F_n r_eff` and the
spin decays at the **exact constant rate**

    α = dω/dt = −τ_max/I = −(5/4)·μ_r·g/R

(the same closed form as `bench_rolling_decay`'s constant model and `bench_twisting_friction`).
`sweep.py` fits the saturated slope past the brief wind-up and compares to α across
μ_r ∈ {0.05, 0.10, 0.20}. Match: **0.00 %**.

## Running

```bash
# standalone (the over-damped elastic case)
cargo run --release --example bench_sds_rolling --no-default-features --features precision-double -- examples/bench_sds_rolling/config.toml

# full sweep + validation (build, run all cases, PASS/FAIL, plots)
$BENCH_PYTHON examples/bench_sds_rolling/sweep.py            # = generate + start + graph
```

The driver exits 0 on PASS. The initial spin ω0 is passed to the binary via the `SDS_OMEGA0`
env var (default 8.0 rad/s) so one binary serves every case without a source edit.

## Pass criteria

| check | tolerance | measured |
|---|---|---|
| elastic ω(t) vs exact damped oscillator | ≤ 1.5 % of ω0 | 0.10 % / 0.56 % |
| springless (k=0) control must be worse | ≥ 3× the SDS error | 2.4 % / 131 % |
| Coulomb-cap saturated slope vs (5/4)μ_r g/R | ≤ 3 % rel | 0.00 % |
| off-axis spin `omega_perp` | ≤ 0.1 % of ω0 | 0 |
| lateral drift | ≤ 10 µm | 0 |

## Reference

The SDS rolling model is the elastic–plastic spring-dashpot rolling resistance of
J. Ai, J.-F. Chen, J.M. Rotter, J.Y. Ooi, *"Assessment of rolling resistance models in
discrete element simulations"*, Powder Technology 206 (2011) 269–282, and is implemented in
LAMMPS `pair_granular` (rolling `sds`, `doc/src/pair_granular.rst`). Both the elastic
damped-oscillator dynamics and the `μ_r F_n r_eff` Coulomb cap are the model's own defining
behaviour — exactly what this benchmark validates.

## Outputs

- `data/roll_<case>.csv`, `data/sweep.csv` — per-case time series + validation summary (gitignored)
- `plots/sds_rolling_elastic.png` — elastic ω(t): DIRT vs exact vs springless
- `plots/sds_rolling_cap.png` — saturated spin-down vs α = (5/4)μ_r g/R
