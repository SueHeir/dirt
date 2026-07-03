# bench_twisting_friction — twisting (torsional) friction spin-down

Validates DIRT's **twisting-friction couple** — for both the `constant` and the
`sds` (spring-dashpot-slider) twisting models — against the *exact* analytical
spin-down of a sphere twisting on an enduring contact.

Twisting friction resists relative rotation about the **contact normal** n̂ (the
torsional degree of freedom). It is the last of the contact friction trio
(sliding, rolling, twisting) and, unlike sliding and rolling, was previously
exercised by no benchmark — this closes that gap on the cleanest possible path.

## Setup

Two identical spheres (R = 5 mm) are stacked along +z:

- **Lower sphere** — a frozen anchor (`[[freeze]]`, so it can neither translate
  nor spin: a true immovable contact partner).
- **Upper sphere** — rests on the anchor under gravity, seated at the **static
  Hertz overlap** so the normal contact force is `F_n = m g` from step 0 (no
  settling ring-down), and is given a **pure twisting spin** `ω = (0, 0, ω₀)`
  about the contact normal n̂ = +ẑ.

A pure twist about n̂ produces **zero** relative surface velocity at the contact
point (it lies on the spin axis): `v_rel = ω × (∓R n̂) = 0`. Hence there is no
sliding and no rolling — the **only** contact torque is the twisting-friction
couple, applied purely along n̂. The spin therefore decays about a fixed axis
while ω_x, ω_y and the horizontal position stay identically zero (checked).

The upper sphere's seating and spin (there is no `omega` knob in
`[[particles.insert]]`) are applied once at the first step in `main.rs`.

## Theory (exact)

DIRT's twisting couple opposes the twist with magnitude

```
τ_tw = μ_tw · F_n · r_eff
```

For the resting sphere, vertical equilibrium gives `F_n = m g`; for two equal
spheres the effective contact radius is `r_eff = R₁R₂/(R₁+R₂) = R/2`; and a
sphere has `I = (2/5) m R²`. The twist rate then obeys

```
dω/dt = − τ_tw / I = − (μ_tw · m g · R/2) / ((2/5) m R²)
      = − (5/4) · μ_tw · g / R                          ← exact, constant
```

The mass cancels: the spin-down rate depends only on `μ_tw`, `g`, and `R`.

- **`constant` model** — the couple is applied at full magnitude every step, so
  ω_z(t) decays linearly at the rate above the entire way down.
- **`sds` model** — a torsional spring winds up to the same cap
  `τ_max = μ_tw F_n r_eff` within a few steps, then the slider holds it (the
  accumulated twist keeps growing, so it stays saturated), giving the identical
  constant deceleration past a brief elastic wind-up. `sweep.py` fits the slope
  past that transient, so both models validate against the same α.

## What `sweep.py` checks (PASS/FAIL gated)

For each model ∈ {`constant`, `sds`} and each `μ_tw` ∈ {0.05, 0.10, 0.20}:

1. **Spin-down rate** — linear-fit `ω_z(t)` over `0.15 ω₀ … 0.85 ω₀` (above the
   near-zero discretization chatter, below the SDS wind-up) and require the
   fitted deceleration to match `α = (5/4) μ_tw g / R` within **3 %**.
2. **Twist purity** — the off-axis spin `|ω_perp|` must stay < 0.1 % of ω₀ and
   the lateral drift < 10 µm (confirming the contact drives *only* the torsional
   DOF, with no spurious sliding, rolling, or translation).

Exit code is 0 only if every case passes (the runner treats exit 0 = PASS).

## Running

```bash
source ~/projects/.build-env
# one-shot: generate configs, build, run all cases, validate + plot
"$BENCH_PYTHON" examples/bench_twisting_friction/sweep.py
# or the regression runner (records PASS/FAIL):
~/projects/automation/bin/run-bench.sh examples/bench_twisting_friction
```

Single standalone case (default constant model, μ_tw = 0.1):

```bash
cargo run --release --example bench_twisting_friction \
  --no-default-features --features precision-double -- \
  examples/bench_twisting_friction/config.toml
```

## Latest result

```
      model  mu_tw     a_fit    a_pred  rel_err  max_perp   drift
   constant  0.050   122.625   122.625    0.00%  0.00e+00  0.00e+00
   constant  0.100   245.250   245.250    0.00%  0.00e+00  0.00e+00
   constant  0.200   490.500   490.500    0.00%  0.00e+00  0.00e+00
        sds  0.050   122.625   122.625    0.00%  0.00e+00  0.00e+00
        sds  0.100   245.250   245.250    0.00%  0.00e+00  0.00e+00
        sds  0.200   490.500   490.500    0.00%  0.00e+00  0.00e+00
RESULT: PASS
```

Because `F_n = m g` holds exactly at the static-overlap seating and the twist
stays pure, the spin-down ODE is `dω/dt = const`, which the leapfrog integrator
reproduces to round-off — so the fitted rate matches theory to < 1e-6 relative
(the 3 % tolerance is comfortable headroom, not a fudge).

![twist spin-down](plots/twist_spindown.png)
![spin-down vs mu_tw](plots/spindown_vs_mu_tw.png)

## Reference

The constant-directional-torque family of rolling/twisting resistance models:
J. Ai, J.-F. Chen, J.M. Rotter, J.Y. Ooi, *"Assessment of rolling resistance
models in discrete element simulations"*, Powder Technology **206** (2011)
269–282 (model A, "constant directional torque"), here applied to the twisting
(normal-axis) degree of freedom.
