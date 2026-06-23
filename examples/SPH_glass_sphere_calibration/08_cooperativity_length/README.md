# sphcal_cooperativity_length

DEM measurement of the **granular cooperativity length** ξ — the closure for the
MUD SPH model's *nonlocal* (NGF-style) branch.

## Why

MUD's contact branch is local μ(I): rigid below yield. Real dense beds **creep
below yield** when a neighbouring region flows (Henann–Kamrin's Nonlocal Granular
Fluidity). MUD keeps the one irreplaceable NGF idea inside its granular-temperature
model — a cooperativity length that **diverges at yield** and is *athermal* (it
survives as T → 0, where the kinetic-theory conduction length vanishes):

```
ξ(μ) = A · d / √(μ − μ_s)        # the closure this rig calibrates: amplitude A
```

and it drives the fluidity from temperature via the **Zhang–Kamrin bridge**
(PRL 2017, "the fluidity field is a measure of velocity fluctuations"):

```
g = γ̇/μ  ∝  √T                  # this rig also validates this relation
```

`MUD/crates/mud_constitutive` consumes `A` as `MaterialParams::coop_amplitude`
(0 = nonlocal off). This rig measures it.

## What it does

Same Lees–Edwards homogeneous-shear rig as `bench_lebc_shear` (so μ, p, T, Φ, I are
measured identically), plus a **spatial velocity-fluctuation correlation** along the
vorticity (z) direction:

```
C(Δz) = ⟨δv'(0)·δv'(Δz)⟩ / ⟨|δv'|²⟩ ,   δv' = v − γ̇(y − y_c) x̂
```

restricted to near-columns (|Δx|, |Δy| < d) so the Lees–Edwards xy tilt is
irrelevant. ξ is the integral of C(Δz) to its first zero crossing. As μ → μ_s
(slower shear, denser/closer to jamming) the fluctuations become long-range and ξ
grows — that growth, fit to `A d/√(μ−μ_s)`, gives `A`.

## Run

Single representative case (γ̇ = 50 s⁻¹):
```bash
cargo run --release --example sphcal_cooperativity_length --no-default-features -- \
    examples/SPH_glass_sphere_calibration/08_cooperativity_length/config.toml
```
Full closure sweep + fit (generates a γ̇ grid, runs each, fits A and g∝√T):
```bash
python3 examples/SPH_glass_sphere_calibration/08_cooperativity_length/sweep.py --run
```
Outputs `data/.../cooperativity_results.csv` (per-sample μ, p, T, Φ, I, g, ξ, and the
correlation bins C0…C11), `data/calibration.yaml` (the `A` to drop into MUD), and
`plots/cooperativity.png`.

## Notes / limitations (v0)

- **Single-rank**: the pair correlation is computed on local atoms with minimum
  image; run serial (the configs are single-process).
- ξ is measured as an *integral* length of C(Δz); `sweep.py` can be extended to an
  exponential-fit length if the decay is cleaner.
- The closure is exactly the one MUD needs, but the *interpolation hypothesis*
  ℓ_nl² = ℓ_KT²(T) + ξ(μ)² (kinetic-theory conduction length crossing over to the
  athermal cooperativity length) should be checked by comparing ξ(μ→μ_s) against
  the KT conduction length at the same T — a follow-on analysis on the same CSVs.
