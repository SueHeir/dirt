# bench_granular_conductivity — granular-temperature conductivity (vibro-fluidized bed)

Measures the one kinetic-theory transport coefficient the homogeneous LEBC
rheometer (`bench_lebc_shear`) can't see: the **conductivity of granular
fluctuation energy, κ**. It's the missing ingredient for a granular-temperature
SPH model's energy balance — and the rig doubles as a **de-fluidization** testbed,
the regime of interest for landing.

## Rig

A bed of glass beads under gravity rests on an **oscillating base wall** (periodic
in x,z). The vibrating base injects fluctuation (granular-temperature) energy at
the bottom; it conducts upward and is removed by inelastic dissipation, giving
steady **Φ(y)** and **T(y)** profiles — the canonical inhomogeneous KT benchmark,
and a gas-free fluidized bed. No mean shear, so the bulk energy balance is pure
conduction vs dissipation, which isolates κ.

## Run

```bash
cargo run --release --example bench_granular_conductivity --no-default-features -- examples/bench_granular_conductivity/config.toml
# or the driver (build + run + analyze):
python3 examples/bench_granular_conductivity/sweep.py
```

The recorder streams horizontal-slab profiles to `data/conductivity_profiles.csv`:
per y-bin the solid fraction Φ(y), the granular temperature
`T(y) = ⅓⟨|v − v̄_bin|²⟩` (per-bin streaming velocity removed, so the coherent
vibration isn't counted), and the kinetic fluctuation-energy flux
`q_y(y) = (1/V_bin) Σ ½m|v'|² v'_y`.

## κ extraction (printed + `plots/kappa_of_phi.png`)

In steady state with no mean shear, `dq_y/dy = −Γ(Φ,T)` and `q_y = −κ dT/dy`, with
`Γ = (12/d√π)(1−e²) ρ_s Φ² g₀ T^{3/2}`. Integrating dissipation from the top down
gives the upward flux `q_y(y)`, and `κ(y) = q_y(y)/(−dT/dy)`. Since Φ varies with
height, the single column sweeps **κ(Φ)**, compared to the KT (Lun/Gidaspow)
conductivity `κ* = κ/(ρ_s d √T)` (whose dilute limit gives κ*/η* = 15/4, as KT
requires). The **directly measured** kinetic heat flux `q_y` gives an independent
κ estimate (a lower bound, since it omits the collisional flux that dominates at
dense Φ) — plotted as open markers.

## De-fluidization (same rig)

Set the base `oscillate.amplitude = 0` (restart from a fluidized state): the bed
cools and consolidates — T(y,t) decays, Φ(y,t) rises, a consolidation front
descends. This is the landing-relevant transient; the homogeneous limit is the
`bench_*_haff_cooling` decay.

## Notes
- Frictionless, `e = 0.9` for a clean KT comparison; gravity in −y; base vibration
  `a(2πf)² ≈ 178 m/s²` (Γ_vib ≈ 18, well fluidized).
- κ measured here in a Φ-gradient setting is representative of how it's *used*
  (de-fluidization fronts have Φ gradients). The dense-Φ heat flux also carries a
  `∇Φ` contribution beyond `−κ∇T`; the energy-balance route captures the total.

## References
- Lun, Savage, Jeffrey & Chepurniy, *JFM* **140** (1984); Gidaspow, *Multiphase
  Flow and Fluidization* (1994) — KT conductivity.
- Companion: `bench_lebc_shear` (viscosity, dissipation, pressure); the two
  together calibrate the granular-temperature SPH (de-fluidization) model.
