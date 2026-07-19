# bench_granular_conductivity — granular-temperature conductivity (vibro-fluidized bed)

Measures the one kinetic-theory transport coefficient the homogeneous LEBC
rheometer (`bench_lebc_shear`) can't see: the **conductivity of granular
fluctuation energy, κ**. The rig is a DEM transport benchmark and also provides
a controlled cooling-and-consolidation testbed.

## Rig

A bed of glass beads under gravity rests on an **oscillating base wall** (periodic
in x,z). The vibrating base injects fluctuation (granular-temperature) energy at
the bottom; it conducts upward and is removed by inelastic dissipation, giving
steady **Φ(y)** and **T(y)** profiles — the canonical inhomogeneous KT benchmark,
and a gas-free fluidized bed. No mean shear, so the bulk energy balance is pure
conduction vs dissipation, which isolates κ.

## Run

```bash
python3 examples/bench_granular_conductivity/sweep.py         # BOUNDED smoke gate (PASS/FAIL) — the harness default
python3 examples/bench_granular_conductivity/sweep.py smoke   # same as no-arg
python3 examples/bench_granular_conductivity/sweep.py full    # full scientific run (config.toml): profiles + κ(Φ)
# a single full run directly:
cargo run --release --example bench_granular_conductivity --no-default-features --features precision-double -- examples/bench_granular_conductivity/config.toml
```

### Bounded smoke gate (CI)

The full run (3000 grains, 45 mm column, 4M steps at dt=1e-7 ≈ 0.4 s of drive)
legitimately overran the 1800 s automation cap every hourly run — so it validated
nothing. `sweep.py` with **no argument** now runs a **bounded gate** on the
shallower [`config.smoke.toml`](config.smoke.toml) (800 grains, 30 mm column,
800k steps at dt=2e-7 ≈ 0.16 s): it reaches a steady fluidized state and asserts
the measured Fourier-law conductivity `κ*(Φ)` is positive, finite, and within
order-unity of kinetic theory (median `κ*_EB/κ*_KT` ∈ [0.4, 6.0]). It completes
in ~50 s and prints `ALL CHECKS PASSED`/`CHECKS FAILED`. The full κ*(Φ)-vs-KT
scientific run (`config.toml`) is **unchanged** and still run via `sweep.py full`.

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

![Dimensionless conductivity vs kinetic theory](plots/kappa_of_phi.png)

*Measured dimensionless conductivity `κ*(Φ)` from the DEM column against the
Lun/Gidaspow kinetic-theory reference. The shaded band is the smoke gate's
order-unity acceptance window (0.4-6x KT); `sweep.py smoke` prints PASS/FAIL.*

![Steady vibro-fluidized profiles](plots/profiles.png)

*Steady solid fraction, granular temperature, and kinetic heat-flux profiles used
to extract `κ*(Φ)` from the conduction-dissipation balance.*

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
- Companion: `bench_lebc_shear` covers viscosity, dissipation, and pressure in
  homogeneous shear; this rig isolates the inhomogeneous conductivity term.
