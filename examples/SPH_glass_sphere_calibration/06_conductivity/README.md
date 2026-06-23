# sphcal_conductivity — granular-temperature conductivity κ(Φ)

Calibrates the **conductivity closure κ(Φ)** for the granular-temperature SPH
de-fluidization energy balance, using the canonical glass-bead material. This is
the one kinetic-theory transport coefficient a homogeneous (LEBC) rheometer can't
see — the conduction of granular fluctuation energy that sets how fast a
de-fluidizing bed cools and consolidates.

## Physics

A bed of glass beads under gravity rests on an **oscillating base wall** (periodic
in x,z). The vibrating base injects fluctuation (granular-temperature) energy at
the bottom; it conducts upward and is removed by inelastic + frictional
dissipation, giving steady **Φ(y)** and **T(y)** profiles. With no mean shear the
bulk energy balance is pure conduction vs dissipation, which isolates κ:

    d q_y/dy = −Γ(Φ,T),   q_y = −κ dT/dy,   Γ = (12/d√π)(1−e²) ρ_s Φ² g₀ T^{3/2}.

Integrating dissipation from the top down gives the upward energy flux q_y(y), and
**κ(y) = q_y(y) / (−dT/dy)**. Since Φ varies with height, the single column sweeps
**κ(Φ)**, compared against the KT (Lun/Gidaspow) conductivity
`κ* = κ/(ρ_s d √T)`. The directly recorded kinetic heat flux gives an independent
κ estimate (a lower bound; it omits the collisional flux that dominates at dense Φ).

**Why friction matters here:** with friction (μ_p = 0.16) the extra dissipation
steepens the T(y) gradient enough to measure κ. A frictionless near-elastic bed
stays essentially isothermal, so the canonical material's friction is load-bearing
for this calibration.

## Material Properties

Canonical glass bead (hertz contact):

| property      | value   |
|---------------|---------|
| youngs_mod    | 7.0e7   |
| poisson_ratio | 0.245   |
| restitution   | 0.926   |
| friction      | 0.16    |
| density       | 2500.0  |

Beads: uniform radius 0.225–0.275 mm, count 3000.

## Parameter Sweep

A single column run. The sweep over **Φ** is *spatial*: the steady bed spans dilute
(top) to dense (base) solid fractions, so one run yields κ(Φ) across the Φ range.
The base drive is fixed at amplitude 0.00025 m, frequency 160 Hz
(a(2πf)² ≈ 253 m/s², Γ_vib ≈ 26, well fluidized).

## Validation Criteria

The conduction–dissipation route assumes all dissipation is the KT collisional
term Γ(Φ,T); frictional dissipation (which is what makes the gradient measurable)
is not in that closure, so the extracted κ* is an **effective** conductivity that
lumps the friction-steepened gradient into the closure — exactly the quantity the
SPH energy balance consumes. The KT (Lun/Gidaspow) curve is therefore a **reference
overlay**, not a strict PASS/FAIL target. Sanity checks: κ* > 0, monotone-ish in Φ,
within ~1 decade of the KT reference over the resolved Φ band.

## How to Run

```bash
# single-case standalone run:
cargo run --release --example sphcal_conductivity --no-default-features -- examples/SPH_glass_sphere_calibration/06_conductivity/config.toml

# or the driver (build + run + analyze):
python3 examples/SPH_glass_sphere_calibration/06_conductivity/sweep.py start
python3 examples/SPH_glass_sphere_calibration/06_conductivity/sweep.py graph
python3 examples/SPH_glass_sphere_calibration/06_conductivity/sweep.py          # both
```

The recorder streams horizontal-slab profiles to `data/conductivity_profiles.csv`:
per y-bin Φ(y), the granular temperature `T(y) = ⅓⟨|v − v̄_bin|²⟩` (per-bin
streaming velocity removed, so the coherent vibration isn't counted), and the
kinetic fluctuation-energy flux `q_y(y) = (1/V_bin) Σ ½m|v'|² v'_y`.

## Expected Plots

- `plots/profiles.png` — steady Φ(y), T(y), q_y(y): a settled bed at the base, T
  highest near the vibrating base and decaying upward, finite upward heat flux.
- `plots/kappa_of_phi.png` — DEM κ*(Φ) (energy-balance total, and kinetic-flux-only
  lower bound) vs the KT (Lun/Gidaspow) reference curve.

## Assumptions

- No interstitial gas (gas-free granular bed); gravity in −y.
- Steady state: profiles are time-averaged over the last half of snapshots.
- κ extracted in a Φ-gradient setting is representative of how it's *used*
  (de-fluidization fronts have Φ gradients). The dense-Φ heat flux also carries a
  ∇Φ contribution beyond −κ∇T; the energy-balance route captures the total.
- LAMMPS overlay omitted: the deliverable is the κ(Φ) *closure* extracted from the
  energy balance, and the KT reference is the meaningful comparison; a LAMMPS bed
  would reproduce the profiles but not change the extracted closure.

## References

- Lun, Savage, Jeffrey & Chepurniy, *JFM* **140** (1984); Gidaspow, *Multiphase
  Flow and Fluidization* (1994) — KT conductivity.
- Companion calibration stages under `SPH_glass_sphere_calibration/` together
  calibrate the granular-temperature SPH (de-fluidization) model.
