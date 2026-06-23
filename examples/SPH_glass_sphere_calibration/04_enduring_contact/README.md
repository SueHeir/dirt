# sphcal_enduring_contact — enduring-contact (rate-independent) stress branch

Tier-4 of the SPH glass-bead calibration. A triperiodic box of glass beads
(gravity off) is sheared at a constant rate γ̇ using DIRT's **native Lees–Edwards
simple shear** — a genuine triclinic box driven by the `[deform] xy` style, so the
flow is homogeneous and parallelizes across MPI ranks. This is the **frictional
dense** branch: across a Φ range, in steady state, the recorder reports the full
stress tensor, the granular temperature `T`, and the solid fraction `Φ`. From
these the sweep extracts the **enduring-contact stress branch** σ_contact(Φ) — the
rate-independent contact-network stress the granular-temperature de-fluidization
model consumes.

## Physics

The total DEM pressure splits into a collisional (rate-dependent, kinetic-theory)
part and a rate-independent enduring-contact part. Collisional kinetic theory
(Lun et al. 1984) predicts the pressure from the **measured** granular temperature
and solid fraction:

```
p_KT(Φ, T) = ρ_s · T · p*(Φ, e),    p*(Φ, e) = Φ [ 1 + 2(1+e) Φ g₀(Φ) ],
g₀(Φ) = (2 − Φ) / (2 (1 − Φ)³)      (Carnahan–Starling).
```

The **enduring-contact stress branch** is the dense-regime residual — what kinetic
theory cannot explain given the observed agitation:

```
σ_contact(Φ) ≈ p_DEM(Φ) − p_KT(Φ, T_measured).
```

It is ≈0 below Φ≈0.4 (collisions carry the stress; KT is accurate) and grows
toward jamming as a persistent contact network forms (this gap was seen opening at
Φ≈0.43 in `bench_lebc_shear`). The stress tensor is the Love–Weber contact virial
(`VirialStress`/V, normal + tangential) **plus** the kinetic term `Σ m v'⊗v' / V`
with the streaming velocity `v̄(y)=γ̇·y` removed; `T = ⅓⟨|v−v̄(y)|²⟩` likewise.

## Material Properties

Canonical glass-bead material (`[[dem.materials]]`):

| Property | Value |
|---|---|
| `youngs_mod` | 7.0e7 Pa (softened from ~65 GPa glass; rigid-grain limit) |
| `poisson_ratio` | 0.245 |
| `restitution` | 0.926 (measured glass–glass COR) |
| `friction` (μ_p) | 0.16 (measured glass–glass sliding) |
| `density` (ρ_s) | 2500 kg/m³ |
| diameter `d` | ≈ 0.5 mm, uniform radius ∈ [0.225, 0.275] mm |

## Parameter Sweep

A single **frictional** family. Each case compresses a loose insert (Φ≈0.30) to a
target Φ in the fixed 0.006 m box, then shears at γ̇ = 100 s⁻¹ to total strain ≈ 30.
The Φ grid brackets the onset of the enduring-contact branch:

```
Φ ∈ {0.30, 0.40, 0.45, 0.50, 0.55, 0.58, 0.60}
```

## Validation Criteria

`graph` prints a PASS/FAIL table and exits non-zero on FAIL:

1. **Sane stress** — all `p_DEM > 0`.
2. **Dilute floor** — at the lowest Φ, `σ_contact/p_DEM < 0.20` (collisions
   dominate; KT explains the pressure, residual ≈ 0).
3. **Branch growth** — at the densest Φ, `σ_contact` clearly exceeds the dilute
   floor and `σ_contact/p_DEM > 0.20` (the enduring-contact branch has opened).

Steady-state drift (relative change of `p` across the averaging window) is
reported and warned above 20%, but is not a hard failure: the densest near-jamming
cases legitimately plateau slowly. See `plots/convergence.png`.

## How to Run

Single representative case:

```bash
cargo run --release --example sphcal_enduring_contact --no-default-features -- examples/SPH_glass_sphere_calibration/04_enduring_contact/config.toml
```

Full sweep (generate configs → run → validate + plot):

```bash
python3 examples/SPH_glass_sphere_calibration/04_enduring_contact/sweep.py            # all three
python3 examples/SPH_glass_sphere_calibration/04_enduring_contact/sweep.py generate
python3 examples/SPH_glass_sphere_calibration/04_enduring_contact/sweep.py start
python3 examples/SPH_glass_sphere_calibration/04_enduring_contact/sweep.py graph
```

## Expected Plots

| Path | Contents |
|---|---|
| `plots/sigma_contact.png` | the deliverable: σ_contact(Φ) = p_DEM − p_KT, ≈0 at low Φ, rising toward jamming |
| `plots/stress_decomposition.png` | p_DEM (measured) vs collisional p_KT(Φ,T) vs the residual σ_contact, all vs Φ |
| `plots/convergence.png` | p vs strain per case with the averaging window shaded (steady-state check) |

`data/sigma_contact.csv` (gitignored) tabulates `phi, T, p_dem, p_kt,
sigma_contact, frac_contact`.

## Assumptions

- **Gravity is off** and the box is fixed (fixed-Φ sweep) — the simplest valid
  route to a Φ-parameterized residual; pressure control is a later refinement.
- p_KT uses the standard Lun (Carnahan–Starling g₀) EOS evaluated at the
  **measured** T and Φ. Near jamming, KT's own assumptions (instantaneous binary
  collisions, molecular chaos) break down, which is precisely why the residual is
  the *definition* of the dense, enduring-contact regime rather than an error bar.
- **Timestep** `dt ≈ 2e-7 s` is set from the Rayleigh criterion for `E = 7e7`,
  `d = 0.5 mm`; re-derive if `E` or `d` changes.
- Spheres only; mild polydispersity suppresses crystallization.

## No LAMMPS overlay

The deliverable is a *residual against a specific KT closure* evaluated at the
measured T, not a raw stress that maps onto a LAMMPS observable. LAMMPS reproduces
the underlying frictional shear (as in `bench_lebc_shear`), but the σ_contact
decomposition is a post-processing definition tied to the Lun EOS here, so no
cross-code overlay is added. The raw frictional stress vs Φ is cross-checked
against LAMMPS in `bench_lebc_shear`.

## References

- Lun, Savage, Jeffrey & Chepurniy, *JFM* **140** (1984) — kinetic-theory pressure EOS.
- GDR MiDi, *Eur. Phys. J. E* **14** (2004); da Cruz et al., *PRE* **72** (2005) — frictional dense shear.
- Companion: `examples/bench_lebc_shear` (the basis recorder + KT closure).
