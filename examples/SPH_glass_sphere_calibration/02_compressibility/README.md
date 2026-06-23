# Isotropic compression — bulk modulus K(Φ)

Deliverable **#2** of the SPH glass-sphere calibration campaign. It produces the
**bulk compressibility closure** the MUD SPH continuum model consumes to close the
pressure–density relation: the modulus `K(Φ)` and the equation of state `P(Φ)`,
measured from a quasi-static isotropic compression of a triperiodic pack of glass
beads.

## Physics

A fully periodic box of glass beads (gravity off) is compressed quasi-statically
and **isotropically**: the `[deform]` velocity driver pushes all three box faces
inward at equal speed, so the solid fraction Φ rises while the pack stays
homogeneous (no preferred direction, no mean shear). The recorder reads the
Love–Weber contact virial `Res<VirialStress>` and `Res<Domain>` and streams each
thermo interval:

    p = −trace(VirialStress) / (3·V)              (mean pressure; + small kinetic term)
    Φ = Σ (π/6) dᵢ³ / V                            (solid fraction)

As the pack is squeezed through random close packing (`Φ_c ≈ 0.58–0.64`), enduring
contacts percolate and the pressure climbs steeply — the **equation of state**
`P(Φ)`. The bulk modulus is its log-slope:

    K(Φ) ≈ ΔP / (ΔΦ/Φ) = dP/d(lnΦ) = Φ · dP/dΦ.

`sweep.py` fits the jammed branch to the contact-mechanics EOS form
`P(Φ) = A·(Φ − Φ_c)^α` (Hertzian grains near jamming give `α ≈ 1.5`) and reports
`K(Φ)` at representative solid fractions.

## Material Properties

Canonical glass-sphere material (shared by the whole calibration campaign):

| property        | value    | note                                         |
|-----------------|----------|----------------------------------------------|
| youngs_mod      | 7.0e7 Pa | softened from ~65 GPa glass (rigid-grain limit) |
| poisson_ratio   | 0.245    | glass                                        |
| restitution     | 0.926    | measured glass–glass COR                     |
| friction        | 0.16     | measured glass–glass sliding μ_p             |
| density         | 2500 kg/m³ | set on the particle insert (config schema) |
| diameter        | 0.45–0.55 mm uniform | mild polydispersity, d̄ = 0.5 mm |

## Parameter Sweep

A single representative case: one loose pack (Φ ≈ 0.30) compressed continuously up
the dense branch (Φ ≈ 0.63). The *sweep* over Φ is produced **within the single
run** — every thermo sample is one (Φ, p) point on the EOS — rather than as
separate fixed-volume cases. Contact model: Hertz.

## Validation Criteria

`sweep.py graph` runs `validate()` and exits non-zero on FAIL:

- **Monotonic EOS**: on the jammed branch (`p ≥ 1 kPa`), `p(Φ)` rises monotonically
  in ≥ 80 % of consecutive samples.
- **EOS fit quality**: `P = A·(Φ − Φ_c)^α` fits the jammed branch with `R² ≥ 0.95`
  and a physically sensible exponent `1 ≤ α ≤ 3`.

The extracted `K(Φ)` at Φ = 0.60, 0.62, 0.63 is the deliverable emitted into the
shared `calibration.yaml`.

## How to Run

```bash
# all three stages
python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py

# or individually
python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py generate
python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py start
python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py graph

# single standalone case
cargo run --release --example sphcal_compressibility --no-default-features -- \
    examples/SPH_glass_sphere_calibration/02_compressibility/config.toml
```

## Expected Plots

- `plots/eos_p_of_phi.png` — pressure vs solid fraction (log p), with the jammed
  branch highlighted, the fitted EOS `P = A(Φ − Φ_c)^α`, and `Φ_c`.
- `plots/bulk_modulus.png` — the bulk modulus `K(Φ) = dP/d(lnΦ)` from the fitted
  EOS, marked at the reported solid fractions.

## Assumptions

- **Quasi-static**: light viscous damping drains the work the moving faces inject,
  so the measured pressure tracks the rate-independent EOS branch rather than a
  dynamic (Bagnold) one.
- **Soft-grain modulus**: `youngs_mod` is softened to the rigid-grain limit, so the
  absolute scale of `K` reflects the simulated contact stiffness, not real glass;
  the EOS *shape* (Φ_c, α) is the transferable closure.
- **Mild polydispersity** suppresses crystallization, keeping the random close
  packing (`Φ_c`) physical.
- **No LAMMPS overlay**: LAMMPS's `granular` package has no apples-to-apples
  isotropic box-deform-while-measuring-virial EOS rig; the closure is checked
  against the contact-mechanics EOS form instead of a cross-code overlay.

## References

- O'Hern, Silbert, Liu, Nagel, *Jamming at zero temperature and zero applied
  stress: the epitome of disorder*, Phys. Rev. E 68, 011306 (2003) — the
  `P ∝ (Φ − Φ_c)^α` jamming EOS for soft repulsive grains.
- Makse, Gland, Johnson, Schwartz, *Granular packings: nonlinear elasticity, sound
  propagation, and collective relaxation dynamics*, Phys. Rev. E 70, 061302 (2004).
