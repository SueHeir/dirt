# SPH glass-sphere calibration · 01 — shear rheology, μ(I)

Tier-1 of the SPH glass-sphere calibration campaign. A triperiodic box of glass
beads (gravity off) is sheared at a constant rate γ̇ using DIRT's **native
Lees–Edwards simple shear** — a genuine triclinic box driven by the `[deform] xy`
style, with the box tilt and streaming-velocity remap handled in fractional
(lamda) coordinates so the flow is homogeneous and parallelizes across MPI ranks.
In steady state the recorder reports the full stress tensor, the granular
temperature, and the solid fraction; from these the sweep forms the inertial
number `I`, the effective friction `μ(I)`, and `Φ(I)` — the closure the MUD SPH
solver consumes, together with the critical density `ρ_c = Φ_max·ρ_s`.

This is the **frictional/production** path (canonical glass `μ_p = 0.16`). The
frictionless kinetic-theory cross-check (a separate sub-sweep validated against
Lun et al. 1984) lives in `examples/bench_lebc_shear`.

## Physics

Driving Lees–Edwards shear: `[deform] xy = { style = "erate", rate = γ̇ }`
(flow = x, gradient = y, vorticity = z), enabled only in the `shear` stage.

The stress tensor combines the contact (Love–Weber) virial and the kinetic term
with the streaming velocity `v̄(y)=γ̇·y` removed:

- `σ_ij = (VirialStress_ij)/V + Σ m v'_i v'_j / V`, `v' = v − v̄(y) x̂`
- pressure `P = ⅓ tr σ`; shear stress `σ_xy`; normal-stress differences
  `N₁ = σ_xx − σ_yy`, `N₂ = σ_yy − σ_zz`
- granular temperature `T = ⅓⟨|v − v̄(y)|²⟩`; solid fraction `Φ`

The dimensionless rheology closure (GDR MiDi 2004; da Cruz et al. 2005):

- inertial number `I = γ̇ d √(ρ_s / P)`
- effective friction `μ = |σ_xy| / P`
- the fit `μ(I) = μ_s + (μ₂ − μ_s)/(I₀/I + 1)`, the dilatancy law `Φ(I)`,
  and `ρ_c = Φ_max·ρ_s` from the densest (quasi-static) point.

## Material Properties

Canonical glass beads (set identically in `config.toml` and `sweep.py`):

| property | value |
|---|---|
| Young's modulus `E` | 7.0e7 Pa (softened from ~65 GPa, rigid-grain limit) |
| Poisson ratio `ν` | 0.245 |
| restitution `e` | 0.926 |
| sliding friction `μ_p` | 0.16 |
| density `ρ_s` | 2500 kg/m³ |
| mean diameter `d` | 0.5 mm (uniform radius 0.225–0.275 mm) |

## Parameter Sweep

Each case inserts a **loose** pack (random insertion saturates near Φ≈0.38) and
**compresses** it to the target box (→ target Φ) with a velocity-style `[deform]`
ramp before shearing, so the dense end of the Φ range is reachable. Shear
duration is set by total **strain** (`TARGET_STRAIN`), not a fixed step count, so
the low-γ̇ (low-I) cases run longer.

- target solid fraction `Φ ∈ {0.40, 0.50, 0.55, 0.58, 0.60}`
- shear rate `γ̇ ∈ {10, 30, 100, 300} s⁻¹`

The grid spans quasi-static to inertial flow, populating the μ(I)/Φ(I) curve.

## Validation Criteria

`graph` prints a table and returns **PASS/FAIL** (non-zero exit on FAIL):

1. **Sanity** — every steady case has `P > 0`, finite `σ_xy`, and `0 < Φ < 1`.
2. **Steady state** — `drift` (change of `P` across the averaging window, last 50%
   of strain) is flagged when ≥ 15%; `plots/convergence.png` shows stress and T
   vs strain with that window shaded.
3. **Closure ordering** — the fit is physical: `μ_s < μ₂` and `I₀ > 0`.
4. **Bagnold** — `P/γ̇²` spread per Φ reports inertial-regime consistency.

The fitted `(μ_s, μ₂, I₀, Φ_max, ρ_c)` are written to `data/calibration.yaml`
for `mud_constitutive::MaterialParams`. Glass-bead anchors: `μ_s ≈ 0.38`,
`μ₂ ≈ 0.64`, `I₀ ≈ 0.28` (possibly lower for spheres).

## How to Run

Single representative case:

```bash
cargo build --release --no-default-features --example sphcal_shear_rheology
./target/release/examples/sphcal_shear_rheology examples/SPH_glass_sphere_calibration/01_shear_rheology/config.toml
```

Full sweep (generate configs → run → validate + fit + plot):

```bash
python3 examples/SPH_glass_sphere_calibration/01_shear_rheology/sweep.py            # all three
# or individually:
python3 examples/SPH_glass_sphere_calibration/01_shear_rheology/sweep.py generate
python3 examples/SPH_glass_sphere_calibration/01_shear_rheology/sweep.py start
python3 examples/SPH_glass_sphere_calibration/01_shear_rheology/sweep.py graph
```

## Expected Plots

| Path | Contents |
|---|---|
| `plots/mu_of_I.png` | DEM `μ` vs `I` with the fitted `μ(I) = μ_s + (μ₂−μ_s)/(I₀/I+1)` |
| `plots/phi_of_I.png` | dilatancy law `Φ(I)` with the `Φ_max` line |
| `plots/convergence.png` | `σ_xy` and `T` vs strain per case, averaging window shaded |

## Assumptions

- **Gravity is off** and the box is fixed (fixed-Φ sweep) — the simplest valid
  route to μ(I)/Φ(I); convert Φ→P post hoc via the measured `σ_yy`. Pressure
  control (a servo on the gradient direction) is a later refinement.
- **Timestep** `dt ≈ 2e-7 s` is set from the Rayleigh criterion for `E = 7e7`,
  `d = 0.5 mm`; re-derive if `E` or `d` changes.
- **No LAMMPS overlay here.** The physics is reproducible in LAMMPS
  (`fix deform xy erate … remap v`), but the kinetic-theory / cross-code
  validation leg lives in `examples/bench_lebc_shear`; this example is the
  production μ(I) extraction for the SPH closure and validates against the
  empirical μ(I) law plus internal steadiness/sanity checks.

## References

- GDR MiDi, *Eur. Phys. J. E* **14** (2004) — μ(I) from DEM simple shear.
- da Cruz, Emam, Prochnow, Roux & Chevoir, *PRE* **72** (2005) — inertial-number scaling.
- Lun, Savage, Jeffrey & Chepurniy, *JFM* **140** (1984) — kinetic-theory transport coefficients (frictionless cross-check in `bench_lebc_shear`).
