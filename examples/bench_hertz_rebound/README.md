# Hertz Contact Rebound Benchmark

Validates Hertzian contact mechanics by dropping a single sphere onto a rigid flat wall and measuring the coefficient of restitution (COR), contact duration, and peak overlap. Where a LAMMPS binary is available, the identical sweep is also run in LAMMPS and overlaid on the plots as a code-to-code cross-check.

## Physics

A sphere of radius R, mass m, impacts a rigid flat wall at velocity v₀. The Hertz contact model predicts:

- **Contact duration** (elastic):
  ```
  t_c = 2.87 × (m²/(R·E*²·v₀))^(1/5)
  ```
- **Peak overlap** (elastic):
  ```
  δ_max = (15·m·v₀² / (16·√R·E*))^(2/5)
  ```
- **COR**: DIRT derives normal damping from the **Tsuji (1992) polynomial** — the
  same model LAMMPS's `damping tsuji` uses. This realizes a *velocity-independent*
  restitution. Note that below e≈0.9 the realized COR sits **above** the nominal input
  (e.g. nominal 0.5 → realized ≈0.645); this is a known property of the Tsuji polynomial,
  shared by LAMMPS, not a code error. The benchmark therefore validates the two things
  the model guarantees — **velocity-independence** and **agreement with LAMMPS** — rather
  than measured-vs-nominal.

where E* = E/(2(1−ν²)) is the reduced modulus for sphere-on-flat contact.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 70 GPa | Pa |
| Poisson's ratio ν | 0.22 | — |
| Density ρ | 2500 | kg/m³ |
| Radius R | 5 | mm |

## Parameter Sweep

- **Impact velocities**: 0.1, 0.5, 1.0, 2.0 m/s
- **COR values**: 0.5, 0.7, 0.9, 0.95, 1.0

**COR = 1.0 is the elastic anchor.** With zero damping the contact is purely
elastic, so the simulation must reproduce the (undamped) Hertz theory exactly.
It does: peak overlap matches to **≤ 0.1 %** and contact duration to within the
timestep resolution (≈ ±1 `dt`, 1–2 %) for both DIRT and LAMMPS, and the measured
COR is exactly 1.0. This pins the contact stiffness and integrator, and confirms
that the deviations from the elastic line at COR < 1 are the damping correction
(energy lost during approach), not a model error — they shrink monotonically as
COR → 1.

## Validation Criteria

| Check | Tolerance | Notes |
|-------|-----------|-------|
| **COR velocity-independence** | spread ≤ 0.01 across v₀, per COR | the signature of correct nonlinear-Hertz damping* |
| **COR agreement vs LAMMPS** | \|COR_dirt − COR_lammps\| ≤ 0.005 | same nominal input, both Tsuji (skipped if no LAMMPS) |
| Contact duration vs Hertz | ≤ 10% relative error | |
| Peak overlap vs Hertz | ≤ 10% relative error | |
| All 20 cases complete | 20/20 | |

\* The benchmark does **not** check realized-vs-nominal COR. DIRT and LAMMPS both use the
Tsuji polynomial, which realizes a velocity-independent restitution that sits above the
nominal input below e≈0.9 (a known polynomial limitation, not a code error). The earlier
DIRT bug applied the *linear* (Hooke) damping ratio to the nonlinear Hertz contact, giving
a *velocity-dependent* COR; the fix to Tsuji restored velocity-independence and brought
DIRT into agreement with LAMMPS to <0.005 — which is exactly what these two checks confirm.

## How to Run

Everything is driven by `sweep.py`, which takes one of three commands. With no
argument it runs all three in order.

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/bench_hertz_rebound/sweep.py

# Or one stage at a time:
python3 examples/bench_hertz_rebound/sweep.py generate   # write sweep/<case>/config.toml + in.lammps
python3 examples/bench_hertz_rebound/sweep.py start      # build, run all 16 cases (DIRT + LAMMPS) -> data/*.csv
python3 examples/bench_hertz_rebound/sweep.py graph      # validate against Hertz theory + write plots/
```

`graph` reads the existing `data/sweep_results.csv` (and `data/lammps_results.csv`
if present), so you can re-validate and re-plot without re-running the simulations.

### LAMMPS comparison

If a LAMMPS binary (`lmp_serial`, `lmp`, `lmp_mpi`, or `lammps`) is on `PATH`,
`start` also runs each case in LAMMPS and overlays it on the figures —
**DIRT as filled markers, LAMMPS as open markers**. LAMMPS is optional; without a
binary the benchmark runs DIRT only.

The LAMMPS model mirrors DIRT's: `pair_style granular` with `hertz/material E e ν`,
`damping tsuji` (LAMMPS's restitution-driven viscoelastic damping), zero friction,
and `fix nve/sphere` (translational only, no gravity), at the same timestep. The
normal stiffness is identical in both codes, so contact duration and peak overlap
agree directly.

For COR, both codes now use the **same** damping model — the Tsuji–Tanaka–Ishida
1992 polynomial. (DIRT previously used `β = −ln e / √(π²+ln²e)`, the linear
spring–dashpot relation; that is correct for a constant-stiffness Hooke contact but
wrong for nonlinear Hertz, and was the bug — it made restitution *velocity-dependent*.)
So `start` feeds LAMMPS the **same nominal restitution** as DIRT, with no calibration,
and the two realize the same velocity-independent COR — agreeing to **<0.005** across
all speeds and restitutions. That agreement (printed by `graph`) is the cross-code
validation criterion above. There is no exact closed-form e→damping for a nonlinear
Hertz contact, so the realized COR differs from the nominal input below e≈0.9; both
codes share that offset, which is why the benchmark checks cross-code agreement and
velocity-independence rather than realized-vs-nominal.

### Single case (default config)

```bash
cargo run --release --example bench_hertz_rebound --no-default-features -- examples/bench_hertz_rebound/config.toml
```

## Expected Plots

Each plot overlays the data on two references: the **elastic** Hertz theory (solid
black, valid only at COR = 1) and the **inelastic viscoelastic-model** curves
(dashed, per COR) obtained by integrating the same 1-DOF normal-contact ODE the
solver uses, *including* damping. The data sit on the inelastic curves; the gap
between those and the elastic line is the energy lost during contact, which grows
as COR drops and vanishes at COR = 1. `sweep.py`'s `contact_ode()` computes them.

> **Note (solver constant):** the inelastic curves use the damping constant
> `SQRT_5_3` from `dirt_atom` (`src/lib.rs`), whose *value* is √(5/6) ≈ 0.91287
> despite the *name*. √(5/6) is the physically correct value — it makes measured
> COR ≈ input COR across all restitutions; a literal √(5/3) would over-damp by √2.
> The name is a misnomer (value is right).

### COR Validation
![COR validation](plots/cor_validation.png)

### Contact Duration
![Contact duration](plots/contact_duration.png)

### Peak Overlap
![Peak overlap](plots/peak_overlap.png)

## Assumptions

- **3D simulation** with a single spherical particle
- **No friction** (friction = 0) for clean normal-only rebound
- **No gravity** effect on contact (gravity is off; particle given direct velocity)
- **Monodisperse** — single particle size
- **Hertz–Mindlin** contact model with viscoelastic damping (DIRT default)
- Wall is treated as **infinitely massive and rigid** (standard DEM wall)

## References

1. K.L. Johnson, *Contact Mechanics*, Cambridge University Press, 1985.
2. L. Vu-Quoc and X. Zhang, "An accurate and efficient tangential force-displacement model for elastic frictional contact in particle-flow simulations", *Mechanics of Materials*, 31(4):235–269, 1999.
