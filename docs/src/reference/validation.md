# Validation & Benchmarks

DIRT ships a suite of validation examples under `examples/`. Each couples a
small simulation to a reference — a closed-form result, an empirical
correlation, or a LAMMPS run of the same problem — and checks measured
quantities against it with explicit tolerances. Every example's `sweep.py graph`
step prints a PASS/FAIL verdict, so the suite is a regression net, not a gallery
of passing runs.

## The validation discipline

The benchmarks are written to be useful *and* honest. Each one states its
result, then says plainly where it is weak: an idealization, an empirical fit, a
check that is really self-consistent (a model returning its own input), or a
regime the simulation never reaches. They sort into three **evidence tiers**, in
decreasing strength:

- **Analytical** — agreement with a closed-form reference (Hertz contact
  duration, Euler–Bernoulli beam deflection, Haff's cooling law).
- **Cross-code** — agreement with LAMMPS. This tests implementation consistency
  under a *shared* contact model, not correctness against physical reality.
- **Empirical / law / qualitative** — agreement with a scaling exponent, a
  functional form, or a correlation with fitted constants (Beverloo discharge,
  angle of repose vs. friction).

There is no direct comparison to experimental data in the suite; references are
analytical, empirical, the experimentally-established Maw oblique-impact curve
(used as a theory curve), or LAMMPS. The authoritative, continuously-updated
write-up — what each figure shows and exactly where each test is weak — lives in
[`examples/VALIDATION.md`](https://github.com/SueHeir/dirt/blob/master/examples/VALIDATION.md).

## The benchmarks

| Example | What it checks | Tier |
|---|---|---|
| `bench_hertz_rebound` | COR, contact duration, peak overlap for a sphere on a flat wall | Analytical + cross-code |
| `bench_oblique_impact` | tangential rebound (Maw curve) | Analytical |
| `bench_sliding_friction` | slip-to-roll transition on a frictional wall | Analytical |
| `bench_rolling_decay` | rolling-resistance velocity decay | Analytical |
| `bench_sphere_haff_cooling` / `bench_rod_haff_cooling` / `bench_clump_haff_cooling` | Haff's `T_g ∝ t⁻²` cooling for spheres, rods, and clumps | Analytical |
| `bench_lebc_shear` | steady-shear rheology under Lees–Edwards boundaries | Analytical / law |
| `bench_angle_of_repose` | repose angle vs. friction coefficient | Empirical |
| `bench_column_collapse` | granular column runout scaling | Empirical |
| `bench_hopper_beverloo` | hopper discharge `W ∝ (D − kd)^{3/2}` | Empirical |
| `bench_jkr_adhesion` | JKR force–separation and pull-off | Analytical |
| `bench_fiber_crossover` | bonded-fiber response with explicit bonds | Analytical |
| `bench_plate_sinkage` | pressure–sinkage response | Empirical |
| `bench_granular_conductivity` | granular heat conduction | Law / qualitative |
| `perf_mpi_scaling` | strong/weak MPI scaling (performance, not physics) | — |

These tests catch real bugs. The oblique-impact validation alone drove two
contact-model fixes (a tangential damping-sign error injecting energy, and the
requirement that a frozen contact partner also have its rotation frozen), and the
rebound benchmark surfaced a mislabeled damping constant.

## Calibration note: nominal vs. realized restitution

The input `restitution` in `[[dem.materials]]` is the **target** coefficient of
restitution, and for a binary collision DIRT realizes it. `build_pair_tables`
inverts the exact Hertz `COR(β)` curve numerically (bisection over a
once-integrated, velocity-independent collision) to find the damping ratio `β`
that reproduces the requested `e`; the Hooke path uses the closed-form
`β = −ln(e)/√(π² + ln²e)`. See [Materials & the MaterialTable](materials.md) for
the inversion.

This was not always true. An older Tsuji *polynomial* fit overshot below the
elastic limit — nominal 0.95 realized ≈ 0.965 — which is why earlier docs warned
that the input was *not* the realized COR. That bias was removed (commit
`c3ecd67`). With the exact inversion, `bench_hertz_rebound` now realizes a COR
essentially equal to the nominal input:

| Nominal `restitution` | Realized COR |
|---|---|
| 0.95 | ≈ 0.950 |
| 0.90 | ≈ 0.902 |
| 0.70 | ≈ 0.717 |

The Hooke path uses the exact analytic `β`, so it realizes the nominal COR by
construction. The residual offset at lower `e` on the Hertz path is the
discretization of a stiff, short contact, not a model bias. This matters for any
quantitative calibration: near the elastic limit the granular temperature scales
as `T* ∝ 1/(1 − e²)`, so even a small input-vs-realized gap would throw stresses
and cooling rates off — which is precisely why the inversion is exact.
