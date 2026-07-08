# Haff Cooling Benchmark — Rod Clumps

Validates free cooling of a granular gas of **rigid rod-shaped clumps** against
Haff's law. Each rod is 4 sub-spheres in a line, giving a highly asymmetric
inertia tensor (Ix ≪ Iy ≈ Iz) that exercises the angular-momentum integration. A
periodic box of rods is given a random velocity field and left to cool through
inelastic collisions.

## Physics

For inelastic particles with a **velocity-independent** restitution (constant
`e`), the granular temperature follows Haff (1983):

```
T(t) = T0 / (1 + t/tc)^2          →   late-time log-log slope = -2
```

DIRT's contact gives a constant `e`, so the rod gas obeys this `t⁻²` law. The
diagnostic for "is this Haff's law?" is that **1/√T is linear in t** (this holds
across the whole decay), whereas the bare log-log slope only reaches −2
asymptotically at `t ≫ tc`. A dilute gas cools slowly, so a finite run reaches a
moderate `t/tc` and a shallower slope; the R² of the 1/√T fit is the robust
validation. The 4-sphere rod, with its larger `tc` and slower cooling, now runs long
enough to reach `t/tc ≈ 13`; the fitted late-window slope is about −1.84, inside the
unchanged `−2.3 < slope < −1.6` gate while still short of the formal `t/tc → ∞`
limit.

If this benchmark reports a slope near -1.60 at only `t/tc ≈ 5`, it is running
the older 700k-step setup. That run length is too short for the rod case's
finite-window slope gate even though the linear Haff fit is already good; the
current benchmark uses 1.6M steps so the same unchanged gate is reached honestly.
For scheduled harness failures, also check that the harness checkout is actually
on current `origin/main`; a stale detached checkout can still run the old 700k-step
constants after the benchmark fix has merged.

## Setup

| Property | Value |
|----------|-------|
| Rod | 4 sub-spheres, r_sub = 0.5 mm, centers at ±0.4, ±1.2 mm along x (half-length 1.7 mm) |
| Count | 500 rods in a 40 mm periodic cube |
| Material | E = 70 MPa, ν = 0.245, e = 0.926, μ = 0.16, no rolling friction |
| Initial field | random, σ = 0.5 m/s per component |
| Contact | Hertz normal + Mindlin tangential, viscoelastic (tsuji) damping |

## How to Run

```bash
python3 examples/bench_rod_haff_cooling/sweep.py            # generate → start → graph
python3 examples/bench_rod_haff_cooling/sweep.py generate   # write data/rod.mol + in.lammps
python3 examples/bench_rod_haff_cooling/sweep.py start      # build, run DIRT + LAMMPS
python3 examples/bench_rod_haff_cooling/sweep.py graph      # validate vs Haff + plot
```

### LAMMPS cross-check

If a LAMMPS binary with the needed packages is on `PATH`, `start` also runs the
same gas in LAMMPS as **rigid multisphere** (`fix rigid/small molecule` with an
auto-generated rod molecule template, `pair_style granular` with matched Hertz +
Mindlin + tsuji damping, intra-clump neighbor exclusion) and overlays its cooling
curve. LAMMPS is optional — without it, or when the local binary lacks the
MOLECULE/rigid support, the benchmark validates DIRT alone.

**Caveats** (a cooling-*law* comparison, not point-by-point): the codes use
different rigid-body integrators and clump-contact handling, and a many-body gas
is chaotic. The LAMMPS total granular temperature comes from the total clump
kinetic energy (for a rigid body the summed sub-sphere KE equals body
translational + rotational KE); curves are normalized and compared by the Haff
fit and the −2 slope.

## Validation Criteria

| Check | Pass condition |
|-------|----------------|
| Finite, non-negative temperatures | all `T` finite and ≥ 0 |
| Cooling | `T_final < T_initial` |
| No energy growth | `max(T) < 1.5 T₀` |
| Haff's law holds | `1/√T` linear in `t`, R² > 0.99 |
| Late-time decay | log-log slope ∈ [−2.3, −1.6] (→ −2 as `t/tc` grows) |

## Expected Plot

![Haff cooling](plots/haff_cooling.png)

*Left:* the cooling law **past the rotational-equilibration transient**. DIRT
starts at `T_rot=0` and LAMMPS starts already spinning (from the rigid-body
velocity projection), so the early transient is discarded and the equilibration point
is treated as a fresh start (time re-zeroed, `T` re-normalized there). The DIRT
curve follows the Haff fit to `t/tc ≈ 13`, where the fitted late-time slope is
about −1.84 and still approaching the −2 asymptote. If the optional LAMMPS run is
available, it is overlaid on the same panel. *Right:* DIRT's full energy partition
(translational and rotational), showing the start-up transient that is skipped on
the left.

## References

1. P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", *J. Fluid Mech.* 134 (1983) 401–430.
