# Coarse-19 Yarn Tensile Pull — DIRT reproduction of the LAMMPS `bpm/guo` run

Reproduces the LAMMPS 1-ply coarse-19 yarn tensile pull
(`fiber_config/3mdcp_coarse_grain/pull_..._1ply_..._gl10000.0_seed1337`) natively
in DIRT: 5330 bonded spheres (2 bond types), Guo BPM beam bonds with exact
`strain_weibull` breakage, `hertz/material` contact, intra-molecule contact
exclusion, and constant-velocity z-end grips. Then (Phase 2) the same yarn is
rebuilt from single-body sphero-cylinders to check it gives the same curve.

## Input data (not committed)

The start geometry + bonds come from the generated LAMMPS data file. `data/` is
git-ignored; copy the file in before running:

```bash
mkdir -p examples/bench_yarn_pull_coarse19/data
cp ../fiber_config/3mdcp_coarse_grain/ALL_START_DATA/generated_coarse19_yarn_gl10000.0_single_ply \
   examples/bench_yarn_pull_coarse19/data/yarn_gl10000_1ply.data
```

The `[[particles.insert]]` (Atoms, `bpm/sphere`: `id mol type diameter density x y z`)
and `[bonds]` (Bonds, 2 types) both read this one file. Atom IDs are preserved as
tags so the Bonds section matches the correct atoms. The file's placeholder
density (1) is overridden by the config `density = 9.93245e-07`.

## Parameter provenance (all trace to the LAMMPS input + generator)

| Quantity | Value | LAMMPS source |
|----------|-------|---------------|
| Contact | hertz/material E=60.784, ν=0.4868, e=0.6, μ=0.4, tangential linear_nohistory | `pair_coeff` |
| Bond stiffness | Nstiff=60.784, Tstiff=11.743, MNstiff=116.892, MTstiff=11.743 | `bond_coeff` |
| Bond damping | β = 0.141 (normal/shear/twist/bend) | `bond_coeff` |
| Breakage | `guo_strain_weibull` σ₀=0.0135, m=18, L_calib=10000, ρ_shear=1, seed=1337 | `bond_coeff break_model` |
| Grips | type-2 bonds unbreakable; z-end groups moved ±0.0572 | `bond_coeff` (no break) + `fix move` |
| Timestep | 0.001888180053590586 | `timestep` |

**Deferred:** the soft cylindrical container wall (LAMMPS uses a *hooke* wall
while particles are *hertz*; DIRT's contact model is global). Its radial
confinement is negligible for an axial pull — add only if the curve disagrees.
Heat is also deferred (the LAMMPS `heat/flow constant` is the Cp integrator, not a
source; the field stays ~300 K, so it does not affect the stress-strain curve).

## Run

```bash
cargo run --release --example bench_yarn_pull_coarse19 --no-default-features \
    --features precision-double -- examples/bench_yarn_pull_coarse19/config.toml
```

Writes `data/stress_strain.csv` (`step, strain, stress_zz, virial_zz, bond_count`).
Strain is reported as in LAMMPS (`0.000216·step / 9000`). Stress is the Love–Weber
contact+bond virial `σ = −VirialStress/V` (= LAMMPS `compute stress/atom` summed).

## Validation

Overlay `data/stress_strain.csv` against the LAMMPS
`..._pull_atomstress-strain.dat` (strain, stress_zz, nbond). Expect a linear
elastic rise, then bond breakage near the ~1.2–1.5 % Weibull mean failure strain,
then load drop. With the same seed the per-bond Weibull thresholds are identical
bond-for-bond (see `dirt_bond::breakage::guo_weibull_u`).
