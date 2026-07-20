# LAMMPS independent reference

`in_biaxial.lmp` is an independently implemented, deterministic 2-D granular
compression path (LAMMPS 22Jul2025) used to test DIRT's normalized response.
It uses the paper's 197-grain count and interparticle friction (`tan phi=0.40`)
but does **not** claim to recreate the unpublished Fig. 10 coordinates or
stages. Its periodic lateral boundary and prescribed vertical deformation are
explicitly different from DIRT's one-grain-deep walled slice.

Regenerate with:

```bash
cd examples/bench_cundall_strack_biaxial/reference/lammps
lmp -in in_biaxial.lmp
```

The generated raw `lammps_results.csv` is committed as a solver-produced
external comparison record; it is not synthesized by the Python evaluator.
Its SHA-256 is
`0e13ab5c2c2295a2e72822ab793a57291340211ba8b8b7d4d5a1d08e53497b18`;
the benchmark verifies that fingerprint before plotting or measuring it.
