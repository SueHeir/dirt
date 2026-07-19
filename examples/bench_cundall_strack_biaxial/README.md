# Cundall--Strack (1979) source-admission audit

This is **not** a DIRT replication benchmark. The previous 3-D walled-sphere
surrogate was removed because it is not observable-equivalent to the source's
2-D, force-controlled disk-and-beam assembly.

The audit checksum-authenticates the primary PDF and verifies affirmative source
statements showing that its verification is qualitative, geometry was digitised,
pre-load geometry is unknown, and its numerical control differs from the
original experiment. It then deliberately returns `REPLICATION_UNAVAILABLE`.

## Run

```bash
source ~/projects/.build-env
$BENCH_PYTHON -m pip install -r examples/bench_cundall_strack_biaxial/requirements.txt
$BENCH_PYTHON examples/bench_cundall_strack_biaxial/source_admission.py path/to/cundall-strack-1979.pdf
$BENCH_PYTHON -m unittest examples/bench_cundall_strack_biaxial/test_source_admission.py
```

Exit status `2` is the expected authenticated result, not a validation pass.
The source lacks a registered common state, boundary history, stress/deviatoric
and dilatancy trajectories, and contact/fabric evolution. A replication claim
requires an immutable trajectory-bearing reference package plus a
source-equivalent 2-D implementation or independent 2-D DEM parity.

The audit also inventories the source's six published horizontal/vertical
wall-force snapshots (A/B and four BALL states). They are deliberately emitted
as `WALL_FORCE_SNAPSHOTS_ONLY`, not a response curve: the paper does not
register them to a common strain, volume, geometry, or wall-position state.
They must not be interpolated, fitted, or used as an error envelope for a DIRT
stress--strain result.

`source_admission_receipt.txt` records the authenticated result from the
read-only archived primary PDF. It is a reproducibility receipt, not reference
trajectory data. LAMMPS was checked as an available independent solver, but is
not used as a scored comparator: its documented granular styles cannot repair
the missing source state registration or turn a 3-D surrogate into the source's
2-D disk-and-finite-beam experiment.
