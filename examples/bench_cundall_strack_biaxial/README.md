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
python -m pip install -r examples/bench_cundall_strack_biaxial/requirements.txt
python examples/bench_cundall_strack_biaxial/source_admission.py path/to/cundall-strack-1979.pdf
```

Exit status `2` is the expected authenticated result, not a validation pass.
The source lacks a registered common state, boundary history, stress/deviatoric
and dilatancy trajectories, and contact/fabric evolution. A replication claim
requires an immutable trajectory-bearing reference package plus a
source-equivalent 2-D implementation or independent 2-D DEM parity.
