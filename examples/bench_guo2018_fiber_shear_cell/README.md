# Guo flexible-fibre normal-stress-fixed shear cell

Guo et al., *AIChE Journal* 65 (2019), doi:10.1002/aic.16397, compare rubber
cord experiments with a DEM reduction.  The experiment uses a Schulze ring
shear tester; the DEM does **not** simulate that annulus.  It uses a periodic
planar cell with `Lx = 64 mm`, `Lz = 36 mm`, a mobile weight-loaded lid, and a
translating lower wall.  The source specifies 4-mm upper and 2-mm lower blades
at 8-mm pitch.  Table 2 supplies the 17-bead, 2.4-mm-diameter, 21.6-mm rubber
cord representation, 500 fibres, 6.28 MPa modulus, 1157.5 kg/m3 density, and
a 2.3e-7 s timestep.

`prepare.py` now creates that *rectangular periodic* topology rather than the
old circular-cup surrogate.  The 96-mm case uses 750 fibres to preserve the
64-mm areal loading; this is a transparent derived sensitivity input, not a
reported paper population.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --width-mm 64 --output /tmp/guo64
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --audit /tmp/guo64
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi --require-runnable
```

## Execution and validation boundary

The primary paper states that its plates and blades are *rigidly connected
spheres*. The current DIRT input instead uses smooth finite plane walls. A
plane array can reproduce nominal blade locations but not the reported rough
wall contact topology, so the source-contract checker blocks solver execution.

The missing capability is a rigid assembly of wall spheres whose lower plate
can translate and whose upper plate can be force-servoed while retaining the
individual sphere contacts and reporting the plate resultant. That must exist
and be independently checked before a Fig. 6/7 campaign can be run.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --require-runnable
```

This command is expected to fail until that representation exists. No solver
history, comparison plot, or replication PASS is committed. The retained
validator is deliberately fail-closed and accepts only solver-receipted
histories with measured normal-load qualification, a post-drive steady-strain
window, all three loads, and the independent Fig. 6/7 digitization plus the
64/96-mm sensitivity check.
