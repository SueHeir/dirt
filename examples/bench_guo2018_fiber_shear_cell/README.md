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

## Present limitation (fail closed)

The config records planar plates and the required blade dimensions, but DIRT's
current wall model has only independent planes: it cannot yet attach the full
eight-position blade arrays to the translating base and vertically mobile lid
as rigid assemblies.  A single bounded-plane pair is not an array, so
`source_contract.py --require-runnable` deliberately rejects this input.  The
runner therefore cannot produce a history that `validate.py` would compare
with Fig. 6 (steady shear stress) and Fig. 7 (solid fraction).

This is an explicit software capability boundary, not a failed numerical fit:
no solver history, comparison plot, or replication PASS is claimed.  The next
implementation must add a reusable rigid wall-assembly/motion facility (or an
equivalent source-faithful wall geometry), then run independent realizations
and compare measured lid reaction and physical cord volume fraction to the
digitized experimental points.  The existing validator remains fail-closed:
it accepts only solver-receipted histories, measured load, all three loads,
and the external Fig. 6/7 data.
