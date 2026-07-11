# Guo 2019 rubber-cord shear-cell acceptance gate

This directory is an acceptance contract, not a bundled physics result.  Guo
et al., *AIChE Journal* **65** (2019), doi:10.1002/aic.16397, report the
rubber-cord measurements digitized in `data/guo_2019_rubber_cord.csv`: Fig. 6
steady shear stress and Fig. 7 solid fraction at 651, 1735, and 3470 Pa.

The old review history tried to make a replication claim from a short replay
whose lid had zero reaction while the packing fell.  That is neither a loaded
cell nor an experimental comparison.  This replacement makes that failure
impossible to label PASS: `validate.py` accepts only solver-written histories
with a gravity-settle stage, at least 40 measured-load samples within 15% of
the requested pressure, and at least 40 driven-shear samples at engineering
strain >= 0.50.  The final shear window must retain the measured load.

Run a complete campaign with a solver that records the required CSV columns,
then validate its six histories (the 96-mm cases contain 750 17-bead fibres):

```bash
python3 examples/bench_guo2018_fiber_shear_cell/validate.py \
  --case 651:64:generated/p651_w64/cell_history.csv \
  --case 1735:64:generated/p1735_w64/cell_history.csv \
  --case 3470:64:generated/p3470_w64/cell_history.csv \
  --case 651:96:generated/p651_w96/cell_history.csv \
  --case 1735:96:generated/p1735_w96/cell_history.csv \
  --case 3470:96:generated/p3470_w96/cell_history.csv
```

The 64-mm observations are compared to the external experimental digitization
(30% relative shear-stress and 0.05 absolute solid-fraction tolerances).  The
96-mm observations are an independent finite-size check with the same bands;
they are not fitted to a second reference curve.  The validator has no
hand-entered DIRT results and fails when a history is missing, incomplete, or
does not reach the predeclared state.

Digitization error, finite size, wall roughness, and representing rubber cord
as bonded spheres remain limitations.  A passing gate supports this specified
model-to-experiment comparison; it does not establish material calibration or
general rubber-cord mechanics.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/validate.py --self-test
```
