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

`dirt_wall` supports rigid named plane-wall assemblies: follower planes
translate with the displacement and velocity of one driver, including a
force-servo driver.  The four blade-face entries in `config.toml` declare
eight 8-mm-pitch copies each; the generic input loader materializes all 32
finite faces before integration. Lower blades follow the translating base and
upper blades follow the mobile lid. This is a wall geometry/motion facility,
not a fitted rheology closure.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/run_campaign.py --output /tmp/guo-campaign --ranks 1
python3 examples/bench_guo2018_fiber_shear_cell/validate.py \
  --case 651:64:/tmp/guo-campaign/p651_w64 \
  --case 1735:64:/tmp/guo-campaign/p1735_w64 \
  --case 3470:64:/tmp/guo-campaign/p3470_w64 \
  --case 651:96:/tmp/guo-campaign/p651_w96 \
  --case 1735:96:/tmp/guo-campaign/p1735_w96 \
  --case 3470:96:/tmp/guo-campaign/p3470_w96
```

No solver history, comparison plot, or replication PASS is committed by this
revision.  The validator accepts only solver-receipted histories with measured
normal-load qualification, a post-drive steady-strain window, all three loads,
and the independent Fig. 6/7 digitization plus the 64/96-mm sensitivity check.
The authorship assistance for this implementation is AI-assisted; the absence
of a completed full campaign is a validation limit, not evidence of agreement.
