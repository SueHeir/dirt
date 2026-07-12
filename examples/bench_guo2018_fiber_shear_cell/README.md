# Guo flexible-fibre normal-stress-fixed shear-cell campaign

This directory contains a source-parameterized DIRT shear-cell executable and
a fail-closed acceptance contract; it does not contain a bundled physics
result. Guo
et al., *AIChE Journal* **65** (online 2018; issue 2019),
doi:10.1002/aic.16397, report the
rubber-cord measurements digitized in `data/guo_2019_rubber_cord.csv`: Fig. 6
steady shear stress and Fig. 7 solid fraction at 651, 1735, and 3470 Pa.

The prepared 64-mm and 96-mm-diameter circular-cell populations and the executable's material,
bond, and timestep fields use the reported rubber-cord Table 2
representation: 17 beads of 2.4-mm diameter at 1.2-mm spacing, 21.6-mm fibre
length, 1157.5-kg/m³ density, 6.28-MPa modulus, and a 2.3e-7-s time step.
Prepare and audit one without running a solver:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py \
  --width-mm 64 --output /tmp/guo64
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --audit /tmp/guo64
```

To create an auditable DIRT input for one of the six cases, use the executable
path rather than hand-editing a topology or requested force. The default run
uses a cylindrical lateral boundary, an explicit gravity-settle stage, a measured
normal-load qualification before it enables the 20-mm/s lower-wall drive, and
records live lid reactions and solid fraction to `cell_history.csv`.

**Geometry correction.** Earlier revisions treated 64/96 mm as one side of a
64-by-36-mm periodic rectangle. That invented a depth and inflated the initial
64-mm physical-cord fraction to 0.408. The executable now uses DIRT's static
y-axis cylindrical sidewall and circular planform (`A = pi D^2/4`) with the
stated diameter. The initial 64-mm physical-cord fraction is 0.292, derived
only from the Table-2 cord volume, 500 fibres, the 64-mm diameter, and the
50-mm initial lid gap—not from Fig. 7. The 96-mm finite-size case uses 1,125
fibres, preserving that fraction by area scaling. This removes the geometry
rejection, but is not a completed replication: all six DIRT histories and the
external comparison still have to be run and retained.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/run_case.py \
  --pressure-pa 651 --width-mm 64 --output /tmp/guo-p651-w64 --ranks 8 --prepare-only
# Omit --prepare-only to execute the solver; the gamma >= 0.5 case is long.
```

For the complete campaign, `run_campaign.py` writes a campaign manifest and a
per-case manifest before invoking DIRT. The normal wall direction remains
local. The runner builds DIRT with its default `mpi_backend`;
it does not launch a serial-only binary under `mpirun`.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/run_campaign.py \
  --output /scratch/guo2018-campaign --ranks 8
```

The old replay had zero lid reaction while the packing fell.  That is neither
a loaded cell nor an experimental comparison.  `validate.py` makes that
failure impossible to label PASS: it accepts only solver-written histories
with a gravity-settle stage, at least 40 measured-load samples within 15% of
the requested pressure, and at least 40 driven-shear samples at engineering
strain >= 0.50.  The final shear window must retain the measured load.

Run all six solver cases with `run_case.py` (the 96-mm cases contain 1,125
17-bead fibres), then validate their histories:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/validate.py \
  --case 651:64:generated/p651_w64 \
  --case 1735:64:generated/p1735_w64 \
  --case 3470:64:generated/p3470_w64 \
  --case 651:96:generated/p651_w96 \
  --case 1735:96:generated/p1735_w96 \
  --case 3470:96:generated/p3470_w96
```

The 64-mm observations are compared to the external experimental digitization
(30% relative shear-stress and 0.05 absolute solid-fraction tolerances).  The
96-mm observations are an independent finite-size check with the same bands;
they are not fitted to a second reference curve. As a separate external audit,
the stored Fig. 6 rubber-cord shear points must agree within 3% with the
paper's printed experimental fit, `tau = 0.83 sigma_yy - 60 Pa`; this detects
an inconsistent digitization without using any DIRT output. The validator accepts a case
directory, not a free-standing CSV: it requires the `run_case.py` input
manifest and the receipt written only after DIRT exits successfully, then
checks SHA-256 digests that bind the history to those exact inputs. It fails
when a history is missing, incomplete, records a requested rather than
measured load, changes the global particle population, or does not reach the
predeclared state. This is provenance checking, not a substitute for retaining
the full independently reproducible solver run.

Digitization error, finite size, wall roughness, and representing rubber cord
as bonded spheres remain limitations.  In addition, the current Cartesian
prototype is blocked on source geometry and bead-to-cord volume equivalence;
it must not be presented as the paper's Schulze ring-shear geometry. A passing
future gate supports only the specified model-to-experiment comparison, not
material calibration or general rubber-cord mechanics. Until a DIRT run
supplies all six histories, there is no replication result or result graph to
claim. The deterministic starting placement is a documented solver initial
condition, not a claim that the paper's unreported packing microstate was
recovered.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/validate.py --self-test
```
