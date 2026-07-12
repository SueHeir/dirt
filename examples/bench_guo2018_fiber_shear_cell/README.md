# Guo flexible-fibre normal-stress-fixed shear-cell campaign

This directory contains a source-parameterized DIRT fibre-topology audit and
an explicitly blocked surrogate protocol; it does not contain a bundled
physics result. Guo
et al., *AIChE Journal* **65** (online 2018; issue 2019),
doi:10.1002/aic.16397, report the
rubber-cord measurements digitized in `data/guo_2019_rubber_cord.csv`: Fig. 6
steady shear stress and Fig. 7 solid fraction at 651, 1735, and 3470 Pa.

The prepared 64-mm and 96-mm-diameter cup populations and the material,
bond, and timestep fields use the reported rubber-cord Table 2
representation: 17 beads of 2.4-mm diameter at 1.2-mm spacing, 21.6-mm fibre
length, 1157.5-kg/m³ density, 6.28-MPa modulus, and a 2.3e-7-s time step.
Prepare and audit one without running a solver:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py \
  --width-mm 64 --output /tmp/guo64
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --audit /tmp/guo64
```

The current configuration is **not a valid Guo experimental protocol**. The
publisher-indexed abstract identifies the apparatus as a **Schulze ring shear
tester**. DIRT's configuration instead has a translating planar base inside a
cylindrical cup. Those boundary motions are not interchangeable: a numerical
agreement would not validate the experiment. `source_contract.py` queries the
independent Crossref DOI record and `run_case.py` blocks every solver run; it
permits `--prepare-only` solely for the Table-2 topology audit.

**What remains valid.** The preparation script auditablely reproduces the
declared Table-2 bead/bond/material/timestep representation. Its cup planform,
initial gap, and 64/96-mm population scaling are only legacy surrogate inputs,
not recovered ring-shear geometry. They cannot be used to infer Fig. 6 or 7.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi
python3 examples/bench_guo2018_fiber_shear_cell/run_case.py \
  --pressure-pa 651 --width-mm 64 --output /tmp/guo-p651-w64 --ranks 8 --prepare-only
```

An actual campaign is blocked until (1) the ring's annular dimensions and
rotating-member motion are independently recovered from the paper or its
supplement and (2) DIRT can impose that annular rotational shear while measuring
the correct torque/area stress. This is a missing capability/source-data issue,
not a calibration parameter to tune.

`validate.py` is retained as a future fail-closed history validator, but cannot
be invoked honestly until the apparatus contract above is implemented. Its
stored figure points are not a substitute for ring geometry or a DIRT result.

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

Limitations: the full paper figures/supplement are not available locally for
independent geometric recovery; digitized points have not been sufficient to
establish the annular kinematics; and bonded overlapping spheres still differ
from physical cord. There is no solver result, comparison plot, acceptance
PASS, or claim of replication. The deterministic starting placement is only a
documented topology input, not a recovered experimental microstate.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/validate.py --self-test
```
