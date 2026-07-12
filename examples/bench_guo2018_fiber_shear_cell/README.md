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

The source makes an important distinction. The **experiment** uses a Schulze
ring shear tester, but the paper's **DEM comparison** deliberately replaces
the whole ring with a small planar control cell: periodic x/z boundaries, a
weight-loaded vertically mobile lid, a lower wall translating in x, and 4-mm
upper/2-mm lower blade arrays (Computational Set-up, pp. 5--7). Therefore a
rotating annulus is not a requirement for reproducing the published numerical
protocol. The previous branch incorrectly treated the experimental apparatus
as the DEM boundary contract.

The current configuration remains blocked, for a different and source-faithful
reason: it is a fixed-boundary cylindrical cup with unbladed planes, so it is
neither the SRST nor the paper's periodic planar DEM reduction. The contract
checker rejects it before DIRT runs. `--prepare-only` remains available solely
for the Table-2 topology audit.

**What remains valid.** The preparation script auditablely reproduces the
declared Table-2 bead/bond/material/timestep representation. Its cup planform,
initial gap, and 64/96-mm population scaling are only legacy surrogate inputs,
not recovered ring-shear geometry. They cannot be used to infer Fig. 6 or 7.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi
python3 examples/bench_guo2018_fiber_shear_cell/run_case.py \
  --pressure-pa 651 --width-mm 64 --output /tmp/guo-p651-w64 --ranks 8 --prepare-only
```

An actual campaign is blocked until a DIRT input implements the paper's
periodic planar control-cell boundary conditions and source-described blade
arrays, and until its population/cell dimensions are independently recovered
rather than inherited from this legacy circular-cup surrogate. This is a
protocol/source-data issue, not a calibration parameter to tune.

`validate.py` is retained as a future fail-closed history validator, but cannot
be invoked honestly until the numerical control-cell contract above is
implemented. Its stored figure points are not a substitute for a DIRT result.

Limitations: the primary paper is available locally and establishes the DEM
boundary protocol, but the present topology audit does not recover its control
cell population or dimensions. Digitized figure points are not a substitute for
that source recovery, and bonded overlapping spheres still differ from physical
cord. There is no solver result, comparison plot, acceptance PASS, or claim of
replication. The deterministic starting placement is only a documented topology
input, not a recovered experimental microstate.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/validate.py --self-test
```
