# bench_curtis_wet_fiber_breakage — wet fiber agglomerate impact

This validation checks the central trend from Yang, Bello, Buettner, Guo,
Wassgren, and Curtis (AIChE Journal 65(8), 2019): a wet flexible-fiber
agglomerate impacting a plane breaks more as impact velocity, and therefore
modified Weber number, increases.

The DIRT case is deliberately small enough for regression. Three bonded fibers
are initialized as a compact wet agglomerate, each fiber is represented by five
bonded DEM beads, the Willett pendular liquid-bridge model holds neighboring
fibers together, and BPM breakage is active on the fiber bonds. The recorder
measures the paper's two breakage observables: liquid-bridge fiber-fiber contact
breakage ratio and minimum largest-fragment mass ratio. It counts only neighbor
pairs for which DIRT's configured Willett liquid-bridge force is active, then
forms fragments from those active fiber-fiber contacts plus intact intra-fiber
BPM bonds.

```bash
source ~/projects/.build-env
$BENCH_PYTHON examples/bench_curtis_wet_fiber_breakage/sweep.py
```

![breakage vs impact velocity](plots/breakage_vs_impact_velocity.png)

*DIRT breakage ratio and minimum largest-fragment mass ratio over impact
velocity. The quantitative gate checks every point against a digitized
Yang/Curtis Fig. 13 modified-Weber master curve and requires at least one BPM
bond to break. Latest run: PASS.*

![modified Weber trend](plots/weber_trend.png)

*Same measurements against modified Weber number with the digitized Yang/Curtis
Fig. 13 reference. Latest run: PASS.*

## Gate

The gate uses a committed reference table,
[`data/yang_curtis_fig13_digitized.csv`](data/yang_curtis_fig13_digitized.csv),
read from Yang et al. (2019) Fig. 13:

* every DIRT point is linearly compared with the digitized Fig. 13 master curve
  at the same modified Weber number,
* the breakage-ratio absolute error must be `<= 0.30`,
* the largest-fragment mass-ratio absolute error must be `<= 0.40`,
* the breakage ratio must not decrease and largest-fragment ratio must not
  increase over this ordered sweep,
* BPM bond breakage must be active (`max(bonds_broken) > 0`).

The latest comparison is written to
[`data/reference_comparison.csv`](data/reference_comparison.csv). This remains a
small regression-scale agglomerate rather than the full Yang/Curtis Table 2
`Np = 664`, centripetally prepared case; the reference overlay is used to keep
the smoke case tied to the published modified-Weber trend without copying paper
figures into the repository.
