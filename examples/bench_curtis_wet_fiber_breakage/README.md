# bench_curtis_wet_fiber_breakage — wet fiber agglomerate impact

This validation reproduces the central trend from Yang, Bello, Buettner, Guo,
Wassgren, and Curtis (AIChE Journal 65(8), 2019): a wet flexible-fiber
agglomerate impacting a plane breaks more as impact velocity, and therefore
modified Weber number, increases.

The DIRT case is deliberately small enough for regression. Three bonded fibers
are initialized as a compact wet agglomerate, each fiber is represented by five
bonded DEM beads, the Willett pendular liquid-bridge model holds neighboring
fibers together, and BPM breakage is active on the fiber bonds. The recorder
measures the paper's two breakage observables: bridge-contact breakage ratio and
minimum largest-fragment mass ratio.

```bash
python3 examples/bench_curtis_wet_fiber_breakage/sweep.py
```

![breakage vs impact velocity](plots/breakage_vs_impact_velocity.png)

*DIRT breakage ratio and minimum largest-fragment mass ratio over impact
velocity. The quantitative gate requires the breakage ratio to rise, the largest
fragment to shrink, both to correlate with modified Weber number, and at least
one BPM bond to break. Latest run: PASS.*

![modified Weber trend](plots/weber_trend.png)

*Same measurements against modified Weber number, following Yang et al.'s
energy-ratio framing. Latest run: PASS.*

## Gate

The trend gate is quantitative but avoids committing copyrighted figure data:

* breakage-ratio span across the velocity sweep must be at least 0.22,
* largest-fragment mass-ratio drop must be at least 0.28,
* both metrics must have `R2 >= 0.75` against modified Weber number,
* each curve may have at most one small off-trend step,
* BPM bond breakage must be active (`max(bonds_broken) > 0`).

These checks encode the paper's reported trends without storing its plotted
points in this repository.
