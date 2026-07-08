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
source ~/projects/.build-env
$BENCH_PYTHON examples/bench_curtis_wet_fiber_breakage/sweep.py
```

![breakage vs impact velocity](plots/breakage_vs_impact_velocity.png)

*DIRT breakage ratio and minimum largest-fragment mass ratio over impact
velocity. The quantitative gate checks every point against coarse Yang/Curtis
Fig. 13 modified-Weber reference bands and requires at least one BPM bond to
break. Latest run: PASS.*

![modified Weber trend](plots/weber_trend.png)

*Same measurements against modified Weber number with the Yang/Curtis Fig. 13
low-We* and high-We* reference bands shaded. Latest run: PASS.*

## Gate

The gate uses a small committed reference table,
[`data/yang_curtis_reference_bands.csv`](data/yang_curtis_reference_bands.csv),
read from the low- and high-modified-Weber envelopes in Yang et al. (2019)
Fig. 13:

* low `We* <= 150` cases must remain mostly intact
  (`breakage_ratio <= 0.25`, `largest_fragment_mass_ratio >= 0.75`),
* high `225 <= We* <= 500` cases must be in the post-transition breakage band
  (`breakage_ratio >= 0.45`, `largest_fragment_mass_ratio <= 0.65`),
* every DIRT point must fall inside its Yang/Curtis band,
* the breakage ratio must not decrease and largest-fragment ratio must not
  increase over this ordered sweep,
* BPM bond breakage must be active (`max(bonds_broken) > 0`).

The table stores only coarse pass bands with citation/provenance, not copied
paper figures.
