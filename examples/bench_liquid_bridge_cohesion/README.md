# Pendular Liquid-Bridge Cohesion

Validates the opt-in `liquid_bridge_model = "willett2000"` normal cohesion force.
The single-contact arm separates two equal spheres and compares DIRT's tensile
normal force to Willett et al.'s closed-form pendular bridge expression,

`F = 2*pi*R*gamma*cos(theta) / (1 + 1.05*s_hat + 2.5*s_hat^2)`,

with `s_hat = s*sqrt(R/V)` and an explicit rupture distance.
The force arm runs with `skin_fraction = 1.0`, so bridge-only samples require
DIRT's DEM cutoff radius to include the configured rupture distance.

![liquid bridge force](plots/bridge_force.png)

*Measured DIRT force versus the Willett et al. (2000) closed form. Latest run:
PASS, maximum relative error `4.24e-13` inside the `1.0e-9` gate; force drops at
the configured rupture distance under the default neighbor skin.*

The bulk arm reuses the lifted-cylinder angle-of-repose protocol at small scale
and sweeps liquid bridge volume. The check is deliberately a trend check:
pendular liquid content should raise the static angle of repose relative to the
dry case, as reported for wet granular piles by Hornbaker et al. (1997) and
Tegzes et al. (1999).

![wet repose trend](plots/wet_repose_trend.png)

*Static angle of repose for dry, low-liquid, and higher-liquid heaps. Latest run:
PASS, high-liquid angle exceeds dry by `2.71 deg` against a `2.00 deg` trend
gate.*

The sweep also runs a dry identity guard: the default dry config and an explicit
`liquid_bridge_model = "off"` config with nonzero liquid parameters produce
byte-identical force traces.

Run:

```bash
python3 examples/bench_liquid_bridge_cohesion/sweep.py
```

References: Willett, Adams, Johnson, Seville (2000), "Capillary bridges between
two spherical bodies"; Rabinovich, Esayanur, Moudgil (2005), "Capillary forces
between two spheres with a fixed volume liquid bridge"; Hornbaker et al. (1997)
and Tegzes et al. (1999) for the wet-granular angle-of-repose trend.
