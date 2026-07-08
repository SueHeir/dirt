# Rod-Like Shear Aspect-Ratio Benchmark

This benchmark is a compact replication of the glued-sphere rod trends reported by
Guo, Wassgren, Ketterhagen, Hancock, James, and Curtis, "A numerical study of
granular shear flows of rod-like particles using the discrete element method",
JFM 713 (2012). DIRT runs frictionless glued-sphere rods in a periodic
Lees-Edwards shear cell at fixed equivalent-volume diameter and dilute solid
fraction.

The regression gate checks two dimensionless published trends for elongated rods:
the Bagnold-normalized pressure `p/(rho d_v^2 gamma_dot^2)` and shear stress
`|sigma_xy|/(rho d_v^2 gamma_dot^2)` decrease as aspect ratio increases in the
dilute regime. The apparent friction panel overlays an approximate,
self-digitized trace from Guo et al. Fig. 18 for glued-sphere rods near
`nu = 0.05`; it is shown as context with the measured long-axis alignment.

```bash
python3 examples/bench_rod_shear_aspect_ratio/sweep.py
```

![Rod shear aspect-ratio trends](plots/rod_shear_aspect_ratio.png)

*DIRT glued-sphere rod shear compared with Guo et al. (2012) trends. Latest run:
PASS for decreasing dilute Bagnold-normalized pressure and shear stress with
aspect ratio for elongated rods.*
