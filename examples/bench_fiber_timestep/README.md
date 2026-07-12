# `bench_fiber_timestep` — bonded-sphere fiber timestep

This is a solver-backed two-mode check of the timestep recommendation in Guo,
Wassgren, Hancock, Ketterhagen, and Curtis, *Powder Technology* 249 (2013),
Eq. 41:

`dt_Guo = sqrt(2/3) L sqrt(rho/E)`.

The DIRT probe uses a 16-sphere fixed-free bonded fiber, rather than a
two-sphere oscillator. It applies small free-tip axial and transverse impulses
and records the free-tip displacement for the actual DIRT simulation. The
short-wavelength axial chain mode is the physical mesh-scale wave that Eq. 41
limits. For each mode the reference is independently assembled from the stated
element energy: axial `k_n(u_j-u_i)^2/2`, and transverse
`k_s(z_j-z_i+L(theta_i+theta_j)/2)^2/2 + k_b(theta_j-theta_i)^2/2`.
The largest generalized eigenfrequency of this fixed--free lattice gives the
central-difference bound `2/omega_max`; the transverse reference includes the
free translation and rotation coupled by the shear lever arm.

The acceptance gate requires an empirical stable/failure bracket that straddles
each spectral limit within 2%, and requires the axial 16-sphere lattice limit
to agree with Guo Eq. 41 within 1%. It therefore tests DIRT against an external
theoretical fastest-mode reference rather than fitting a threshold from DIRT.
"Explosive" means an actual
free-tip displacement exceeding 1000 m in this millimetre-scale, 1 mm/s pulse;
this six-decade numerical-growth diagnostic is intentionally far beyond the
physical response and is not a calibrated safety factor. It does not relabel a
one-bond Verlet threshold as Guo's production recommendation.

![Measured DIRT stability versus independently derived lattice limits](plots/fiber_timestep_stability.png)

Run from the repository root:

```bash
$BENCH_PYTHON examples/bench_fiber_timestep/sweep.py
```

The committed CSV and PNG are regenerated together by that command. This is a
small, undamped, fixed-free propagation experiment, not the paper's full
fiber-divergence study. It supports the axial elastic bonded-element timestep
bound and shows that the DIRT transverse mode can be stricter. It does not
determine a universal safety factor for contact-rich, damped, plastic, or
heterogeneous fiber networks.
