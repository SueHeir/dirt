# Independent SDS oracle and feasibility boundary

This calibration is checked against LAMMPS (22 Jul 2025 Update 4), not merely
against another DIRT calculation.  The relevant independent implementation is
`src/GRANULAR/gran_sub_mod_rolling.cpp`, `GranSubModRollingSDS::calculate_forces`,
and the documented model is `doc/src/pair_granular.rst`, “rolling sds”.  It
defines a rolling *pseudo-force*

```
F_roll = -k_roll xi_roll - gamma_roll v_roll,
v_roll = -R_eff ((omega_i - omega_j) x n),
|F_roll| <= mu_roll |F_n|,
tau_roll = R_eff (n x F_roll).
```

For a pair indexed as DIRT `i,j`, DIRT stores `n_DIRT=(x_j-x_i)/r`, while
the LAMMPS source's contact normal is `n_LAMMPS=-n_DIRT`.  DIRT therefore uses
the LAMMPS pseudo-velocity and history after normal conversion, and applies
`tau_i=-R_eff (n_DIRT x F_roll)` for the final SDS couple.  This conversion is
covered for both Hertz and Hooke contact paths, including the capped-history
reconstruction; it is a software-parity check, not a repose calibration.

Thus `k_roll` has units N/m and `gamma_roll` N s/m.  The supplied DIRT and
LAMMPS inputs use the same numerical pseudo-force coefficients.  `sweep.py
external` is intentionally an adversarial receipt: it executes LAMMPS, stores
hashes of its generated input and final dump, and applies the same geometry
extractor.  It is not an experimental reference and cannot close the glass
calibration because its `fix pour` formation differs from DIRT's pre-filled,
lifted cylinder.

There is also a model-independent static boundary.  With no cohesion, a grain
sliding on a plane is supported only while `tan(theta) <= mu_p`.  Here
`mu_p = 0.16`, so the sliding-only limit is `atan(0.16) = 9.09 deg`; the retained
experimental interval begins at 22 deg (`tan(22 deg) = 0.404`).  An SDS rolling
couple resists relative rotation but is capped separately; it is not an
additional tangential Coulomb force.  Therefore a 22--26 degree result from this
violent lift-and-collapse case must be demonstrated by qualified solver output
and the external oracle, not produced by changing the fit, target, or tolerance.

AI-assisted analysis and implementation.  The LAMMPS comparison is an
independent software check, not experimental validation; no protocol-matched
glass experiment or successful calibration is claimed by this file.
