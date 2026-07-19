# Retired SPH glass angle-of-repose calibration

## Status: not a DIRT validation claim

The former `examples/SPH_glass_sphere_calibration/03_angle_of_repose`
workflow is no longer part of DIRT.  The SPH suite, its executable, its input
decks, and its archived plots were removed from `main` in DIRT pull request
#163.  Consequently, DIRT has no runnable surface on which to calibrate an
SPH rolling-friction value, and this repository makes **no** claim that a
value was calibrated, transferred, or validated.

This record deliberately does not convert an archived result into a passing
benchmark.  In particular, figures or CSVs produced on the retired branch are
not evidence for the current DIRT revision: they cannot be regenerated from
the current source tree and were obtained under a different solver/protocol
history.

## Preserved scientific contract

If a maintained SPH implementation later proposes this calibration again, its
validation plan must retain all of the following requirements rather than
retuning them around a preferred result:

1. Run the declared multi-value rolling-friction sweep with the required
   independent replicates, using solver-produced deposits and retained
   formation history.
2. Demonstrate the predeclared monotonic trend and the predeclared
   replicate-spread bound, and show at least one measured angle in the
   22–26 degree glass target band.
3. Publish a complete, hash-bound ledger of inputs, solver version, seeds,
   raw outputs, and fitted angles; regenerate the measured-versus-reference
   plot from that ledger.
4. Supply an inspectable, protocol-matched external glass reference.  A
   LAMMPS comparison may test implementation parity, but it is not a physical
   reference and cannot substitute for one.
5. Resolve or explicitly bound any cross-code discrepancy before transferring
   a coefficient across the DIRT/SPH boundary.

The 22–26 degree band, monotonicity, and spread requirements above describe
the historical acceptance contract; this page is not an executable pass gate
and does not assert that those requirements have been met.

## What would constitute new evidence

A future effort belongs with the maintained SPH solver, not with this retired
DIRT example.  It should provide a primary experimental or otherwise
protocol-matched numerical record with sufficient apparatus, material,
particle, and formation details to assess comparability.  It should then
report both favorable and adverse controls, retain raw solver output, and make
the calibration decision reproducible independently of the fitting/plotting
code.

This retirement note and the decision to withhold a coefficient were prepared
with AI assistance.  They are a provenance boundary, not experimental work,
an external-reference validation, or a substitute for an independent review.
