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

## Adversarial evidence disposition

The table below distinguishes evidence that can establish a result from checks
that can only rule out a bad claim.  It is intentionally not a checklist whose
completion may be inferred from this page.

| Question | Minimum independent evidence | Current disposition |
| --- | --- | --- |
| Is the cited literature record the claimed publication? | A resolvable DOI and bibliographic metadata checked independently of the solver | **Negative control only.** Crossref resolves Zhou, Xu, Yu, and Zulli, *Rolling friction in the dynamic simulation of sandpile formation*, *Physica A* 269 (1999), DOI [10.1016/S0378-4371(99)00183-1](https://doi.org/10.1016/S0378-4371(99)00183-1). Bibliographic identity does not establish material or protocol comparability. |
| Does that source establish the 22–26 degree glass target for this protocol? | Recoverable observations plus matching glass material, particle/surface state, vessel, formation/release procedure, and estimator | **No.** The source is a different sandpile simulation and is retained only as an incompatibility control; it is not a numerical target. |
| Does a solver reproduce a retained campaign? | A complete, hash-bound input/seed/raw-output/fitted-angle ledger and a regenerated plot | **No current DIRT surface.** The executable and its ledger were removed with DIRT #163, so archived branch plots cannot answer this question. |
| May a fitted coefficient cross from DIRT to an SPH solver? | The preceding evidence, unchanged monotonicity/spread/band gates, and a qualified independent-code reconciliation | **No.** There is no calibrated coefficient to transfer. |

This separation is a deliberate adversarial check: a DOI lookup, a unit-level
contact-law parity result, or a visually plausible archived heap must never be
promoted into experimental provenance or a calibration closure.  A future
campaign must be capable of failing each row independently.

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
The DOI identity check above was performed from Crossref metadata on 2026-07-19;
it is not a reading of the paper's numerical observations and does not reduce
the evidence requirements in this record.
