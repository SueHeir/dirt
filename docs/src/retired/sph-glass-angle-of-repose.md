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

## Frozen acceptance contract and additional evidence needed

The assigned calibration retains its acceptance criterion. A future maintained
SPH implementation must run the declared multi-value rolling-friction sweep
with independent replicates, demonstrate the existing monotonic trend and
replicate-spread gate, and report at least one calibrated coefficient whose
measured angle lies in the **22–26 degree glass band**. Those requirements are
not met here and have not been relaxed or replaced.

The review history also shows why an exit status alone would be inadequate:
historical executions produced unresolved flanks, angles below the required
band, and an adverse DIRT/LAMMPS replay. Any restarted campaign must retain its
complete formation history and publish a hash-bound ledger of inputs, solver
revision, seeds, raw outputs, fitted angles, and the command that rebuilds the
measured-versus-reference figure. The figure must include the unchanged band
and replicate spread so that the criterion can fail visibly.

An angle of repose also depends on material and particle/surface state,
apparatus, deposition/release history, estimator, and model. Before any
coefficient may be represented as physically calibrated or transferred across
a solver boundary, the campaign additionally needs an inspectable,
protocol-matched external glass reference and an independent check that is not
the fitted campaign itself (for example a separately implemented replay or a
blinded geometry-based measurement of retained snapshots). A LAMMPS replay may
test implementation parity but is not a physical reference. These are
additional evidence requirements; they do not alter the existing 22–26 degree,
monotonicity, or spread acceptance gates.

This page is not an executable pass gate and does not assert that any of the
requirements have been met.

## Adversarial evidence disposition

The table below distinguishes evidence that can establish a result from checks
that can only rule out a bad claim.  It is intentionally not a checklist whose
completion may be inferred from this page.

| Question | Minimum independent evidence | Current disposition |
| --- | --- | --- |
| Is the cited literature record the claimed publication? | A resolvable DOI and bibliographic metadata checked independently of the solver | **Negative control only.** Crossref resolves Zhou, Wright, Yang, Xu, and Yu, *Rolling friction in the dynamic simulation of sandpile formation*, *Physica A* 269 (1999), DOI [10.1016/S0378-4371(99)00183-1](https://doi.org/10.1016/S0378-4371(99)00183-1). Bibliographic identity does not establish material or protocol comparability. |
| Does that source establish the 22–26 degree glass target for this protocol? | Recoverable observations plus matching glass material, particle/surface state, vessel, formation/release procedure, estimator, and uncertainty | **No.** The source is a different sandpile simulation and is retained only as an incompatibility control; it is not sufficient external provenance for the unchanged acceptance band. |
| Does a solver reproduce a retained campaign? | A complete, hash-bound input/seed/raw-output/fitted-angle ledger and a regenerated plot | **No current DIRT surface.** The executable and its ledger were removed with DIRT #163, so archived branch plots cannot answer this question. |
| May a fitted coefficient cross from DIRT to an SPH solver? | The preceding evidence, the unchanged monotonicity/spread/band gates, a complete result ledger, and a qualified independent-code reconciliation | **No.** There is no calibrated coefficient to transfer. |

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

## Relocation check (2026-07-19)

Retiring the DIRT surface does **not** relocate this calibration automatically.
As an independent cross-tier check, the tracked `dev_soil_sph` `origin/main`
snapshot `b4997642678bf8072baa4d98be60429a4dfc59a9` was inspected separately.
Its `MaterialParams` has continuum `mu_s`, `mu_2`, and `i0` parameters, but no
rolling-contact coefficient; its tracked `examples/` tree contains no
angle-of-repose or calibration executable.  Its campaign document lists
rolling friction as `mu_r = 0.0` for v0 and explicitly defers adding it.

Thus there is presently no maintained target to which the deleted DIRT command
can be moved.  A continuum static-friction parameter is not a calibrated DEM
rolling-friction coefficient, and mapping one to the other without a
model-specific, independently validated study would be an unsupported
cross-substrate transfer.  This check rules out a misleading relocation claim;
it does not validate DIRT, dev_soil_sph, or a glass angle of repose.

The removal itself was independently checked from Git objects, not inferred
from this narrative: DIRT merge `f7fe1a4` removes 141 SPH campaign files
(20,284 deleted lines) relative to its first parent, and neither current DIRT
`main` nor this PR head contains a path named
`examples/SPH_glass_sphere_calibration`.  That is a scope/availability check,
not a scientific result.  It explains why the frozen exit-zero criterion
cannot be exercised in DIRT without restoring a removed solver surface or
formally re-scoping the goal.

To reproduce this disposition against that exact snapshot:

```bash
ci/verify-retired-sph-repose.sh \
  --soil-sph ~/projects/dev_soil_sph --online
```

The verifier checks the removal and current-path absence directly, then checks
that the pinned `dev_soil_sph` snapshot has no angle-of-repose/calibration
executable and still records `mu_r = 0.0`.  With `--online` it separately
checks the Crossref bibliographic identity. Its deliberately narrow checks do
not treat continuum `mu_s`, `mu_2`, or `i0` as substitutes for the deleted
contact coefficient. None is a target, a calibration, or a replacement pass
criterion.

This retirement note and the decision to withhold a coefficient were prepared
with AI assistance.  They are a provenance boundary, not experimental work,
an external-reference validation, or a substitute for an independent review.
The DOI identity check above was performed from Crossref metadata on 2026-07-19;
it is not a reading of the paper's numerical observations and does not reduce
the evidence requirements in this record.
