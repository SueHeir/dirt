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

## Re-derived admission contract

An angle of repose is an outcome of the material, particle and surface state,
apparatus, deposition/release history, measurement rule, and the model.  A
numeric band cannot therefore be admitted merely because it appeared in an
archived input deck or an earlier review.  If a maintained SPH implementation
later proposes this calibration, it must establish the following *before* it
chooses a coefficient or an acceptance interval:

1. Obtain an inspectable, protocol-matched external glass reference with the
   material, particle/surface condition, vessel, formation/release procedure,
   estimator, and uncertainty recorded.  Only that record may define the
   numerical target and its uncertainty interval.
2. Pre-register the parameter domain, independent replicate count, estimator,
   stopping/arrest rule, and decision rule without using results from the
   proposed sweep.  The rule must report every eligible candidate, including
   failures; it may not substitute a fitted interval for the external one.
3. Run the declared sweep from solver-produced deposits while retaining the
   complete formation history.  Publish a hash-bound ledger of inputs, solver
   revision, seeds, raw outputs, fitted angles, and the command that rebuilds
   every plot from that ledger.
4. Check at least one independent failure mode that the fitting routine cannot
   make pass: for example, a blinded manual/geometry-based angle measurement
   on retained snapshots, and a separately implemented replay or
   implementation-parity comparison.  Neither is a physical reference.
5. Resolve or quantitatively bound any cross-code discrepancy before a
   coefficient is transferred across a solver or substrate boundary.

The former 22–26 degree range is an **unadmitted archival claim**, not a
preserved target: the available citation check does not connect it to a
protocol-matched glass experiment.  Likewise, monotonicity and a replicate
spread limit are useful pre-registered diagnostics, but no particular form or
tolerance is adopted here.  This page is not an executable pass gate and does
not assert that any of the requirements have been met.

## Adversarial evidence disposition

The table below distinguishes evidence that can establish a result from checks
that can only rule out a bad claim.  It is intentionally not a checklist whose
completion may be inferred from this page.

| Question | Minimum independent evidence | Current disposition |
| --- | --- | --- |
| Is the cited literature record the claimed publication? | A resolvable DOI and bibliographic metadata checked independently of the solver | **Negative control only.** Crossref resolves Zhou, Wright, Yang, Xu, and Yu, *Rolling friction in the dynamic simulation of sandpile formation*, *Physica A* 269 (1999), DOI [10.1016/S0378-4371(99)00183-1](https://doi.org/10.1016/S0378-4371(99)00183-1). Bibliographic identity does not establish material or protocol comparability. |
| Does that source establish the historical 22–26 degree glass claim for this protocol? | Recoverable observations plus matching glass material, particle/surface state, vessel, formation/release procedure, estimator, and uncertainty | **No.** The source is a different sandpile simulation and is retained only as an incompatibility control; it is not a numerical target. |
| Does a solver reproduce a retained campaign? | A complete, hash-bound input/seed/raw-output/fitted-angle ledger and a regenerated plot | **No current DIRT surface.** The executable and its ledger were removed with DIRT #163, so archived branch plots cannot answer this question. |
| May a fitted coefficient cross from DIRT to an SPH solver? | The preceding evidence, a pre-registered result ledger, and a qualified independent-code reconciliation | **No.** There is no calibrated coefficient to transfer. |

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

To reproduce this disposition against that exact snapshot:

```bash
git -C ~/projects/dev_soil_sph ls-tree -r --name-only \
  b4997642678bf8072baa4d98be60429a4dfc59a9 -- examples | \
  rg -i 'repose|calibrat|glass' || true
git -C ~/projects/dev_soil_sph grep -n -i -E \
  'rolling friction|angle.of.repose|mu_r|μ_r' \
  b4997642678bf8072baa4d98be60429a4dfc59a9 -- \
  Cargo.toml crates examples README.md docs ':!docs/dem-campaign-dirt.md'
```

The first command returns no candidate validation surface.  The second finds
only the v0/deferred rolling-friction statement and unrelated `mu_ref` names;
neither is a target, a calibration, or a replacement pass criterion.

This retirement note and the decision to withhold a coefficient were prepared
with AI assistance.  They are a provenance boundary, not experimental work,
an external-reference validation, or a substitute for an independent review.
The DOI identity check above was performed from Crossref metadata on 2026-07-19;
it is not a reading of the paper's numerical observations and does not reduce
the evidence requirements in this record.
