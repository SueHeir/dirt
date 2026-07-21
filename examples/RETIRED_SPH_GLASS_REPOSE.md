# Retired SPH glass angle-of-repose campaign

The requested SPH glass-sphere rolling-friction / angle-of-repose calibration is
not part of DIRT.  DIRT is a DEM repository and commit `f7fe1a4` removed the
former executable campaign at
`examples/SPH_glass_sphere_calibration/03_angle_of_repose/`.

This is a scope decision with a scientific consequence: no current DIRT result
may be labelled an SPH calibration, a glass-bead reproduction, or a validated
value of rolling friction `mu_r`.  The retired case had neither a presently
maintained SPH solver nor a protocol-matched, independently measured target with
the particle distribution, particle--particle and particle--wall friction,
rolling law, preparation method, and uncertainty needed for such a claim.

## Evidence boundary

`verify_retired_sph_glass_repose.py --online` makes two deliberately narrow,
falsifiable checks:

1. Git history derives the complete file manifest of the historical campaign
   from the predecessor of the removal commit, confirms that the removal commit
   deleted every one of those files, and confirms neither a historical path nor
   a byte-identical historical Git blob has returned at the checked-out HEAD.
   This catches a restored or relocated config, pin, plot, or helper as well as
   a README, runner, or model entry point. It does not judge a newly authored
   implementation under another name.
2. Crossref and OpenAlex independently identify two superficially relevant
   citations.  Elekes & Parteli, doi:10.1073/pnas.2107965118, is a cohesive
   granular-material theory; Zhou et al., doi:10.1016/S0378-4371(99)00183-1, is a
   simulation study.  Both are rejected as a dry-glass-bead validation target.

The audit fails closed if either live catalogue is unavailable or disagrees on
identity, or if any historical campaign file has been restored.  A passing audit
does **not** validate a repose angle, a constitutive law, any solver, or any
material coefficient.  It only establishes that there is no admissible external
target behind the historical DIRT campaign path.  It cannot prove that a future,
renamed implementation is scientifically valid; that would require the separate
protocol and evidence prerequisites below.

## Deliberate separation from green validation CI

These evidence-boundary commands are intentionally **not** run by
`examples/ci_validation.py`.  A green CI validation summary must mean that its
listed DEM validation sweeps passed; appending a successful retirement audit to
that summary would create a misleading, self-referential route to an apparent
calibration pass.  The commands remain runnable as explicit, networked,
fail-closed provenance checks, and their output must be reported as
``INCONCLUSIVE BY DESIGN``.  Neither a passing audit nor its absence from CI
changes the frozen monotonicity, replicate-spread, 22--26 degree, raw-ledger,
closure, or criterion-visible-graph requirements.

For a cross-substrate check, run the same audit with a named maintained SPH
checkout:

```bash
python3 examples/verify_retired_sph_glass_repose.py --online \
  --soil-sph ~/projects/dev_soil_sph
```

It reads the other repository's committed `HEAD` (not uncommitted local files)
and fails if an obvious repose/calibration surface is present. It intentionally
does not treat generic friction code as a replacement and a clean result is only
a bounded source-tree observation, not evidence that no SPH solver or experiment
exists elsewhere. A candidate is a review trigger, not a calibration pass.

`audit_retired_sph_claim.py` is a separate adversarial check. It reads the README
stored in the Git tree immediately before retirement, rather than a present-day
paraphrase. It derives the numerical band and every numbered, quoted bibliography
title from that immutable source; Crossref discovers each exact-title DOI and
Crossref and OpenAlex must agree with every archived title. It also establishes a
strictly narrower source fact: the line that makes the numerical claim has no
inline numeric, DOI, or URL citation. Thus the bibliography is identity-checked,
but none of its records is silently promoted to the source of the number.

This deliberately avoids embedding a favoured DOI, author list, title, or
literature classification as a local passing fixture. Catalogue metadata and
titles cannot establish dry-bead apparatus, population, wall friction,
preparation, estimator, or uncertainty; they cannot reject or admit a primary
measurement. The audit therefore demonstrates only that the archived repository
did not connect its numerical claim to a directly inspectable source. It does
not establish that no appropriate measurement exists anywhere, nor does it
validate an SPH or DEM result.

## What would be required to reopen this work

This is the exact user-owned product/scope decision required before any new work
can claim to satisfy `dirt-sph-glass-mur-repose-calibration`:

1. Choose whether the calibration is retired, or is moved to a named, maintained
   SPH solver/repository. It must not be recreated as an SPH feature inside the
   DEM-only DIRT repository merely to preserve this goal identifier.
2. If moved, name the owning repository, executable, and material/contact model;
   select one primary dry-glass-bead experiment whose apparatus, bead population,
   particle--particle and particle--wall friction, preparation/deposition method,
   angle estimator, and reported uncertainty can be matched and inspected.
3. Preregister the resulting campaign before execution: independent replicate
   count and seed ledger, the complete rolling-friction sweep, and the unchanged
   frozen decision gates: monotonic $\theta_r(\mu_r)$, the existing per-$\mu_r$
   spread bound, and at least one result in the 22--26 degree band. The selected
   $\mu_r$, raw ledger, and a criterion-visible measured-versus-reference plot
   must be retained in the new study's closure documentation.

Only after that decision could a separate study run replicates and compare their
uncertainty against the selected primary source. No parameter may be selected
after observing the target angle. This retirement record makes no choice among
those alternatives and does not authorize a tolerance, protocol, case-set, or
target change.

## Validation-ledger status

This note is linked from `examples/VALIDATION.md` solely to make the boundary
discoverable. It is not a validation-ledger entry and carries no result figure:
there is no current SPH executable or solver-produced measurement to plot. A
future result figure is admissible only after the scope decision and evidence
requirements above are met.

## Authorship and review

This retirement audit and note were drafted by an AI coding agent during PR
rescue.  The code is mechanically exercised below, but the scientific judgement
is limited to rejecting inadmissible evidence; it has not been endorsed as a
physics validation.  A fresh human/independent review of the new head is still
required.
