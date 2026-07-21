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
   deleted every one of those files, and confirms none has returned at the
   checked-out HEAD.  This catches a restored config, pin, plot, or helper as
   well as a restored README, runner, or model entry point.
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

`audit_retired_sph_claim.py` is a separate adversarial check.  It reads the
README stored in the Git tree immediately before retirement, rather than a
present-day paraphrase. It extracts the first quoted reference title from that
immutable source and has Crossref discover a unique exact-title DOI. It then
uses that DOI to obtain direct Crossref and OpenAlex records, which must agree
with the archived title. The source-derived title describes the work as a
simulation, so it is not admitted as a primary dry-glass-bead measurement. This
avoids embedding a favoured DOI, author list, or title as a local passing
fixture. It demonstrates the narrower, important fact: the archived repository
claim did not cite a primary dry-glass-bead measurement.
It does not establish that no such measurement exists anywhere, nor does it
validate an SPH or DEM result.

## What would be required to reopen this work

A user-owned decision would be needed to select and maintain an SPH solver and
to define the product scope.  Only then could a new, separate study preregister a
single apparatus-matched protocol and independently measured inputs, run
replicates, and compare uncertainty against a primary experimental source.  No
parameter may be selected after observing the target angle.

## Authorship and review

This retirement audit and note were drafted by an AI coding agent during PR
rescue.  The code is mechanically exercised below, but the scientific judgement
is limited to rejecting inadmissible evidence; it has not been endorsed as a
physics validation.  A fresh human/independent review of the new head is still
required.
