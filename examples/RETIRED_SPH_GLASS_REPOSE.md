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

1. Git history confirms that the removal commit deleted the historical README,
   runner, and model entry point, and that the case is absent from the checked-out
   HEAD.
2. Crossref and OpenAlex independently identify two superficially relevant
   citations.  Elekes & Parteli, doi:10.1073/pnas.2107965118, is a cohesive
   granular-material theory; Zhou et al., doi:10.1016/S0378-4371(99)00183-1, is a
   simulation study.  Both are rejected as a dry-glass-bead validation target.

The audit fails closed if either live catalogue is unavailable or disagrees on
identity.  A passing audit does **not** validate a repose angle, a constitutive
law, any solver, or any material coefficient.  It only establishes that there
is no admissible external target behind a restored DIRT claim.

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
