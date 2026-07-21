# Retired SPH glass angle-of-repose request

This is a scope record, not a validation result.

DIRT is a DEM repository.  Commit
[`f7fe1a4`](../f7fe1a44b3d744eef3b7e068c42b97c8c10ad2dc) removed the former
SPH campaign at
`examples/SPH_glass_sphere_calibration/03_angle_of_repose/`.  Consequently this
repository contains no maintained SPH executable, generated trajectory, measured
angle, selected rolling-friction value, or calibration claim for that request.

The historical README proposed a 22--26 degree band, but did not put a primary
measurement citation on the line making that numerical claim.  Its bibliography
cannot repair that missing provenance: title metadata does not establish a
dry-glass-bead apparatus, particle distribution, particle--particle or
particle--wall friction, deposition protocol, angle estimator, or uncertainty.
We therefore do not promote the archived number, a generic repose paper, or a
cross-code comparison to a target after the fact.

### External-reference check

This boundary was also checked against the archived README's first DOI,
[`10.1016/S0378-4371(99)00183-1`](https://doi.org/10.1016/S0378-4371(99)00183-1),
using the independent [Crossref record](https://api.crossref.org/works/10.1016/S0378-4371(99)00183-1)
on 2026-07-21.  Crossref identifies the title as *Rolling friction in the dynamic
simulation of sandpile formation* and lists Zhou, Wright, Yang, Xu, and Yu.  This
does not match the archived README's author list (which names Xu and Zulli in
different positions), and the Crossref record supplies no abstract or experimental
protocol.  That discrepancy is a bibliographic warning, not a finding about the
paper or its scientific merit.  In particular, neither record establishes that
the paper measured dry glass beads, nor does either one justify the 22--26 degree
band.  The DOI is retained here only so a reviewer can independently reproduce
this limited provenance check; it is not admitted as a calibration reference.

## Independent reproduction of this boundary

The following commands inspect committed Git objects, rather than this note or a
test fixture.  They should show the removal commit's deleted files and no current
file at the old campaign path:

```bash
git diff --name-status f7fe1a4^ f7fe1a4 -- \
  examples/SPH_glass_sphere_calibration/03_angle_of_repose
git ls-tree -r --name-only HEAD -- \
  examples/SPH_glass_sphere_calibration/03_angle_of_repose
```

The first command must contain the historical campaign files as deletions; the
second must be empty.  These are narrow repository facts.  They neither prove
that no suitable experiment exists elsewhere nor validate a solver, an angle, or
`mu_r`.

## Reopening the study

Any future study needs an explicit owner and a maintained SPH solver outside
DIRT, plus an independently inspectable primary dry-glass-bead experiment.  Before
running it, record the matching apparatus and bead population; particle--particle
and wall contact laws; deposition/preparation method; estimator and uncertainty;
the complete parameter sweep; and independent replicate seeds.  The comparison
must retain raw outputs and report failures as well as successes.  Do not choose a
parameter from an observed target angle or weaken a pre-registered criterion to
create closure.

No such study has been run or reviewed here.  This record was drafted with AI
assistance during PR rescue; its only claim is the auditable repository-scope
boundary above.  A new head still requires independent human review.
