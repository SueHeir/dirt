# Source and measurement boundary

## Reproduction source

The read-only primary source used for the transcription is the local research
library copy named `[Géotechnique 1979-mar vol. 29 iss. 1] A discrete numerical
model for granular assemblies ... {18139602} libgen.li.pdf`.  Its bibliographic
identity is Cundall & Strack, *Géotechnique* **29**(1), 47--65 (1979), DOI
10.1680/geot.1979.29.1.47.  The two CSV values below are transcribed from the
labels printed in Fig. 10 on p. 60; they are not digitized curve samples.

P. A. Cundall and O. D. L. Strack, *A discrete numerical model for granular
assemblies*, Geotechnique 29(1), 47--65 (1979), DOI 10.1680/geot.1979.29.1.47,
reports force-vector plots for a 197-disc photoelastic experiment (pp. 60--63).
Figure 10 labels two configurations with `F_H/F_V = 0.39` and `0.33`; Figs. 11
and 13 show the authors' BALL simulations.  The paper does not provide a
loading history, strain coordinates, controller history, or published
coordinates/wall positions that map its two images to points on a continuous
loading trace.  It does give a narrow, auditable A-to-B fact: vertical force
was increased by 4.348% while horizontal force was decreased by 4.348%.  It
also gives the 197-disc radius census and density.  Those facts are transcribed
and machine-audited in `data/cundall_strack_protocol.csv`; they are protocol
provenance, not a substitute for the missing geometry/state registration.

DIRT now records the directly comparable *kind* of quantity: the mean reaction
on the two x walls divided by the mean reaction on the moving top platen.  The
y walls merely confine the one-particle-deep slice and are excluded.  The sweep
independently recomputes the ratio from its emitted reaction columns; this is a
recorder-integrity check, not an external material-validation tolerance.

The source values remain provenance, not executable targets.  Assigning them
to preselected DIRT strain windows would manufacture the missing state
registration.  The geometry is also a 3-D sphere slice with a Hertz contact
law and velocity-driven walls, whereas the source is a 2-D disk/beam-wall
experiment.  Consequently this example establishes an auditable walled-cell
capability, not a successful numerical replication of Fig. 10.

`data/external_evidence_inventory.csv` records the four evidence classes the
goal requires: a state map, stress/deviatoric path, volumetric/dilatancy path,
and contact/fabric evolution. Each is marked absent with its source location
and limitation. The runner audits this inventory before external mode, so a
future snapshot registration cannot disguise the still-missing trajectories as
a completed replication.

## Enforced limitation

The driver now exposes a separate `external` mode. It fails closed without a
reviewed `data/source_state_registration.csv` that maps both Fig. 10 stages by
source and DIRT state coordinates and records an independent basis. No such
file is committed: the source does not supply the necessary state information.
This makes it impossible for the measurement-smoke PASS to be re-labelled as
external agreement by choosing a convenient DIRT force-ratio window.

## Independent solver countercheck

The `reference/lammps/` directory contains a LAMMPS 22Jul2025 input deck and
the raw `lammps_results.csv` it produced. It independently implements a 2-D
197-grain, `tan(phi)=0.40` compression trajectory; it does not use DIRT output
as input and is not a reconstruction of unpublished Fig. 10 geometry. Against
the deterministic DIRT run, the fixed-grid comparison reports its values in
the generated figure. No tolerance or selected window is introduced. Thus this
is independent negative evidence, not a claim
that LAMMPS validates the source paper or that DIRT reproduces it.
