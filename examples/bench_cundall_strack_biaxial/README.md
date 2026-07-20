# Cundall--Strack-style quasi-2-D walled assembly

This deterministic DIRT case is a quasi-2-D, 197-grain dense-assembly
compression benchmark.  It records only the dense late-loading part of a
loose-insertion-to-compaction path, rather than interpreting the sparse
insertion transient as a biaxial specimen. It is motivated by the wall-resultant apparatus in
Cundall & Strack (1979), but it does not claim to reconstruct that paper's
unpublished coordinates or controller history. The original paper is
provenance for particle count and friction only; it is not a trajectory oracle.

```bash
$BENCH_PYTHON examples/bench_cundall_strack_biaxial/sweep.py
```

![Measured DIRT wall-cell diagnostics](plots/stress_volume_response.png)

The top panel exposes DIRT's volumetric and contact/fabric observables; it does
not imply that the primary paper published corresponding curves. The lower
panel is the recorder-integrity check described below. It deliberately contains
no external curve: no available artefact measures the same cell and all
acceptance observables.

The executable check verifies forward axial and lateral compression, a dense
late-stage contact network, a positive platen reaction, finite output, and—independently from the precomputed
column—that `F_H/F_V` equals the recorded x-wall mean divided by the recorded
platen mean. It separately audits the committed primary-source transcription
of Fig. 10 A=`0.39` and B=`0.33`. It does not compare a selected DIRT time
window with those numbers.

"Dense" is an executable integrity condition, not a label: every recorded
state must have `phi >= 0.50` and coordination `>= 4.0`. These floors protect
against the prior loose-insertion transient (which could have only a few
contacts); they are not fitted to the Cundall--Strack force snapshots and do
not turn the unavailable source trajectory into an external pass condition.

## Comparison and limitation

Fig. 10 provides two force-network configurations (`F_H/F_V = 0.39` and
`0.33`). The paper additionally says that A→B increased vertical force and
decreased horizontal force by 4.348%, and gives a 197-disc radius census; the
driver audits that transcription. It does not publish the digitized disc
coordinates or wall positions used by BALL, nor a source strain/time coordinate
that identifies which DIRT state should be A or B. A fixed late-loading window
therefore cannot be a defensible Fig. 10 target.
The source also used 2-D discs and beam-wall control; this case uses a
one-grain-deep 3-D sphere slice, Hertz contacts, and prescribed wall velocity.

The primary source does not publish a state-registered trajectory. Therefore
its two Fig. 10 snapshots are never used as a curve, interpolation target, or
tolerance. A previously bundled LAMMPS trace was removed because its periodic
2-D virial measurement and preparation do not match this finite-wall 3-D
resultant protocol. Plotting two unlike traces as an "independent comparison"
would make the missing common observable look like a numerical disagreement.

The executable also inventories the full acceptance evidence rather than
checking only snapshot ratios. The primary source has no registered state
coordinate, stress/deviatoric path, volumetric/dilatancy path, or
contact/fabric-evolution series for this assembly. A state map alone would
therefore still not satisfy the replication contract.

## External-admission gate

`sweep.py external` evaluates the cited source before a response is read or a
tolerance is considered. A candidate must provide a registered state map plus
stress/deviatoric, dilatancy, and fabric/contact series at common states. The
primary source lacks all four; the command exits INELIGIBLE. This does not meet
the frozen acceptance criterion: a protocol-comparable external trajectory
with stress, dilatancy, and fabric/contact paths remains required.

For a future positive comparison, `candidate_package.py` additionally requires
SHA-256-bound independent solver input, protocol ledger, and response CSV
artifacts.  An all-`yes` protocol ledger alone is insufficient: the response
must expose monotone registered states plus stress-ratio, volumetric-strain,
and fabric-anisotropy series.  This prevents a prose assertion of equivalence
from being scored as a replication.

## Authorship and validation boundary

This benchmark and its evidence-contract code were AI-assisted. The committed
DIRT execution demonstrates only the stated solver behavior on this host; it
is not independent experimental confirmation or a substitute for the
unavailable primary-source trajectory. Any future positive claim requires an
independently traceable external data set and review of its protocol equivalence;
this code supplies no tolerance or criterion that can turn the present evidence
into a pass.

The cell observables use the instantaneous named x-wall positions, rather than
the fixed domain/decomposition bounds. This matters for the diagnostic itself:
the prescribed x walls move during loading, so volumetric strain and solid
fraction must follow their live separation. It improves the internal
measurement but does not repair the missing external trajectory.
