# Guo flexible-fibre normal-stress-fixed shear cell

Guo et al., *AIChE Journal* 65 (2019), doi:10.1002/aic.16397, is the intended
external reference. Crossref independently confirms the DOI/title, and the
local primary PDF is now hash-bound in `data/reference_provenance.json`.
That receipt authenticates the cited source pages and the approximate manual
Fig. 6/7 digitization; it does **not** authenticate the existing DIRT geometry.

The source describes a periodic planar numerical control cell, but explicitly
builds its walls and blades from rigidly connected spheres and applies the
normal stress through the *weight* of the upper wall (pp. 5-6). The retained
DIRT input instead uses plane-wall assemblies and a force-feedback lid servo.
Both changes alter the published boundary-value problem. More importantly, the
paper does not report the wall-sphere diameter/layout or upper-wall mass, so a
sphere-built, gravity-loaded body cannot be uniquely reconstructed from the
primary PDF. It is therefore **not a source-equivalent Guo case** and is
deliberately non-runnable as a replication.
The configuration is retained only as an implementation sketch; no solver
history, comparison plot, or replication PASS is committed.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/source_geometry_audit.py --verify
python3 examples/bench_guo2018_fiber_shear_cell/reconstruction_readiness.py --verify
```

`data/reconstruction_ledger.json` is a separate, page-receipted transcription
of the direct Table-2/cell facts used by the candidate.  Its audit checks those
facts against `prepare.py`, marks the 96-mm/750-fibre input as a derived
sensitivity choice, and reports the missing body inputs without turning a
default into source evidence.  Given the local primary artifact, run:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/reconstruction_readiness.py \
  --verify --source-pdf /path/to/guo-aic-16397.pdf
```

`AUDIT OK` means only that this independent transcription and the topology
constants agree.  It is explicitly not a physics result or a replication pass.

## Execution and validation boundary

The source receipt identifies what must be implemented, and also exposes the
current mismatch. The runner therefore rejects the plane-wall draft before it
writes a solver input. `prepare.py` can only write a separately labelled
candidate topology after an explicit `--candidate-only` acknowledgement; it
cannot be represented as a prepared replication case.

Before any future campaign, obtain the omitted wall-body data from the authors
or an archival supplement, then implement and audit a rigid sphere-built
wall/blade assembly and weighted, vertically free upper-wall body. A plausible
tessellation or mass is not source evidence. Then independently test that
boundary machinery before a Fig. 6/7 campaign. The PDF receipt is necessary
evidence, not a waiver of this work.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --require-runnable
```

`run_case.py` and `run_campaign.py`, including `--prepare-only`, are expected
to fail now because the geometry is non-equivalent. A future campaign must pass
the same hash-bound PDF receipt; the validator also rejects an artifact that is
not a PDF. The DOI/title remains independently checked through Crossref; the
receipt is not a claim that bibliography alone proves a DIRT implementation.
No solver history, comparison plot, or replication PASS is committed. The
retained validator is deliberately fail-closed and accepts only solver-receipted histories with
measured normal-load qualification, a post-drive steady-strain window, all
three loads, two independent Fig. 6/7 comparisons, and the 64/96-mm sensitivity
check.

For a topology-only audit (not a solver input or a result), use:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py \
  --width-mm 64 --candidate-only --output /tmp/guo64-candidate
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --audit /tmp/guo64-candidate
```

## Evidence boundary

The six values in `data/guo_2019_rubber_cord.csv` are an approximate,
AI-assisted manual digitization of the experimental rubber-cord series in
Figs. 6(b) and 7 (PDF pp. 30-31). They are hash-bound to the local source PDF
and Fig. 6 is cross-checked against its printed fit. They are not raw data and
must be independently inspected before using them to support a result.

`evidence_contract.py` requires the matching local primary paper, committed
SHA-256/digitizer/date record, page evidence for each method claim and both
figures, and PDF artifact check before `validate.py` will read the reference
CSV. The runner uses the same guard before it writes input.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py \
  --source-pdf /path/to/guo-aic-16397.pdf
```

Do not replace the source artifact with a fitted relation or with DIRT output.

This rescue revision is AI-assisted. Its checks establish only the stated
provenance and protocol boundaries; it makes no physics-correctness claim.
