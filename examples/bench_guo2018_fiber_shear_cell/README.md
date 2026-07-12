# Guo flexible-fibre normal-stress-fixed shear cell

Guo et al., *AIChE Journal* 65 (2019), doi:10.1002/aic.16397, is the intended
external reference. Crossref independently confirms the DOI and title, but
that only proves bibliographic identity. It does **not** authenticate the
inherited numerical geometry, material values, or Figure 6/7 values in this
directory.

The primary PDF is not in the local reference library and the publisher
returned an access challenge during this rescue. Consequently this directory
has **no source-equivalent geometry** and no runnable Guo case. The planar
configuration and `prepare.py` constants are retained only as unverified
historical artifacts; they are not a replication setup.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py --self-test
```

## Execution and validation boundary

The inherited notes say the paper's plates and blades are *rigidly connected
spheres*, while the current DIRT input uses smooth finite plane walls. That
claim itself awaits the local primary source, but either way the plane array is
not established as source-equivalent. The runner first requires a hash-bound
primary PDF and then rejects the smooth-wall surrogate.

The missing capability is a rigid assembly of wall spheres whose lower plate
can translate and whose upper plate can be force-servoed while retaining the
individual sphere contacts and reporting the plate resultant. That must exist
and be independently checked before a Fig. 6/7 campaign can be run.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --require-runnable
```

`run_case.py` and `run_campaign.py`, including `--prepare-only`, are expected
to fail now. They may proceed only after a legitimate local PDF is recorded
with its SHA-256, digitizer, and date in `reference_provenance.json`, and that
PDF is supplied via `--source-pdf`; source-equivalent wall geometry must then be
implemented and independently audited. No solver history, comparison plot, or
replication PASS is committed. The retained validator is deliberately
fail-closed and accepts only solver-receipted histories with measured normal-
load qualification, a post-drive steady-strain window, all three loads, two
independent Fig. 6/7 comparisons, and the 64/96-mm sensitivity check.

## Evidence boundary

The six values in `data/guo_2019_rubber_cord.csv` were inherited without a
primary-PDF hash or digitizer record.  They therefore are **not currently an
external reference** and cannot be used to report this replication.  Agreement
with the printed Fig. 6 line is only an internal transcription check; it does
not authenticate the Fig. 6/7 points.

`evidence_contract.py` requires a local copy of the primary paper and a
committed SHA-256/digitizer/date record in `data/reference_provenance.json`
before `validate.py` will even read the reference CSV. The runner uses the
same guard before it writes input, so an unaudited method transcription cannot
be mistaken for a source-derived campaign. The present manifest intentionally
says `UNVERIFIED`, so this fails closed:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py \
  --source-pdf /path/to/guo-aic-16397.pdf
```

Once a legitimate source copy is available, record its content hash and the
digitization details, retain only values visibly read from Figs. 6/7, and run
the validator with `--source-pdf` for every comparison. Do not replace the
source artifact with a fitted relation or with DIRT output.

This rescue revision is AI-assisted. Its checks establish only the stated
provenance and protocol boundaries; it makes no physics-correctness claim and
does not ask a reviewer to supply one.
