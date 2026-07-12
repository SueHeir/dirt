# Guo flexible-fibre normal-stress-fixed shear cell

Guo et al., *AIChE Journal* 65 (2019), doi:10.1002/aic.16397, is the intended
external reference. Crossref independently confirms the DOI and title, but
that only proves bibliographic identity. It does **not** authenticate the
inherited numerical geometry, material values, or Figure 6/7 values in this
directory.

The primary PDF is not in the local reference library and the publisher
returned an access challenge during this rescue. Consequently this directory
has **no source-authenticated geometry** and no runnable Guo case. Every
numerical detail in the retained configuration and preparer is an unverified
historical artifact, including the cell dimensions, wall shape, blade layout,
material constants, loading, and timestep. It is not a replication setup.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py --self-test
```

## Execution and validation boundary

No claim about the paper's wall, blade, or loading construction is made here:
the branch does not have the source artifact needed to support one. The runner
therefore rejects the historical draft before it writes a solver input.

Before any future campaign, the primary source must be recorded through the
receipt below. Only then may an implementer derive and audit wall construction,
normal-load control, domain geometry, and shear protocol; a source-equivalent
implementation must then be independently tested before a Fig. 6/7 campaign.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --require-runnable
```

`run_case.py` and `run_campaign.py`, including `--prepare-only`, are expected
to fail now. They may proceed only after a legitimate local PDF is recorded
with its SHA-256, digitizer, date, method-page evidence, and Fig. 6/7 pages in
`reference_provenance.json`; the validator also rejects an artifact that is not
a PDF. The DOI/title remains independently checked through Crossref; the
receipt is not a claim that bibliography alone proves any methods. No solver
history, comparison plot, or replication PASS is committed. The retained validator is
deliberately fail-closed and accepts only solver-receipted histories with
measured normal-load qualification, a post-drive steady-strain window, all
three loads, two independent Fig. 6/7 comparisons, and the 64/96-mm sensitivity
check.

## Evidence boundary

The six values in `data/guo_2019_rubber_cord.csv` were inherited without a
primary-PDF hash or digitizer record.  They therefore are **not currently an
external reference** and cannot be used to report this replication.  Agreement
with the printed Fig. 6 line is only an internal transcription check; it does
not authenticate the Fig. 6/7 points.

`evidence_contract.py` requires a local copy of the primary paper, a committed
SHA-256/digitizer/date record, page evidence for each method claim and both
figures, and PDF artifact check before `validate.py` will even read the
reference CSV. The runner uses the same guard before it writes input, so an
unaudited method transcription cannot be mistaken for a source-derived
campaign. The present manifest intentionally says `UNVERIFIED`, so this fails
closed:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py \
  --source-pdf /path/to/guo-aic-16397.pdf
```

Once a legitimate source copy is available, record its content hash,
digitization details, the pages used for each method claim, and the figure
pages; retain only values visibly read from Figs. 6/7, and run the validator
with `--source-pdf` for every comparison. Do not replace the source artifact
with a fitted relation or with DIRT output.

This rescue revision is AI-assisted. Its checks establish only the stated
provenance and protocol boundaries; it makes no physics-correctness claim and
does not ask a reviewer to supply one.
