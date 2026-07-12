# Guo flexible-fibre shear-cell source audit

This directory is an evidence and reconstruction audit for Guo et al.,
*AIChE Journal* 65 (2019), doi:10.1002/aic.16397. It is **not** a DIRT
replication benchmark and contains no solver input, solver history, result
plot, or replication PASS.

The acceptance contract for a future replication is non-negotiable: reproduce
the published normal-stress-fixed cell with rigidly connected-sphere walls and
blades, and a vertically free upper body whose gravity supplies the prescribed
normal load. It must then compare two solver-derived observables against the
paper across the declared cases with visible, pre-stated error bands.

The primary PDF reports the fibre representation, cell dimensions, and that a
prescribed normal stress is applied through the upper wall's weight. Thus the
assigned gravity load is derivable as `mass = normal_stress * Lx * Lz / g`:
the printed 1735-Pa case in the 64-mm x 36-mm cell requires
0.407486238532 kg at `g = 9.81 m/s²`. It does **not** report the wall-sphere
diameter or layout. Those omissions still prevent a unique sphere-built
boundary reconstruction. An analytic plane wall or force-feedback lid would be
a different problem, so this audit deliberately rejects that substitution
rather than running it under the Guo label.

## Reproducible checks

```bash
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --verify-doi
python3 examples/bench_guo2018_fiber_shear_cell/source_contract.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/evidence_contract.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/source_geometry_audit.py --verify
python3 examples/bench_guo2018_fiber_shear_cell/reconstruction_readiness.py --verify
python3 examples/bench_guo2018_fiber_shear_cell/normal_load_audit.py --self-test
```

The hash-bound receipt in `data/reference_provenance.json`, the page-receipted
`data/reconstruction_ledger.json`, and the candidate-only topology generator
make the positive evidence auditable. The 64-mm candidate can be generated and
audited without invoking a solver:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py \
  --width-mm 64 --candidate-only --output /tmp/guo64-candidate
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --audit /tmp/guo64-candidate
```

`AUDIT OK` confirms only Table-2/cell arithmetic and topology counts; it is
not a physics result. The 96-mm/750-fibre population is explicitly marked as a
derived sensitivity candidate, not a reported source input.

## What would unblock replication

Obtain the missing sphere-wall dimensions/layout from the authors, an archival
supplement, or another citable primary artifact. Then add a
separately validated sphere-built, gravity-loaded wall implementation and run
the complete campaign. Do not infer these values from DIRT output or tune them
to the published curves.

The paper metadata is checked through Crossref; the method claims are checked
against the hash-bound local PDF. Approximate Fig. 6/7 values in
`data/guo_2019_rubber_cord.csv` are AI-assisted manual digitizations and are
not raw data. Codex AI assisted this audit and its documentation. Independent
human verification of the primary pages and any future digitization remains
necessary.
