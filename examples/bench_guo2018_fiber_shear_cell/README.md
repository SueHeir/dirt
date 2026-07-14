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
python3 examples/bench_guo2018_fiber_shear_cell/figure_scale_audit.py --verify
python3 examples/bench_guo2018_fiber_shear_cell/figure_scale_audit.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/campaign_preflight.py
python3 examples/bench_guo2018_fiber_shear_cell/campaign_preflight.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/result_evidence_contract.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/reconstruction_readiness.py --verify
python3 examples/bench_guo2018_fiber_shear_cell/normal_load_audit.py --self-test
python3 examples/bench_guo2018_fiber_shear_cell/status_contract.py --verify
python3 examples/bench_guo2018_fiber_shear_cell/status_contract.py --self-test
```

The hash-bound receipt in `data/reference_provenance.json`, the page-receipted
`data/reconstruction_ledger.json`, and the candidate-only topology generator
make the positive evidence auditable. The 64-mm candidate can be generated and
audited without invoking a solver:

`figure_scale_audit.py` separately records an important negative source fact:
Figure 2 is a perspective rendering with no scale bar, camera calibration, or
wall-sphere count. It confirms the sphere-built boundary class but cannot turn
a pixel measurement into an exact wall diameter or lattice.

```bash
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py \
  --width-mm 64 --candidate-only --output /tmp/guo64-candidate
python3 examples/bench_guo2018_fiber_shear_cell/prepare.py --audit /tmp/guo64-candidate
```

`AUDIT OK` confirms only Table-2/cell arithmetic and topology counts; it is
not a physics result. The 96-mm/750-fibre population is explicitly marked as a
derived sensitivity candidate, not a reported source input.

## Sphere-built boundary candidates

`wall_realisation.py` now materializes three predeclared, non-calibrated
sphere-built wall/blade meshes (0.6, 1.2, and 2.4 mm). They enforce the paper's
64 x 36-mm periodic planform, 8-mm blade pitch, and 2/4-mm lower/upper blade
lengths, while preserving the unresolved sphere diameter and lattice as a
non-equivalence. For example:

```bash
python3 examples/bench_guo2018_fiber_shear_cell/wall_realisation.py \
  --diameter-mm 1.2 --output /tmp/guo-wall-1p2
python3 examples/bench_guo2018_fiber_shear_cell/source_geometry_audit.py \
  --verify --wall-realisation /tmp/guo-wall-1p2/wall_realisation.json
```

These artifacts are generated geometry and provenance only: they do not drive
the current DIRT wall API, produce a history, or support a Fig. 6/7 comparison.
They make the missing implementation boundary precise: a future solver must
consume every predeclared candidate through rigid moving assemblies, retain all
histories, and report sensitivity rather than choose a target-fitting mesh.
Any wall selected from an archival artifact or independent measurement must
also provide a locator and SHA-256 receipt; the geometry audit rejects an
unsupported evidence label.

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

## Boundary uncertainty, without curve fitting

The missing wall-sphere diameter and tessellation make an exact source mesh
unidentifiable, but they do not make the reported cell physics unknowable. A
future campaign may therefore use a separately evidenced wall mesh, or a
predeclared resolution-sensitivity ensemble. It must label that campaign as a
non-equivalent reconstruction, retain all candidate runs, and select no mesh
using the Figure 6/7 targets. `source_geometry_audit.py --verify` enforces the
provenance rule for any proposed mesh. This is intentionally not a route to a
replication PASS: the current tree has no such mesh, no solver history, and no
DIRT-vs-reference result.

`data/non_equivalent_sensitivity_campaign.json` now makes the permissible
fallback concrete without pretending it resolves the source omission: it
pre-registers three distinct sphere-wall resolutions, all three digitized
normal-stress cases, both eventual observables, and a fixed strain window.
`campaign_preflight.py --require-results RESULT_ROOT` refuses any later result
claim unless each wall/load case supplies its solver input, history, observable
summary, retained wall manifest, and a SHA-256 receipt that identifies the
solver revision and discloses AI authorship/limits. It deliberately does not compare a result to
Fig. 6/7; that comparison belongs only after the complete, non-selected
ensemble exists. The values are a sensitivity bracket, **not** reported wall
dimensions, not calibration parameters, and not a replication result.

The receipt check makes artifacts tamper-evident and attributable; it does not
validate numerical correctness, convergence, or physical agreement. Those
remain independent review tasks and cannot be satisfied by internally
consistent files or by AI-generated assertions.

`data/replication_status.json` records the scientific status in a
machine-readable form. `status_contract.py` cross-checks it against the
independent source ledger and geometry contract, preventing a solver or
curve-comparison claim while the source wall geometry remains unresolved.
