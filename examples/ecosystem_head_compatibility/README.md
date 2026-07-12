# GRASS → SOIL → DIRT HEAD compatibility matrix

This validation exercises the same temporary Cargo path patches used by
[`ci/ecosystem-head-check.sh`](../../ci/ecosystem-head-check.sh). It tests the
actual source trees rather than the unrelated Git revisions in `Cargo.lock`.

The matrix has two required outcomes. The reviewed GRASS/SOIL/DIRT tuple must
complete `cargo metadata` and the non-MPI `precision-double` workspace check.
Then the same DIRT and GRASS revisions are paired with SOIL `72622c1`, whose
typed particle-extension update requires `AtomData::snapshot` in both
`dirt_granular` and `dirt_bond`; that tuple must fail and show the API drift.
The second failure is a successful detection result, not a compatibility pass.

![GRASS SOIL DIRT compatibility matrix](plots/ecosystem_head_compatibility_matrix.png)

*Actual compatibility matrix: the reviewed tuple passes metadata and the
non-MPI precision-double check; the newer SOIL tuple visibly fails on the two
missing `AtomData::snapshot` implementations. Latest result: PASS—the gate
detects drift instead of silently resolving locked Git sources.*

## Reproduce

Run the matrix (it creates and removes detached worktrees at the revisions in
`config.toml`):

```bash
source ~/projects/.build-env
$BENCH_PYTHON examples/ecosystem_head_compatibility/sweep.py
```

For the developer/CI gate against checkout heads, use the one command documented
in the repository root:

```bash
ci/ecosystem-head-check.sh --grass ../grass --soil ../soil
```
