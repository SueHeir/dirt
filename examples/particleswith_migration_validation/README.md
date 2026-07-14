# `particleswith_migration_validation`

This compatibility check runs representative contact/wall, BPM bond, and clump
paths from scratch in a detached `41cb2f1` baseline worktree and the current
typed-`ParticlesWith` checkout. It compares each emitted CSV by SHA-256 and its
terminal observable, so a regression in either revision is visible to the
runner.

Run `$BENCH_PYTHON examples/particleswith_migration_validation/sweep.py` after
`source ~/projects/.build-env`. The pinned Python plotting dependency is in
`requirements.txt`; install it with `$BENCH_PYTHON -m pip install -r
examples/particleswith_migration_validation/requirements.txt` if needed. The
runner regenerates the CSV and figure and fails if either the full outputs or
the strict `1e-15` comparison changes.

![Before/after compatibility](plots/before_after_compatibility.png)

*PASS: all three representative outputs are byte-identical before and after the
query migration; the red markers show the `1e-15` pass limit.*
