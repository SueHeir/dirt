# `particleswith_migration_validation`

This compatibility check compares the emitted CSVs from `origin/main`
(`41cb2f1`) and the typed `ParticlesWith` migration for representative
contact/wall, BPM bond, and clump paths. Each case was executed in a separate
isolated worktree with non-MPI `precision-double`; the full CSV SHA-256 hashes
match, while the plotted observables show the actual before/after difference.

Run `python3 examples/particleswith_migration_validation/sweep.py` to regenerate
the figure and fail if either the full outputs or the strict `1e-15` comparison
changes.

![Before/after compatibility](plots/before_after_compatibility.png)

*PASS: all three representative outputs are byte-identical before and after the
query migration; the red markers show the `1e-15` pass limit.*
