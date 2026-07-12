# Clump Insertion Determinism

Validates the production `[[clump.insert]]` setup path, including
`ClumpPlugin` and `clump_insert_atoms`. Two runs with the same seed must write
identical inserted atom/body state; a changed seed must alter the state.

Run:

```bash
python3 examples/bench_clump_insertion_determinism/sweep.py
```

The gate byte-compares the full CSV fingerprints for same-seed runs and measures
the maximum numeric divergence for the changed-seed run.

![Clump insertion determinism](plots/clump_insertion_determinism.svg)

Same-seed config-level insertion is byte-identical; changing the seed produces a
non-zero fingerprint divergence. PASS means the actual setup path is seeded.
The graph is generated with Python's standard library (no matplotlib or other
undeclared runtime dependency).
