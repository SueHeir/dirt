# bench_mpi_decomposition — MPI cross-rank DEM correctness

Closes the long-standing validation gap noted in `examples/VALIDATION.md`:
every other benchmark runs `1×1×1`, so DIRT's MPI domain decomposition —
ghost exchange, atom migration across rank boundaries, and per-contact history
transfer — was **untested at the dirt/DEM level**. This bench proves that a
contact-rich run under a `>1` processor decomposition reproduces the single-rank
result.

## What it proves

The same seeded, fully-periodic, gravity-off frictional granular gas (~400 glass
spheres, Hertz–Mindlin contact, restitution 0.9) is run at three decompositions
of the **same** box:

| tag  | grid    | ranks | boundaries crossed by migrating atoms |
|------|---------|-------|---------------------------------------|
| `n1` | `1×1×1` | 1     | none (serial reference)               |
| `n2` | `2×1×1` | 2     | `x = 0.01`                            |
| `n4` | `2×2×1` | 4     | `x = 0.01` and `y = 0.01`             |

Each multi-rank run is asserted against the `1×1×1` reference, sampled from the
MPI-gathered `[dump]` frames (each frame gathers every rank's *local* atoms to
rank 0, so a lost or duplicated atom is caught as a wrong global count / tag set):

- **atom-count / identity conservation** — every gathered frame holds exactly
  `N = 400` atoms with the exact reference tag set (migration + ghost exchange
  never drop or duplicate an atom);
- **global momentum conservation** — with gravity off and a periodic box, total
  linear momentum `P = Σ mᵥ` is an exact invariant; `|P(t) − P(0)|` stays at
  round-off for every decomposition, and the multi-rank `P(t)` matches the
  reference (a dropped ghost-force reverse-comm would drift it);
- **global kinetic-energy agreement** — `KE(t)` of each multi-rank run matches
  the `1×1×1` trajectory at every sample;
- **per-atom trajectory agreement** — at the final step, per-atom velocities and
  minimum-image positions (matched by global tag) agree with `1×1×1` to the
  floating-point associativity floor.

## Why the agreement is the FP floor, not a loosened band

`[neighbor] every = 1`, `check = false`, `sort_every = 0` pins the neighbor
rebuild schedule and disables the cache re-sort, so the **only** thing that
differs between the runs is the domain decomposition — i.e. the order in which
pairwise and reverse-communicated ghost forces are reduced. That difference is
pure floating-point associativity. Measured worst case over 4000 steps:

```
per-atom  : pos 6.0e-17 , vel 7.8e-14
momentum  : drift 9.1e-23 , cross-decomposition match 4.2e-23
KE        : cross-decomposition match 7.5e-16
```

![MPI decomposition deltas](plots/mpi_decomposition_deltas.svg)

*Measured `2×1×1` and `2×2×1` tag/count, momentum, energy, and final per-atom
state deltas against the `1×1×1` serial reference. The dashed line is the
unchanged `1e-9` pass tolerance; the tag/count identity delta is exactly zero.
PASS.*

— four or more orders of magnitude under the `1e-9` acceptance gate (the same
FP-floor band `bench_restart_determinism` uses). The gate is not a relaxed
physics tolerance; the runs agree essentially to machine epsilon.

## Run

```bash
# Standalone single-rank smoke test:
cargo run --release --example bench_mpi_decomposition -- \
    examples/bench_mpi_decomposition/config.toml

# 2-rank MPI run along x:
cargo build --release --example bench_mpi_decomposition
mpiexec -n 2 target/release/examples/bench_mpi_decomposition \
    examples/bench_mpi_decomposition/config.toml

# Full gated benchmark (1×1×1 vs 2×1×1 vs 2×2×1, compares + gates):
python3 examples/bench_mpi_decomposition/sweep.py
```

`sweep.py` exits `0` iff every check passes and prints `ALL CHECKS PASSED`.

The binary is the stock plugin wiring (`CorePlugins` + `GranularDefaultPlugins`)
— no bench-specific code — and is built with the crate's default features
(`mpi_backend` + `precision-double`), so the same binary serves the serial
reference and the `mpiexec` runs. Everything that varies between runs (the
`[comm]` processor grid, step count, dump cadence, output dir) is composed from
the declarative `config.toml` by the driver.
