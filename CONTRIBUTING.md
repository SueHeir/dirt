# Contributing to DIRT

Thanks for your interest in DIRT. This is a small research project, and
contributions — a bug report, a failing test case, a new benchmark, a physics
plugin — are genuinely welcome. This document tells you how to report a problem,
where to get help, and how to land a code change.

By participating you agree to the [Code of Conduct](CODE_OF_CONDUCT.md).

## Where DIRT is developed

DIRT's canonical development happens on a **private Gitea server**, where every
change goes through a pull request, a reviewer, and a regression gate (the
`bench_*` validation suite) before it reaches `main`. Stable `main` is then
**published downstream to the public GitHub mirror** at
<https://github.com/SueHeir/dirt>. GitHub is where the released code lives and
where the outside world interacts with the project.

The practical consequence: **you contribute through GitHub.** Open issues and
pull requests there. A maintainer replays accepted work onto the internal Gitea,
where it runs the same review-and-regression flow before being published back
out. You don't need Gitea access to contribute — the GitHub side is the door.

## Reporting a bug or issue

Open an issue at <https://github.com/SueHeir/dirt/issues>. Before filing, a quick
search for an existing report saves everyone time.

A good bug report gets a fast fix. Please include:

- **What you ran** — the example or config, the exact `cargo run` command, and
  whether it was serial (`--no-default-features`) or multi-rank (`mpirun -np N`).
- **What you expected vs. what happened** — numbers, a panic message, or a
  `PASS`/`FAIL` line from a benchmark's `sweep.py graph` step.
- **Environment** — Rust version (`rustc --version`), OS, and MPI toolchain if
  the build used the `mpi_backend` feature.
- **A minimal reproducer** if you can manage one — the smallest config or
  particle count that still shows the problem.

Physics disagreements are bugs too. If a benchmark's measured quantity drifts
from its reference outside the stated tolerance, that's worth an issue even if
nothing crashed. (Some benchmarks are red *on purpose* — see
[`examples/VALIDATION.md`](examples/VALIDATION.md) before reporting one.)

## Getting help / support

- **Usage questions** ("how do I set up walls?", "which plugin does X?") —
  open a [GitHub Discussion](https://github.com/SueHeir/dirt/discussions) if the
  repository has them enabled, otherwise a GitHub issue labelled `question`.
- **Start with the docs.** The [DIRT book](docs/src/introduction.md) covers
  getting started, the physics plugins, the config runner, and the validation
  reference. [`examples/`](examples/) is full of runnable, commented programs;
  the `bench_*` ones double as worked references.
- **Private matters** (security-sensitive reports, Code of Conduct concerns) —
  email Elizabeth Suehr at <elizabeth.suehr@gmail.com> rather than opening a
  public issue.

## Contributing code

DIRT is a library stack. Most new physics is a **plugin** — a set of systems you
`add_update_system` — and most of it belongs in an **example**, not the library
core. Keep that split in mind and a change reviews quickly.

### Workflow

1. **Fork and branch.** Fork on GitHub, branch off `main`.
2. **Make the change**, keeping the tier boundaries below intact.
3. **Validate it** (next section) and put the evidence in your PR description.
4. **Open a pull request** against `main` at
   <https://github.com/SueHeir/dirt/pulls>. Describe what changed and *why*, and
   paste the `PASS` output or the numbers that back it.

Commits use a `Signed-off-by`/co-author trailer as appropriate; keep them
focused and their messages in the imperative mood, matching the existing log.

### Where code lives — `src` vs. example

Default to putting example-, benchmark-, or test-specific logic **in the
example**, not in a crate's `src`. Local helpers in your example are cheap;
bloating the library core is expensive. Only add to `crates/*/src` when the thing
is **general, reusable across multiple examples, and simple**. A benchmark that
validates physics almost never needs new core API — reach for the example first.
A pull request that grows `src` to serve one example will be asked to move it.

### Respect the tier boundaries

DIRT is the top of a three-tier stack, and each tier stays agnostic to the ones
above it. Changes that reach below DIRT are held to those contracts:

- **DIRT** owns the DEM physics: contact, bonds, walls, clumps, heat, diagnostics.
- **[SOIL](https://github.com/SueHeir/soil)** is the particle substrate and must
  stay method-agnostic — no DEM-specific fields (torque, damage, bonds, Hertzian
  overlap) on the base atom; those are DIRT-side `AtomData` columns.
- **[GRASS](https://github.com/SueHeir/grass)** is the framework and knows nothing
  about particles at all.

If a change assumes particles or a specific force law, it belongs in DIRT, not
below it.

### Validate — and never weaken a check to go green

DIRT's whole trust story is that its benchmarks mean something, so validation is
not optional:

```bash
source ~/projects/.build-env   # sets the Rust + MPI env and BENCH_PYTHON
cargo build --release
cargo test --workspace
```

### Cross-tier HEAD compatibility

From a DIRT checkout, run this one command to check its current commit against
the clean local HEADs of GRASS and SOIL (the default locations are
`~/projects/grass` and `~/projects/soil`):

```bash
ci/ecosystem-head-check.sh
```

It prints each repository's commit and remote, patches Cargo only for that run,
then runs `cargo metadata` and the non-MPI `precision-double` workspace check.
Use `--grass PATH --soil PATH` to check isolated worktrees or another pair of
candidate commits. The Gitea CI workflow runs the same command after checking
out the three `main` heads.

If you touched physics, run the relevant benchmark and cite the verdict:

```bash
$BENCH_PYTHON examples/<bench_dir>/sweep.py run
$BENCH_PYTHON examples/<bench_dir>/sweep.py graph   # prints PASS / FAIL
```

**Do not make a check pass by weakening its validation** — don't loosen a
tolerance, shorten a fitted window, shrink the case set, or back-fit a reference
value. Reference values come from theory, a closed-form result, an independent
code (LAMMPS), or a citation. A green light that was bought by lowering the bar is
worse than an honest red one, and it will be sent back. If a benchmark is red for
a real reason, say so plainly in the PR — that is a valued kind of contribution,
not a failure.

New or changed benchmarks should keep documenting where they are weak — the
idealization, the empirical fit, the regime the run never reaches — in the same
honest voice as [`examples/VALIDATION.md`](examples/VALIDATION.md).

## License

By contributing, you agree that your contributions are licensed under the same
terms as DIRT: **MIT OR Apache-2.0** (see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE)).
