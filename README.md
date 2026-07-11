# DIRT

<!-- disclaimer-banner -->
> **Research-software status:** This ecosystem is AI-authored and under active evaluation. DIRT's DEM claims are accompanied by reproducible analytical, cross-code, or empirical evidence; repositories prefixed `dev_` are experimental method demonstrations outside the author's domain expertise. See [DISCLAIMER.md](DISCLAIMER.md) and [examples/VALIDATION.md](examples/VALIDATION.md).
<!-- /disclaimer-banner -->


**A granular-DEM solver that runs as a complete code today, built on a
composition-oriented simulation stack.**

DIRT — the *Discrete-element Interaction-Resolved Toolkit* — provides contact,
walls, bonds, clumps, integration, diagnostics, and validated DEM examples. Its
unusual property is not simply that those capabilities are Rust plugins. It is
that DIRT shares its state, lifecycle, and schedule model with the rest of the
GRASS ecosystem.

Run DIRT as a standalone DEM application, and add or replace its DEM physics as
ordinary systems. The same underlying GRASS and SOIL composition model is a
foundation for future multi-solver work, but DIRT does not yet ship or validate a
cross-substrate coupling case. Such a case belongs in a dedicated
`dev_couple_*` repository, where its exchange and validation can be assessed on
their own terms.

DIRT is a ground-up Rust reimplementation building on roughly two years of prior
DEM development (formerly MDDEM).

## A simulation is a few plugins you assemble in Rust

You build a DIRT program by composing plugins from the crate prelude. Here is a
complete settling-bed solver:

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)            // I/O, comm, domain, neighbors, run loop
       .add_plugins(GranularDefaultPlugins) // Hertz–Mindlin contact + Velocity Verlet
       .add_plugins(GravityPlugin)          // opt into the physics you want …
       .add_plugins(WallPlugin);            // … one plugin at a time
    app.start();
}
```

That is the whole point of the stack: **a plugin is just a set of systems, and a
system is just a function.** When you need physics that isn't in the box, you
write it as one more function and `add_update_system` it — the scheduler injects
the typed state it asks for (`Res<Atom>`, `ResMut<Walls>`, …) and the substrate
handles decomposition, so the same code runs serial or across MPI ranks. Learn
these seams once and you can write the physics *you* need — and the same mental
model carries to any other solver on the stack. The
[Your First Simulation](docs/src/getting-started/first-simulation.md) walk-through
adds a real custom system (remove a blocker wall when the bed goes quiet) in a
dozen lines.

## Why this architecture matters for DEM

DEM is frequently one half of a larger problem: particles interact with a fluid,
deformable structure, thermal field, electrostatic field, or continuum model. If
the DEM code owns a private timestep driver and private state model, every
coupling begins as an integration project.

DIRT is organized around a common composition layer:

- SOIL owns reusable parallel particle infrastructure.
- DIRT owns DEM physics.
- GRASS owns scheduling, lifecycle, I/O, and solver composition.
- A future dedicated coupling can own the exchange between participating solvers.

Validation remains essential evidence for individual DEM claims, but it is not
evidence of an unbuilt coupling capability. The architecture keeps the relevant
solver seams explicit; whether they support a particular coupling remains to be
demonstrated in a dedicated, validated case.

> Prefer config files to code? A prebuilt driver runs the shipped scenarios and
> parameter sweeps from TOML with no recompile — see
> [Run scenarios and sweeps from config](#run-scenarios-and-sweeps-from-config-no-recompile).
> It's a convenience for canned runs, not the main path.

## Trust: what's actually validated

DIRT's first-class feature is that it tells you where to believe it. It ships a
suite of `bench_*` validation examples under [`examples/`](examples/); each
couples a small simulation to an **independent** reference and checks measured
quantities against it with explicit tolerances. Every benchmark's `sweep.py graph`
step prints a **PASS/FAIL** verdict, so the suite is a regression net — not a
gallery of runs that happen to pass.

Each benchmark is graded by how strong its evidence is:

- **Analytical** — agreement with a closed-form result (Hertz contact duration,
  Euler–Bernoulli beam deflection, Haff's `T_g ∝ t⁻²` cooling).
- **Cross-code** — agreement with **LAMMPS** on the same problem. This tests
  implementation consistency under a *shared* contact model, not correctness
  against physical reality.
- **Empirical / law / qualitative** — agreement with a scaling exponent or a
  fitted correlation (Beverloo discharge, angle of repose vs. friction).

And each one says plainly where it is weak — an idealization, an empirical fit, a
check that is really self-consistent (a model returning its own input), or a
regime the run never reaches. Some benchmarks are still **red on purpose**
(`bench_column_collapse` FAILs its exponent gate today), and that's left visible
rather than tuned away. There is no comparison to raw experimental data in the
suite; the closest tie is `bench_kharaz_oblique`, anchored to Kharaz et al.'s
measured restitution and friction.

**These tests catch real bugs.** The oblique-impact validation alone drove two
contact-model fixes — a tangential damping-sign error that was injecting energy,
and a requirement that a frozen contact partner also freeze its rotation — and
the rebound benchmark surfaced a mislabeled damping constant.

The authoritative, continuously-updated write-up of every figure and every weak
spot lives in [`examples/VALIDATION.md`](examples/VALIDATION.md).

## The physics menu (opt in per plugin)

Everything past the contact force is a plugin you add only if you want it:

- **Contact** — Hertz–Mindlin normal + tangential, rolling/twisting resistance,
  JKR/DMT-SJKR adhesion, rotational dynamics.
- **Walls** — plane / cylinder / sphere / cone / region-surface, with Mindlin
  wall friction and servo-controlled boundaries.
- **Bonds** — bonded-particle beams (normal/shear/twist/bending), breakage,
  plasticity.
- **Clumps** — rigid multisphere composites.
- **Fixes** — gravity, prescribed motion, pin/freeze, viscous damping,
  add/set-force.
- **Diagnostics** — coordination number, fabric tensor, rattlers, measurement
  planes for flux and profiles.

## Run scenarios and sweeps from config (no recompile)

`examples/run` is a prebuilt driver that assembles the full plugin stack —
contact, gravity, walls, fixes, box deformation — and reads every case-specific
detail from a **declarative TOML config**, so you can run the shipped scenarios
(settle, pour, cylinder pour, Lees–Edwards shear, uniaxial compression) without
writing or recompiling any Rust:

```bash
cargo run --release --example run -- examples/run/pour_settle.toml
```

Two things this prebuilt driver gives you that a hand-written `main.rs` does not:

- **Toggle physics without recompiling.** A hand-written program fixes its plugin
  set at compile time; to add walls or deformation you edit Rust and rebuild. The
  driver pre-adds the superset, and each plugin stays inert unless its config
  section is present — so settle → pour → shear are chosen purely by the TOML.
- **Parameter sweeps are one binary + N configs.** Point the same binary at a
  series of TOMLs varying volume fraction, friction, or strain rate and nothing
  recompiles. (This is the pattern the validation suite runs — each `bench_*`
  builds one example binary, then loops it over a config tree.) For a deformation
  run, the declarative `[loading]` block lets the driver derive the step count
  and own the deform loop for you.

The config is declarative throughout: you *describe* geometry, materials,
insertion, walls, body forces, and loading; you never *script* a step sequence.
Reach for a custom plugin instead when you need physics or a measurement the
shipped stack doesn't already have. Details:
[Run from a Config](docs/src/getting-started/run-from-config.md).

## The stack

DIRT is the top tier of a three-repo stack; each tier depends only on the ones
below it and knows nothing about the tiers above:

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling      (no particles)
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists   (no physics)
       └─ DIRT   DEM physics: contact, bonds, walls, clumps   ← you are here
```

## How the three fit together

GRASS gives you the `App`/scheduler/coupling; SOIL turns that into a parallel
particle substrate via one `AtomData` contract; DIRT is the proof that a full
LAMMPS-validated physics tier rides it — and the same seams are open for SPH,
peridynamics, or your own method.

One line per tier, worded identically wherever these three repos describe
themselves:

- **[GRASS](https://github.com/SueHeir/grass)** — Build solvers as composable
  plugins instead of a hand-rolled main loop — explicit time-stepping or a
  single implicit global solve, particles or a mesh — and couple several
  together, in-process or across MPI.
- **[SOIL](https://github.com/SueHeir/soil)** — Write your own particle method
  without hand-writing domain decomposition, halo exchange, migration, and
  neighbor lists — declare your state once, SOIL carries it through all of it.
- **[DIRT](https://github.com/SueHeir/dirt)** — A LAMMPS-validated granular-DEM
  engine, easily extended by composing Rust plugins on the GRASS framework.

**Where to start:** to *run* granular simulations, start at
[DIRT](https://github.com/SueHeir/dirt), the batteries-included physics tier; to
*write your own* particle method or solver, start at
[SOIL](https://github.com/SueHeir/soil) (the particle substrate) or
[GRASS](https://github.com/SueHeir/grass) (the framework). The full walkthrough
of how the tiers compose — one timestep end to end, and where the seams are — is
the canonical [How the stack fits together](https://sueheir.github.io/grass/stack/how-the-stack-fits-together.html)
page in the GRASS book.

## Install

You need **Rust** (stable, 2021 edition or newer; [rustup.rs](https://rustup.rs)).
You do *not* need to check out GRASS or SOIL — DIRT pulls them from GitHub during
the build.

```bash
git clone https://github.com/SueHeir/dirt
cd dirt
cargo run --release --example hello_bed \
  --no-default-features --features precision-double \
  -- examples/hello_bed/config.toml
```

The default feature set is `["mpi_backend", "precision-double"]`.
`--no-default-features` turns off `mpi_backend` and builds a single-process
binary — the fastest way to get running, with no C compiler or MPI library
required — but it also drops `precision-double`, which the solver requires, so
you re-add it explicitly with `--features precision-double`. Drop
`--no-default-features` entirely once you have an MPI toolchain and want
multi-rank domain-decomposed runs:

```bash
cargo build --release           # mpi_backend on by default
mpirun -np 4 ./target/release/examples/hopper examples/hopper/config.toml
```

Full details, including using DIRT as a library dependency, are in
[Installation & Building](docs/src/getting-started/installation.md).

## Crate map (reference)

`dirt_core` is the batteries-included umbrella crate — depend on that and you get
the prelude plus the plugin groups. The rest are the individual physics tiers,
useful when you want to reach for one directly:

| crate | role |
|---|---|
| [`dirt_core`](crates/dirt_core/README.md) | umbrella: `CorePlugins`, `GranularDefaultPlugins`, prelude |
| [`dirt_atom`](crates/dirt_atom/README.md) | per-atom DEM data (`DemAtom`), materials, particle insertion |
| [`dirt_granular`](crates/dirt_granular/README.md) | Hertz/Mindlin contact, rolling/twisting, adhesion, rotational dynamics |
| [`dirt_wall`](crates/dirt_wall/README.md) | plane/cylinder/sphere/cone/region-surface walls, with Mindlin wall friction |
| [`dirt_bond`](crates/dirt_bond/README.md) | bonded-particle model: normal/shear/twist/bending beam, breakage, plasticity |
| [`dirt_clump`](crates/dirt_clump/README.md) | multisphere/clump rigid composites |
| [`dirt_contact_analysis`](crates/dirt_contact_analysis/README.md) | coordination number, fabric tensor, rattlers |
| [`dirt_measure_plane`](crates/dirt_measure_plane/README.md) | measurement planes for flux/profiles |
| [`dirt_fixes`](crates/dirt_fixes/README.md) | DEM group fixes: add/set force, freeze, pin, prescribed motion, viscous damping, gravity |
| [`dirt_test_utils`](crates/dirt_test_utils/README.md) | shared test helpers |

## Citing

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff) (GitHub renders a
"Cite this repository" button from it); please cite the version you ran, with
per-version changes in [`CHANGELOG.md`](CHANGELOG.md). A JOSS paper is planned.

## License

MIT OR Apache-2.0
