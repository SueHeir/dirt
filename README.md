# DIRT

**A Rust granular-DEM engine that resolves every contact individually —
cross-checked against LAMMPS and closed-form theory.**

DIRT — the *Discrete-element Interaction-Resolved Toolkit* — simulates granular
matter by computing each inter-particle contact force directly: Hertz–Mindlin
contact, rotational dynamics, parallel bonds, walls, multisphere clumps, heat
conduction, and contact analysis. It rides the [SOIL](https://github.com/SueHeir/soil)
substrate for the method-agnostic machinery (atom data, domain decomposition,
halo exchange, neighbor lists) and adds only the granular physics. It is a
ground-up Rust reimplementation building on roughly two years of prior DEM
development (formerly MDDEM).

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

> Prefer config files to code? There's a zero-Rust runner too — see
> [below](#dont-want-to-write-rust). It's a convenience, not the main path.

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
- **Heat** — granular conduction through contacts.

## Don't want to write Rust?

If you'd rather describe a run than compile one, DIRT ships a single generic
`run` driver that assembles the standard plugin stack and takes every
case-specific detail from a **TOML config** you pass on the command line:

```bash
cargo run --release --example run -- examples/run/pour_settle.toml
```

The config is **declarative** — you *describe* geometry, materials, insertion,
walls, body forces, and duration; you never *script* a step sequence (the driver
owns the loop). Swap the config path to run a different scenario with the same
binary. Details: [Run from a Config — Zero Rust](docs/src/getting-started/run-from-config.md).

## The stack

DIRT is the top tier of a three-repo stack; each tier depends only on the ones
below it and knows nothing about the tiers above:

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling      (no particles)
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists   (no physics)
       └─ DIRT   DEM physics: contact, bonds, walls, clumps   ← you are here
```

GRASS gives you the `App`/scheduler/coupling; SOIL turns that into a parallel
particle substrate via one `AtomData` contract; DIRT is the proof that a full
LAMMPS-validated physics tier rides it — and the same seams are open for SPH,
peridynamics, or your own method. See the
[GRASS](https://github.com/SueHeir/grass) and [SOIL](https://github.com/SueHeir/soil)
books to write your own tier.

## Install

You need **Rust** (stable, 2021 edition or newer; [rustup.rs](https://rustup.rs)).
You do *not* need to check out GRASS or SOIL — DIRT pulls them from GitHub during
the build.

```bash
git clone https://github.com/SueHeir/dirt
cd dirt
cargo run --release --example hello_bed --no-default-features -- examples/hello_bed/config.toml
```

`--no-default-features` disables the `mpi_backend` feature and builds a
single-process binary — the fastest way to get running, with no C compiler or
MPI library required. Drop it once you have an MPI toolchain and want multi-rank
domain-decomposed runs:

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

## How to cite

If you use DIRT in academic work, please cite the version you used so results
stay reproducible. Machine-readable metadata lives in
[`CITATION.cff`](CITATION.cff) (GitHub renders a "Cite this repository" button
from it); per-version changes are in [`CHANGELOG.md`](CHANGELOG.md).

```bibtex
@software{suehr_dirt_2026,
  author  = {Suehr, Elizabeth},
  title   = {{DIRT — Discrete-element Interaction-Resolved Toolkit}},
  version = {0.1.3},
  year    = {2026},
  url     = {https://github.com/SueHeir/dirt},
  license = {MIT OR Apache-2.0}
}
```

Update `version` (and `year`) to match the release you actually ran.

## License

MIT OR Apache-2.0
