# Introduction

**DIRT** — the Discrete-element Interaction-Resolved Toolkit — is a Discrete
Element Method (DEM) physics engine written in Rust. It resolves every
inter-particle contact individually: Hertz–Mindlin contact, rotational
dynamics, parallel bonds, walls, multisphere clumps, heat conduction, and
contact analysis.

DIRT is the top tier of a three-repo stack. Each tier depends only on the ones
below it:

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling      (no particles)
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists   (no physics)
       └─ DIRT   DEM physics: contact, bonds, walls, clumps   ← you are here
```

- **[GRASS](https://github.com/SueHeir/grass)** is the framework: a Bevy-style
  `App` + `Plugin` container, a dependency-injection scheduler, TOML config, and
  MPI coupling. It knows nothing about particles.
- **[SOIL](https://github.com/SueHeir/soil)** is the substrate: the base `Atom`,
  domain decomposition, halo/ghost communication, atom migration, and neighbor
  lists. It knows nothing about physics.
- **DIRT** adds the granular physics on top.

If you want to *use* a granular DEM engine — fill a hopper, measure an angle of
repose, validate a Hertz rebound — start here. If you want to *write your own*
particle physics (SPH, peridynamics, your own force law), the substrate and
framework books are where that story lives:

- The [SOIL book](https://sueheir.github.io/soil) — how to add per-particle
  state and ride the substrate's communication for free.
- The [GRASS book](https://sueheir.github.io/grass) — how to build a
  plugin-based solver from scratch.

## What DIRT gives you

`dirt_core` is the batteries-included umbrella crate. Depend on it and you get
the whole stack re-exported, plus convenience plugin groups:

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)            // I/O, comm, domain, neighbors, run loop
       .add_plugins(GranularDefaultPlugins); // Hertz–Mindlin contact + Velocity Verlet
    app.start();
}
```

Everything else — walls, bonds, clumps, gravity, measurement planes — is an
additional plugin you opt into. The [DEM Physics](./physics/contact.md) chapters
cover each in turn.

## How to read this book

- **[Getting Started](./getting-started/installation.md)** — build the engine and
  run your first simulation in a few minutes.
- **[The Stack](./stack/overview.md)** — the architecture, and why DIRT is split
  across three repos.
- **[DEM Physics](./physics/contact.md)** — the force laws and how to configure
  them.
- **[Reference](./reference/config.md)** — the full config schema and the
  validation suite.
