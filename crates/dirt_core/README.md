# dirt_core

The batteries-included umbrella crate for **DIRT** — a plugin-based Discrete
Element Method (DEM) simulation framework. Depend on this crate to get
everything.

## What It Does

`dirt_core` is the entry point to DIRT. It re-exports the `dirt_*` physics
crates and the `soil_*` substrate crates, and provides convenience plugin
groups so a complete simulation can be set up in a few lines.

The App / Plugin / dependency-injection scheduler / TOML / multi-stage run
substrate lives in the [grass](https://github.com/SueHeir/grass) framework; the
method-agnostic particle substrate (base `Atom`, `AtomData` registry, domain
decomposition, communication, neighbor lists) lives in
[soil](https://github.com/SueHeir/soil). DIRT adds the DEM physics on top.

## Plugin Groups

| Plugin Group | Purpose |
|---|---|
| `CorePlugins` | Infrastructure: input/TOML config, MPI or single-process communication, domain decomposition, neighbor lists, groups, run loop, and output. Does **not** include time integration. |
| `GranularDefaultPlugins` (from `dirt_granular`) | DEM granular: Hertz–Mindlin contacts, rotational dynamics, particle insertion, and Velocity Verlet integration |

Combine `CorePlugins` with `GranularDefaultPlugins`, or add individual plugins
for finer control.

## Quick Start

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
       .add_plugins(GranularDefaultPlugins);
    app.start();
}
```

Parameters are read from a TOML config file passed on the command line, e.g.
`cargo run -- config.toml`.

## The Prelude

Use `dirt_core::prelude::*` to import the most common types: the `CorePlugins`
group, the framework (`App`, `Plugin`, `PluginGroup`), core types (`Atom`,
`Config`, `RunState`), and the major DEM plugins (`DemAtomPlugin`,
`GranularDefaultPlugins`, `DemBondPlugin`, `ClumpPlugin`, `WallPlugin`,
`ContactAnalysisPlugin`, `MeasurePlanePlugin`, `FixesPlugin`,
`DeformPlugin`).

## Crate Organization

**Framework** ([`grass_*`](https://github.com/SueHeir/grass)): App / Plugin /
scheduler, TOML loading, multi-stage run loop, MPI abstraction, derive macros

**Substrate** ([`soil_*`](https://github.com/SueHeir/soil)): base `Atom` +
`AtomData` registry, domain decomposition, neighbor lists, Velocity Verlet
integration, output, box deformation

**DEM** (`dirt_*`) re-exported here:

| Crate | Description |
|---|---|
| `dirt_atom` | DEM per-atom data, material table, particle insertion, size distributions |
| `dirt_granular` | Hertz normal, Mindlin tangential, rotational dynamics, granular temperature |
| `dirt_bond` | Inter-particle bonds: normal/tangential/bending, auto-bonding, breakage |
| `dirt_clump` | Rigid clump (multisphere) composites for non-spherical particles |
| `dirt_wall` | Walls: plane, cylinder, sphere, cone, region surfaces; wall motion |
| `dirt_contact_analysis` | Coordination number, fabric tensor, rattler detection, per-contact CSV |
| `dirt_measure_plane` | Measurement planes for flux and profile sampling |
| `dirt_fixes` | General-purpose fixes: gravity, addforce, setforce, freeze, movelinear, viscous |

## Feature Flags

| Feature | Default | Description |
|---|---|---|
| `mpi_backend` | yes | MPI-based parallel communication. Disable with `--no-default-features` for single-process runs. |

## License

MIT OR Apache-2.0
