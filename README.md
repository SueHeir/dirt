# DIRT — Discrete Interaction-Resolved Toolkit

Discrete Element Method (DEM) physics in the **GRASS → SOIL → DIRT** stack.

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists
       └─ DIRT   DEM physics: contact, bonds, walls, clumps, thermal   ← you are here
```

DIRT resolves every inter-particle contact individually (vs continuum models),
riding the [SOIL](https://github.com/SueHeir/soil) substrate for all the
method-agnostic machinery (atom data, domain decomposition, halo exchange,
neighbor lists). It adds the granular physics: Hertz–Mindlin contact, rotational
dynamics, parallel bonds, walls, multisphere clumps, heat conduction, and
contact analysis. Target applications: fiber/yarn/weave, regolith, and
fluid–solid coupling.

## Layout

`dirt_core` is the batteries-included umbrella crate (depend on this):

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins).add_plugins(GranularDefaultPlugins);
    app.start();
}
```

| crate | role |
|---|---|
| `dirt_core` | umbrella: `CorePlugins`, `GranularDefaultPlugins`, prelude |
| `dirt_atom` | per-atom DEM data (`DemAtom`), materials, particle insertion |
| `dirt_granular` | Hertz/Mindlin contact, rolling/twisting, adhesion, rotational dynamics |
| `dirt_wall` | plane/cylinder/sphere/cone/region-surface walls |
| `dirt_bond` | bonded-particle model: normal/shear/twist/bending beam, breakage, plasticity |
| `dirt_clump` | multisphere/clump rigid composites |
| `dirt_thermal` | contact heat conduction |
| `dirt_contact_analysis` | coordination number, fabric tensor, rattlers |
| `dirt_measure_plane` | measurement planes for flux/profiles |
| `dirt_fixes` | group fixes: addforce, setforce, move_linear, freeze, pin, viscous, nve_limit, gravity |
| `dirt_test_utils` | shared test helpers |

## Building

DIRT path-depends on [soil](https://github.com/SueHeir/soil) and
[grass](https://github.com/SueHeir/grass) as siblings; clone all three side by
side:

```
GitHub/
  grass/
  soil/
  dirt/   ← here
```

```bash
cargo run --release --example hopper --no-default-features -- examples/hopper/config.toml
```

## License

MIT OR Apache-2.0
