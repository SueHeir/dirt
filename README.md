# DIRT — Discrete-element Interaction-Resolved Toolkit

Discrete Element Method (DEM) physics in the **GRASS → SOIL → DIRT** stack.

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists
       └─ DIRT   DEM physics: contact, bonds, walls, clumps   ← you are here
```

DIRT resolves every inter-particle contact individually,
riding the [SOIL](https://github.com/SueHeir/soil) substrate for all the
method-agnostic machinery (atom data, domain decomposition, halo exchange,
neighbor lists). It adds the granular physics: Hertz–Mindlin contact, rotational
dynamics, parallel bonds, walls, multisphere clumps, heat conduction, and
contact analysis. 

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

## Building

Clone DIRT and build — its [soil](https://github.com/SueHeir/soil) and
[grass](https://github.com/SueHeir/grass) dependencies are pulled from GitHub
automatically, so you don't need to check them out yourself:

```bash
git clone https://github.com/SueHeir/dirt
cd dirt
cargo run --release --example hopper --no-default-features -- examples/hopper/config.toml
```

## License

MIT OR Apache-2.0
