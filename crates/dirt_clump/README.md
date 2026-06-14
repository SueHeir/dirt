# dirt_clump

Multisphere/clump rigid body composites for non-spherical DEM particles.

## What it does

A **clump** is a rigid body composed of multiple overlapping spheres. Each sub-sphere
participates in normal contact detection as an ordinary atom, but its forces are
aggregated onto a rigid body that integrates translational and rotational
(Euler-equation) dynamics. Sub-sphere positions, velocities, and angular velocities
are then derived from the body's state each step.

There is **no phantom parent atom**: body state lives in the `MultisphereBodyStore`
resource, and each sub-sphere atom references its body through `body_id` in `ClumpAtom`.
The inertia tensor is computed analytically (parallel axis theorem) for non-overlapping
spheres or by Monte Carlo for overlapping ones, then diagonalized into principal moments
and a principal-axes quaternion. The plugin also handles ghost-cutoff extension, body
exchange across MPI ranks, and periodic-boundary wrapping of body centers of mass.

## Key types

| Item | Role |
| --- | --- |
| `ClumpPlugin` | Registers clump data, resources, and all integration/aggregation systems. Depends on `dirt_atom::DemAtomPlugin`. |
| `ClumpDef` / `ClumpSphereConfig` | A clump type and its constituent spheres (body-frame offset + radius). |
| `ClumpInsertConfig` / `ClumpTopConfig` | `[clump.insert]` and top-level `[clump]` TOML config. |
| `ClumpAtom` | Per-atom data: `body_id` and body-frame `body_offset`. |
| `ClumpRegistry` | Runtime registry of loaded clump definitions. |
| `MultisphereBody` / `MultisphereBodyStore` | Rigid body state and the resource that owns all bodies (with an O(1) ID→index map). |

Public helpers include `insert_clump`, `compute_inertia_tensor_analytical`,
`compute_inertia_tensor_montecarlo`, `diagonalize_inertia`, `quat_rotate`, `cross`,
`same_body`, and `is_body_atom`.

## Configuration

Clumps live under the top-level `[clump]` TOML section, separate from `[dem]`
(which uses `deny_unknown_fields`). Define clump types under `[[clump.definitions]]`
and insert them with `[[clump.insert]]`:

```toml
[[clump.definitions]]
name = "dimer"
spheres = [
    { offset = [-0.0003, 0.0, 0.0], radius = 0.001 },
    { offset = [0.0003, 0.0, 0.0], radius = 0.001 },
]

[[clump.insert]]
definition = "dimer"
count = 100
density = 2500.0
material = "glass"          # must match a [[dem.materials]] entry
velocity = 0.5             # optional: each component uniform in [-v, +v]
region = { type = "block", min = [0.001, 0.001, 0.001], max = [0.019, 0.019, 0.019] }
```

## Usage

```rust
use dirt_clump::ClumpPlugin;
use grass_app::App;

let mut app = App::new();
app.add_plugin(ClumpPlugin);
// ... add dirt_atom, contact, materials, etc.
app.run();
```

The plugin loads definitions from `[clump]`, inserts clumps from `[[clump.insert]]`,
aggregates sub-sphere forces/torques onto bodies, integrates the bodies, and
reconstructs sub-sphere kinematics each step.

## License

MIT OR Apache-2.0
