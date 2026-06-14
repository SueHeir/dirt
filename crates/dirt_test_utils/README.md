# dirt_test_utils

Shared test fixtures for the DIRT crates.

## What it does

Provides helper functions that build minimal, valid `Atom`, `GroupRegistry`, `CommResource`, `DemAtom`, and `MaterialTable` instances so unit tests can focus on the logic under test instead of setup boilerplate.

## Key functions

| Function | Purpose |
|----------|---------|
| `make_atoms(n)` | `Atom` with `n` atoms along the x-axis at `(i, 0, 0)`, each radius `0.5`, mass `1.0`, tags `0..n`, `dt = 0.001`. |
| `make_group_registry(name, mask)` | `GroupRegistry` with one named group; `mask: Vec<bool>` sets per-atom membership and the count is derived from it. |
| `make_single_comm()` | Single-process `CommResource` backed by `SingleProcessComm` (no MPI). |
| `push_dem_test_atom(atom, dem, tag, pos, radius)` | Pushes a solid sphere (ρ = 2500 kg/m³) into parallel `Atom` and `DemAtom` structures, filling radius, density, inverse inertia, quaternion, omega, angular momentum, torque, and body id. |
| `make_material_table()` | `MaterialTable` with one `"glass"` material (E = 8.7 GPa, ν = 0.3, e = 0.95, μ = 0.4) and pre-built pair tables. |

## Usage

```rust,ignore
use dirt_test_utils::{make_atoms, make_group_registry, make_single_comm};
use dirt_test_utils::{push_dem_test_atom, make_material_table};

// Shared substrate fixtures
let atoms = make_atoms(3);
let groups = make_group_registry("all", vec![true, true, true]);
let comm = make_single_comm();

// DEM fixtures
let mut atom = soil_core::Atom::new();
let mut dem = dirt_atom::DemAtom::default();
push_dem_test_atom(&mut atom, &mut dem, 0, [0.0, 0.0, 0.0], 0.5);

let materials = make_material_table();
```

## License

MIT OR Apache-2.0
