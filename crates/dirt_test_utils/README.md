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
| `ParticleFixture` | Non-empty DEM fixture builder with synchronized atom/DEM rows, safe timestep, material pair tables, CSR neighbours, and optional `App` scheduling. |

## Usage

```rust,ignore
use dirt_test_utils::{ParticleFixture, ParticleSpec};

// Default is a valid overlapping pair; it cannot silently run a zero-atom loop.
let fixture = ParticleFixture::new().build();
assert_eq!(fixture.atom.nlocal, 2);

// A wall/fix test can instead use a single declared particle.
let fixture = ParticleFixture::single(ParticleSpec::new(0, [0.0; 3], 0.001)).build();
```

## License

MIT OR Apache-2.0
