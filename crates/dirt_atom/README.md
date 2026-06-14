# dirt_atom

Per-atom DEM data and material property tables for DIRT.

## What it does

`dirt_atom` provides the core per-atom data structures, material property system, and
particle insertion for Discrete Element Method (DEM) simulations on the SOIL substrate:

- **`DemAtom`** — Per-atom extension data (radius, density, inverse inertia, orientation,
  angular velocity/momentum, torque, rigid-body id), registered via SOIL's `AtomData` derive
  with `#[forward]` / `#[reverse]` / `#[zero]` attributes controlling MPI pack/unpack and
  per-timestep zeroing.
- **`MaterialTable`** — Named materials with mechanical properties (Young's modulus,
  Poisson's ratio, restitution, friction, cohesion/surface energy, rolling/twisting) plus
  precomputed per-pair mixing tables for contact-force evaluation.
- **`DemAtomPlugin`** — Registers `DemAtom` and builds the `MaterialTable` from
  `[[dem.materials]]` config.
- **`DemAtomInsertPlugin`** — Particle insertion from `[[particles.insert]]`: random
  placement (with overlap checking), rate-based trickle insertion, or file-based loading
  (CSV, LAMMPS dump, LAMMPS data).
- **`RadiusSpec`** — Fixed radius or statistical distribution (uniform, gaussian, lognormal,
  discrete).

## Key types and plugins

| Item | Kind | Role |
|------|------|------|
| `DemAtomPlugin` | Plugin | Registers `DemAtom`, builds `MaterialTable`, sets `Atom::ntypes` |
| `DemAtomInsertPlugin` | Plugin | Setup + runtime systems for random/rate/file insertion |
| `DemAtom` | AtomData | Per-atom radius, density, `inv_inertia`, `quaternion`, `omega`, `ang_mom`, `torque`, `body_id` |
| `MaterialTable` | Resource | Per-material + per-pair (`*_ij`) precomputed contact properties |
| `MaterialConfig` / `DemConfig` | Config | Deserialized `[[dem.materials]]` / `[dem]` |
| `InsertConfig` / `ParticlesConfig` | Config | Deserialized `[[particles.insert]]` / `[particles]` |
| `RadiusSpec` / `RadiusDistribution` | Config | Fixed value or sampled radius distribution |
| `same_body(dem, i, j)` | fn | True when atoms `i` and `j` share a rigid body |

## AtomData attributes on `DemAtom`

- **`#[forward]`** — sent from owner to ghost atoms each timestep (`omega`, `body_id`).
- **`#[reverse]`** — ghost contributions summed back into the owner (`torque`).
- **`#[zero]`** — zeroed before each force computation (`torque`).

Fields without attributes are communicated only during atom migration.

## TOML configuration

```toml
[dem]
contact_model = "hertz"        # "hertz" (default) or "hooke"
adhesion_model = "jkr"         # "jkr" (default) or "dmt", used when surface_energy > 0
rolling_model = "constant"     # "constant" (default) or "sds"
twisting_model = "constant"    # "constant" (default) or "sds"

[[dem.materials]]
name = "glass"
youngs_mod = 8.7e9
poisson_ratio = 0.3
restitution = 0.95
friction = 0.4                 # default 0.4
# rolling_friction = 0.1       # default 0.0 (disabled)
# cohesion_energy = 0.05       # SJKR cohesion, J/m² (default 0.0)
# surface_energy = 0.05        # JKR/DMT adhesion, J/m² (default 0.0; mutually exclusive with cohesion_energy)
# twisting_friction = 0.0
# kn = 0.0 / kt = 0.0          # Hooke stiffnesses (0 = use Hertz/Mindlin)

[[particles.insert]]
material = "glass"             # must match a [[dem.materials]] name
count = 100
radius = 0.001                 # or a distribution table
density = 2500.0
```

## Usage

```rust
use dirt_atom::MaterialTable;

let mut table = MaterialTable::new();
let glass = table.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.1, 0.0);
let steel = table.add_material("steel", 200e9, 0.28, 0.8, 0.3, 0.05, 0.0);
table.build_pair_tables();

// Per-pair properties indexed as table_ij[i][j]:
let friction = table.friction_ij[glass as usize][steel as usize];
let e_eff = table.e_eff_ij[glass as usize][steel as usize];
```

Per-pair tables are mixed with **geometric mean** (friction, restitution, rolling/twisting
friction, cohesion/surface energy, SDS damping), **harmonic mean** (Hooke and SDS stiffnesses),
and **effective-modulus** formulas for Hertz (`E*`) and Mindlin (`G*`).

## License

MIT OR Apache-2.0
