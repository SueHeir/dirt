# Configuration Reference

Every DIRT simulation reads its parameters from a TOML file. The `main.rs`
chooses *which* plugins run; the config supplies *every number* they need. Each
plugin reads its own section, so the schema is **additive** — you only write the
sections for the physics you enabled. This page enumerates the sections; for a
narrative walk-through see
[Anatomy of a Config File](../getting-started/config-anatomy.md).

> **TOML conventions.** A single section is `[section]`; a repeated section is
> an array of tables, `[[section]]`. Several sections use
> `deny_unknown_fields`, so a typo'd key is a hard error, not a silent default.

## Infrastructure (`CorePlugins`)

### `[comm]` — domain decomposition grid

```toml
[comm]
processors_x = 1
processors_y = 1
processors_z = 1     # product must equal the MPI rank count; 1×1×1 = serial
```

### `[domain]` — simulation box and boundaries

```toml
[domain]
x_low = 0.0
x_high = 0.04        # box extents [m]
y_low = 0.0
y_high = 0.02
z_low = 0.0
z_high = 0.08
boundary_x = "fixed"
boundary_y = "periodic"   # "fixed" walls off an axis; "periodic" wraps it
boundary_z = "fixed"
```

### `[neighbor]` — neighbor lists

```toml
[neighbor]
skin_fraction = 1.1   # search-radius multiplier (1.0–1.5 typical)
bin_size = 0.005      # neighbor-bin width [m] ≥ largest particle diameter
every = 1             # rebuild interval in steps
```

### `[run]` — run stages

```toml
[[run]]
name = "fill"
dt = 1.0e-5           # timestep [s]
steps = 200000
thermo = 2000         # thermo print interval
```

Use multiple `[[run]]` blocks for staged runs, or a single `[run]`.

### `[output]` / `[vtp]` — dump files

```toml
[output]
dir = "examples/my_run"

[vtp]
interval = 2000
```

## DEM physics (`GranularDefaultPlugins`)

### `[dem]` — global model selection

```toml
[dem]
contact_model = "hertz"     # "hertz" (default) or "hooke"
adhesion_model = "jkr"      # "jkr" (default) or "dmt"; consulted when surface_energy > 0
rolling_model = "constant"  # "constant" (default) or "sds"
twisting_model = "constant" # "constant" (default) or "sds"
track_orientation = false   # integrate per-sphere quaternion (default false; pure overhead for spheres)
```

| Key | Type | Default | Meaning |
|---|---|---|---|
| `contact_model` | string | `"hertz"` | Normal force law. `"hertz"`: nonlinear `F_n = (4/3) E*√(R*δ)·δ`. `"hooke"`: linear `F_n = kn·δ`. |
| `adhesion_model` | string | `"jkr"` | Adhesion law applied when a material's `surface_energy > 0`. `"jkr"` or `"dmt"`. Consulted on the Hertz path only. |
| `rolling_model` | string | `"constant"` | Rolling-resistance model: `"constant"` torque or `"sds"` spring-dashpot-slider. |
| `twisting_model` | string | `"constant"` | Twisting-resistance model: `"constant"` or `"sds"`. |
| `track_orientation` | bool | `false` | Integrate each sphere's orientation quaternion. Inert for pure spheres (forces depend on `ω`, not absolute orientation) — leave `false` unless something downstream reads orientation. |

`[dem]` uses `deny_unknown_fields` — a typo'd key is a hard error.

### `[[dem.materials]]` — named materials

```toml
[[dem.materials]]
name = "glass"
youngs_mod = 8.7e9        # Young's modulus E [Pa]
poisson_ratio = 0.3       # ν (0–0.5)
restitution = 0.95        # target coefficient of restitution (0–1)
friction = 0.4            # sliding (Coulomb) friction μ (default 0.4)
rolling_friction = 0.0    # rolling friction (0 = disabled)
twisting_friction = 0.0   # twisting friction (0 = disabled)
cohesion_energy = 0.0     # SJKR cohesion energy density [J/m²] (0 = disabled)
surface_energy = 0.0      # JKR/DMT surface energy γ [J/m²] (0 = disabled)
# Hooke-model linear stiffnesses (0 = use Hertz/Mindlin):
kn = 0.0                  # linear normal stiffness [N/m]
kt = 0.0                  # linear tangential stiffness [N/m]
# SDS rolling/twisting (0 = use constant model):
rolling_stiffness = 0.0   # k_roll [N·m/rad]
rolling_damping = 0.0     # γ_roll
twisting_stiffness = 0.0  # [N·m/rad]
twisting_damping = 0.0
```

| Key | Type | Req. | Default | Meaning |
|---|---|---|---|---|
| `name` | string | yes | — | Material name; referenced by `[[particles.insert]]` and `[[wall]]`. |
| `youngs_mod` | f64 | yes | — | Young's modulus E [Pa]. |
| `poisson_ratio` | f64 | yes | — | Poisson's ratio ν (0–0.5). |
| `restitution` | f64 | yes | — | **Target** coefficient of restitution (0–1); inverted to a damping ratio (see [Materials](materials.md)). |
| `friction` | f64 | no | `0.4` | Sliding (Coulomb) friction coefficient μ. |
| `rolling_friction` | f64 | no | `0.0` | Rolling-resistance coefficient μ_r (0 = disabled). |
| `twisting_friction` | f64 | no | `0.0` | Twisting friction coefficient μ_tw (0 = disabled). |
| `cohesion_energy` | f64 | no | `0.0` | SJKR cohesion energy density [J/m²] (0 = disabled). |
| `surface_energy` | f64 | no | `0.0` | JKR/DMT surface energy γ [J/m²] (0 = disabled; Hertz path only). |
| `kn` | f64 | no | `0.0` | Linear normal stiffness [N/m] for `contact_model = "hooke"` (0 = use Hertz). |
| `kt` | f64 | no | `0.0` | Linear tangential stiffness [N/m] for the Hooke path (0 = use Mindlin). |
| `rolling_stiffness` | f64 | no | `0.0` | SDS rolling spring stiffness [N·m/rad]. |
| `rolling_damping` | f64 | no | `0.0` | SDS rolling dashpot coefficient. |
| `twisting_stiffness` | f64 | no | `0.0` | SDS twisting spring stiffness [N·m/rad]. |
| `twisting_damping` | f64 | no | `0.0` | SDS twisting dashpot coefficient. |

`[[dem.materials]]` uses `deny_unknown_fields`. `cohesion_energy` and
`surface_energy` cannot both be set on one material — it is a fatal config error
(the process exits). See [Materials & the MaterialTable](materials.md) for the
mixing rules and the restitution→damping inversion.

### `[[particles.insert]]` — particle insertion

```toml
[[particles.insert]]
material = "glass"        # must match a [[dem.materials]] name
count = 200
radius = 0.001            # fixed value or a distribution
density = 2500.0          # [kg/m³]
velocity_z = -1.0         # directional initial velocity
region = { type = "block", min = [0.005, 0.0, 0.055], max = [0.035, 0.02, 0.075] }
seed = 0                  # deterministic insertion RNG (default 0)
```

A block selects one of three modes by its fields. **Random** (default) and
**rate-based** (`rate` present) keys:

| Key | Type | Mode | Default | Meaning |
|---|---|---|---|---|
| `material` | string | all | required | A `[[dem.materials]]` name. |
| `count` | usize | random | required | Particles inserted at setup. |
| `radius` | RadiusSpec | all | required | Fixed `f64` or a distribution table (see below). |
| `density` | f64 | all | required | Particle density [kg/m³]. |
| `velocity` | f64 | random | `0.0` | Random speed magnitude [m/s]. |
| `velocity_x/y/z` | f64 | random | `0.0` | Directional velocity offset [m/s] (additive). |
| `region` | Region | all | domain inset by max radius | Insertion region (see [Insertion](insertion.md)). |
| `seed` | u64 | random | `0` | RNG seed; fixes the global packing. |
| `rate` | usize | rate | — (triggers mode) | Particles inserted per interval. |
| `rate_interval` | usize | rate | `1` | Insert every N steps. |
| `rate_start` / `rate_end` | usize | rate | `0` / never | First / last step for insertion. |
| `rate_limit` | usize | rate | unlimited | Cap on total particles inserted. |

**File-based** (`source = "file"`) keys:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `source` | string | — | Must be `"file"`. |
| `file` | string | — | Path to the input file. |
| `format` | string | — | `"csv"`, `"lammps_dump"`, or `"lammps_data"`. |
| `material` | string | — | Default material for atoms without a `type_map` entry. |
| `density` / `radius` | f64 | — | Defaults; overridden by file columns/data where present. |
| `columns` | table | `{x=0,y=1,z=2,radius=3}` | Column-index map for CSV. |
| `type_map` | table | — | `{ 1 = "glass", 2 = "steel" }` — file integer types → materials. |
| `atom_style` | string | auto | LAMMPS data style: `"atomic"`, `"sphere"`, `"bpm/sphere"`. |

See [Particle Insertion](insertion.md) for the `RadiusSpec` distributions,
region types, the born-in-owner model, and file-format details.

## Fixes & gravity (`FixesPlugin`, `GravityPlugin`)

```toml
[gravity]
gz = -9.81

[[freeze]]
group = "anchor"

[[viscous]]
group = "all"
gamma = 0.1
```

| Section | Required | Optional (default) | Effect |
|---|---|---|---|
| `[gravity]` | — | `gx`/`gy`/`gz` (`0`/`0`/`−9.81`) | Body force `F = m·(gx,gy,gz)`. |
| `[[addforce]]` | `group` | `fx`/`fy`/`fz` (`0`) | `force[i] += (fx,fy,fz)`. |
| `[[setforce]]` | `group` | `fx`/`fy`/`fz` (`0`) | `force[i] = (fx,fy,fz)` (replaces). |
| `[[move_linear]]` | `group` | `vx`/`vy`/`vz` (`0`) | Prescribe constant velocity (translation only). |
| `[[freeze]]` | `group` | — | Zero velocity, force, and (DEM) `omega`/`torque`/`ang_mom`. |
| `[[viscous]]` | `group`, `gamma` | — | `F = −γv` damping. `gamma` is required. |
| `[[nve_limit]]` | `group`, `max_displacement` | — | Cap speed at `max_displacement/dt`. Required key, no default. |

Every fix section uses `deny_unknown_fields`; `gamma` and `max_displacement`
have no defaults. The full set is documented in
[Fixes, Gravity & Damping](../physics/fixes.md). Atom **groups** referenced by
fixes are declared with `[[group]]` blocks.

## Walls (`WallPlugin`)

```toml
[[wall]]
type = "plane"
point_z = 0.0
normal_z = 1.0
material = "glass"
name = "floor"            # optional, for runtime deactivate_by_name
```

Common keys (all wall types) and per-geometry keys:

| Key | Type | Geometry | Default | Meaning |
|---|---|---|---|---|
| `type` | string | all | `"plane"` | `"plane"`, `"cylinder"`, `"sphere"`, `"region"`. |
| `material` | string | all | required | A `[[dem.materials]]` name. |
| `name` | string | all | none | Optional; enables runtime `deactivate_by_name`. |
| `temperature` | f64 | all | none | Stored, never read — hook for external heat transfer. |
| `point_x/y/z` | f64 | plane | `0.0` | A point on the plane. |
| `normal_x/y/z` | f64 | plane | `0.0` | Outward normal (normalized; fatal if zero). |
| `bound_{x,y,z}_{low,high}` | f64 | plane | ±∞ | AABB restricting where the plane is active (finite faces). |
| `velocity` | [f64;3] | plane | none | Constant-velocity motion [m/s]. |
| `oscillate` | table | plane | none | `{ amplitude, frequency }` — sinusoidal along normal. |
| `servo` | table | plane | none | `{ target_force, max_velocity, gain }` — force controller. |
| `axis` | string | cylinder | `"z"` | `"x"`/`"y"`/`"z"`. |
| `center` | [f64;2] / [f64;3] | cyl / sphere | required | Center (2D for cylinder, 3D for sphere). |
| `radius` | f64 | cyl / sphere | required | Cylinder / sphere radius. |
| `lo` / `hi` | f64 | cylinder | ±∞ | Axial bounds. |
| `inside` | bool | cyl/sphere/region | `false` | `true` = particles confined inside. |
| `region` | Region table | region | required | Any `soil_core::Region` shape. |

Motion (`velocity`/`oscillate`/`servo`) is **plane-only**; curved and region
walls are permanently static. Adhesion via `surface_energy` is also plane-only.
Plane, cylinder, sphere, and region walls are documented in full in
[Walls](../physics/walls.md).

## Bonds (`DemBondPlugin`)

```toml
[bonds]
auto_bond = true
youngs_modulus = 1.0e9
shear_modulus  = 4.0e8
beta_normal = 0.05
```

`[bonds]` geometry and stiffness keys:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `auto_bond` | bool | `false` | Bond every pair within `bond_tolerance·(R_i+R_j)` at setup. |
| `bond_tolerance` | f64 | `1.001` | Sum-of-radii multiplier for auto-bond eligibility. |
| `bond_radius_ratio` | f64 | `1.0` | Bond cylinder radius as a multiple of `min(R_i,R_j)`. |
| `ghost_cutoff_multiplier` | f64 | `2.5` | Scales max bond r₀ to extend the MPI ghost cutoff; `0.0` disables (single-rank only). |
| `youngs_modulus` | Option\<f64\> | none | E [Pa]; derives `K_n = E·A/L`, `K_bend = E·I/L`. |
| `shear_modulus` | Option\<f64\> | none | G [Pa]; derives `K_t = G·A/L`, `K_tor = G·J/L`. |
| `normal_stiffness` | f64 | `0.0` | Direct `K_n` [N/m] (used if `youngs_modulus` absent). |
| `shear_stiffness` | f64 | `0.0` | Direct `K_t` [N/m]. |
| `twist_stiffness` | f64 | `0.0` | Direct `K_tor` [N·m/rad]. |
| `bending_stiffness` | f64 | `0.0` | Direct `K_bend` [N·m/rad]. |
| `beta_normal/shear/twist/bending` | f64 | `0.0` | Critical-damping ratios per channel. |
| `normal/shear/twist/bending_damping` | Option\<f64\> | none | Raw `γ` overrides (bypass the β calc). |
| `seed` | u64 | `0` | Per-bond threshold RNG seed (MPI-stable). |
| `file` / `format` | Option\<String\> | none | Explicit bond file; `format = "lammps_data"`. |

`[bonds.breakage]` (absent → bonds never break) takes `kind` (one of nine: from
`"unbreakable"` through `"axial_*"`, `"combined_*"`, `"interaction_linear_*"`)
plus per-channel `ThresholdDistribution` sub-tables (`constant`, `weibull`,
`crack_band`). `[bonds.plasticity.bending]` (`guo_bending` / `guo_trilinear` /
`piecewise`) and `[bonds.plasticity.axial]` (`piecewise`) make individual
channels inelastic. See [Parallel Bonds](../physics/bonds.md) for the full
breakage menu, plasticity variants, and crack-band regularization.

## Clumps (`[clump]`)

```toml
[[clump.definitions]]
name = "dimer"
spheres = [
    { offset = [-0.0003, 0.0, 0.0], radius = 0.001 },
    { offset = [ 0.0003, 0.0, 0.0], radius = 0.001 },
]

[[clump.insert]]
definition = "dimer"
count = 100
density = 2500.0
material = "glass"
```

| Section / Key | Type | Default | Meaning |
|---|---|---|---|
| `[[clump.definitions]].name` | string | required | Clump-type name; referenced by inserts. |
| `[[clump.definitions]].spheres` | array | required | Each `{ offset = [x,y,z], radius }` (body-frame offset from COM, m). |
| `[[clump.insert]].definition` | string | required | Must match a definition `name`. |
| `[[clump.insert]].count` | u32 | required | Number of clumps to insert. |
| `[[clump.insert]].density` | f64 | required | [kg/m³]; drives mass and inertia. |
| `[[clump.insert]].material` | string | required | A `[[dem.materials]]` name. |
| `[[clump.insert]].velocity` | f64 | `0.0` | Optional; each component drawn uniform in [−v, +v]. |
| `[[clump.insert]].region` | Region | domain inset by eff. radius | Optional insertion region. |

The `[clump]` section does *not* use `deny_unknown_fields` (unlike `[dem]`).
Overlapping-sphere clumps have stochastic, non-reproducible inertia (Monte Carlo,
100 000 samples). See [Clumps (Multisphere)](../physics/clumps.md).

## Diagnostics

```toml
[[measure_plane]]
name = "outlet"
point = [0.0, 0.0, 0.03]
normal = [0.0, 0.0, -1.0]
report_interval = 1000

[contact_analysis]
coordination = true
fabric_tensor = true
```

`[[measure_plane]]` keys (`deny_unknown_fields`):

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | string | required | Suffix for thermo keys (`crossings_<name>`, …). |
| `point` | [f64;3] | required | Any point on the plane. |
| `normal` | [f64;3] | required | Direction counted as a crossing (auto-normalized). |
| `report_interval` | usize | `1000` | Rate-averaging window in steps. |

`[contact_analysis]` keys (`deny_unknown_fields`):

| Key | Type | Default | Meaning |
|---|---|---|---|
| `interval` | usize | `0` | Dump per-contact CSV every N steps (`0` = disabled). |
| `coordination` | bool | `false` | Per-atom coordination → `coord_avg/max/min` thermo. |
| `rattlers` | bool | `false` | Detect <4-contact particles (requires `coordination = true`). |
| `fabric_tensor` | bool | `false` | Fabric tensor → six thermo components. |
| `file_prefix` | string | `"contact"` | Prefix for CSV filenames. |

`ContactAnalysisPlugin` must be registered **after** `CorePlugins` and
`GranularDefaultPlugins` (it panics otherwise). See
[Diagnostics](../physics/diagnostics.md) for the thermo keys and caveats.

## Generating a schema

Each plugin can emit its own default-config snippet via the framework's
`default_config` hook; running a binary with the generate-config flag prints the
assembled schema for exactly the plugins that binary added. The
[example configs](https://github.com/SueHeir/dirt/tree/master/examples) are also
heavily commented and remain a useful authoritative reference.
