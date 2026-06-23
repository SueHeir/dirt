# Planning: `dirt_wall` Documentation

Source read: `crates/dirt_wall/src/lib.rs` (2849 lines), `crates/dirt_wall/README.md`,
`crates/dirt_wall/Cargo.toml`. Cross-referenced: `docs/src/physics/walls.md`,
`docs/src/getting-started/first-simulation.md`,
`docs/src/getting-started/config-anatomy.md`, `docs/src/reference/config.md`.

---

## Purpose

`dirt_wall` provides rigid boundary surfaces for DIRT DEM simulations. Four
geometry types — plane, cylinder, sphere, region — enforce particle confinement
or contact with external surfaces. Wall–particle contacts reuse the identical
Hertz–Mindlin machinery as particle–particle contacts (same `MaterialTable`
mixing tables, same damping and friction models). Plane walls additionally
support three prescribed motion modes (constant velocity, sinusoidal oscillation,
servo force-feedback). Any wall can be named and deactivated at runtime without
restarting the simulation.

---

## Public Surface to Document

### Plugin

| Item | Location | Notes |
|---|---|---|
| `WallPlugin` | `lib.rs:578` | Registers `Walls` resource + 3 systems; depends on `DemAtomPlugin` |

### Resource

| Item | Location | Key fields |
|---|---|---|
| `Walls` | `lib.rs:504` | `planes`, `cylinders`, `spheres`, `regions` (parallel `active` flags); `tangential_springs`, `rolling_springs` HashMaps; `time` f64 |

### Runtime Control API

| Method | Location | Signature |
|---|---|---|
| `Walls::deactivate_by_name` | `lib.rs:540` | `(&mut self, name: &str)` — searches all four wall-type vecs, sets parallel active flag to false for any matching name |

There is no `activate_by_name` — deactivation is one-way in the current code
(lib.rs:540–561). This is a doc gap and potential user surprise.

### Geometry Types (runtime structs)

| Struct | Location | Key fields |
|---|---|---|
| `WallPlane` | `lib.rs:351` | `point_{x,y,z}`, `normal_{x,y,z}`, `material_index`, `name`, `bound_{x,y,z}_{low,high}`, `velocity`, `motion`, `origin`, `force_accumulator`, `temperature` |
| `WallCylinder` | `lib.rs:421` | `axis` (usize 0/1/2), `center` ([f64;2] in the cross-section plane), `radius`, `lo`, `hi`, `inside`, `material_index`, `name`, `force_accumulator`, `temperature` |
| `WallSphere` | `lib.rs:460` | `center` ([f64;3]), `radius`, `inside`, `material_index`, `name`, `force_accumulator`, `temperature` |
| `WallRegion` | `lib.rs:484` | `region` (soil_core `Region`), `inside`, `material_index`, `name`, `force_accumulator`, `temperature` |

### Motion Types

| Enum variant | Location | Parameters |
|---|---|---|
| `WallMotion::Static` | `lib.rs:313` | default — no motion |
| `WallMotion::ConstantVelocity` | `lib.rs:315` | uses `WallPlane::velocity` [f64;3] |
| `WallMotion::Oscillate` | `lib.rs:317` | `amplitude` (m), `frequency` (Hz); position = `amplitude * sin(2π * freq * t)` along normal; velocity computed analytically |
| `WallMotion::Servo` | `lib.rs:322` | `target_force` (N), `max_velocity` (m/s), `gain` (m/s per N); `velocity = clamp(gain * (target - measured), ±max_velocity)` |

### Config Structs (TOML-deserialized)

| Struct | Location |
|---|---|
| `WallDef` | `lib.rs:219` — union struct for all wall types |
| `OscillateDef` | `lib.rs:186` |
| `ServoDef` | `lib.rs:204` |

### Systems

| System | Schedule slot | Location | Purpose |
|---|---|---|---|
| `wall_move` | `PreInitialIntegration` | `lib.rs:826` | Updates plane positions from motion mode; advances `walls.time` by dt |
| `wall_zero_force_accumulators` | `PreForce` | `lib.rs:881` | Zeros per-wall `force_accumulator` on all four wall types |
| `wall_contact_force` | `Force` (label `"wall_contact"`) | `lib.rs:1062` | Hertz + damping + adhesion + friction for all active walls |

---

## Config / TOML Schema

All walls are `[[wall]]` array-of-table entries. `material` is always required.

### Common keys (all types)

| Key | Type | Default | Notes |
|---|---|---|---|
| `type` | string | `"plane"` | `"plane"`, `"cylinder"`, `"sphere"`, `"region"` |
| `material` | string | required | Must match a `[[dem.materials]]` name; fatal error if not found |
| `name` | string | None | Optional; enables `deactivate_by_name` at runtime |
| `temperature` | f64 | None | Stored, never read by this crate — hook for external heat-transfer |

### Plane-specific keys

| Key | Type | Default | Notes |
|---|---|---|---|
| `point_x`, `point_y`, `point_z` | f64 | 0.0 | A point on the plane; updated by motion at runtime |
| `normal_x`, `normal_y`, `normal_z` | f64 | 0.0 | Outward normal (normalized at parse time); fatal if zero vector |
| `bound_{x,y,z}_{low,high}` | f64 | ±∞ | AABB restricting where the plane is active; particles outside bounds are skipped |
| `velocity` | [f64;3] | None | Constant velocity mode (m/s) |
| `oscillate` | `OscillateDef` | None | `{ amplitude, frequency }` — sinusoidal along normal |
| `servo` | `ServoDef` | None | `{ target_force, max_velocity, gain }` — proportional force controller |

Only one of `velocity`, `oscillate`, `servo` should be set; parse precedence is
`oscillate` > `servo` > `velocity` (lib.rs:745–766).

### Cylinder-specific keys

| Key | Type | Default | Notes |
|---|---|---|---|
| `axis` | string | `"z"` | `"x"`, `"y"`, or `"z"` (case-insensitive); fatal otherwise |
| `center` | [f64; 2] | required | Center in the 2D plane perpendicular to axis (e.g. [cx, cy] for z-axis) |
| `radius` | f64 | required | Cylinder radius (m) |
| `lo` | f64 | −∞ | Lower axial bound; particles below ignored |
| `hi` | f64 | +∞ | Upper axial bound; particles above ignored |
| `inside` | bool | false | `true` = particles confined inside; contact normal points toward axis |

### Sphere-specific keys

| Key | Type | Default | Notes |
|---|---|---|---|
| `center` | [f64; 3] | required | Sphere center [x, y, z] |
| `radius` | f64 | required | Sphere radius (m) |
| `inside` | bool | false | `true` = particles inside; contact normal points toward center |

### Region-specific keys

| Key | Type | Default | Notes |
|---|---|---|---|
| `region` | Region table | required | Any `soil_core::Region` shape (block, sphere, cylinder, cone, union, intersect, etc.) |
| `inside` | bool | false | `true` = particles live inside the region surface |

Region contact detection delegates to `Region::closest_point_on_surface` for
signed distance and outward normal (lib.rs:1527).

---

## Key Behaviors, Invariants, and Gotchas

### 1. Effective contact mechanics (infinite-mass / infinite-radius wall)

`R* = particle_radius`, `m* = particle_mass` — the wall contributes neither
radius nor mass to the contact. This means Hertz stiffness and Mindlin stiffness
are functions of the particle alone. (lib.rs:1–13, lib.rs:1121–1122)

### 2. Adhesion asymmetry by geometry

Plane walls support JKR (`surface_energy > 0`) including the extended pull-off
regime (`delta < 0`) and DMT (lib.rs:1129–1195). Cylinder, sphere, and region
walls **silently ignore `surface_energy`** — only `cohesion_energy` (SJKR) is
applied (lib.rs:1347–1354, lib.rs:1455–1462, lib.rs:1568–1575). This is
documented in the module-level doc and README but is an easy source of user
confusion: setting `surface_energy` on a material used by a cylinder wall does
nothing.

### 3. Twisting friction on plane walls only

`twisting_friction_ij` is applied exclusively in the plane-wall loop (lib.rs:1203–1217).
Cylinder, sphere, and region walls receive sliding and rolling friction but not
twisting. No warning is emitted if a non-zero twisting friction material is
assigned to a curved wall.

### 4. Motion is plane-only

`wall_move` iterates only `walls.planes` (lib.rs:831). Cylinder, sphere, and
region walls are permanently static. Servo and oscillation are unavailable for
curved/region geometries.

### 5. Deactivation is irreversible (one-way)

`deactivate_by_name` sets `active[i] = false` with no inverse method
(lib.rs:540–561). A deactivated wall cannot be re-enabled without direct field
mutation (`walls.active[i] = true`). This is not documented anywhere in the
current mdBook.

### 6. Tangential- and rolling-spring history: per-step rebuild

`tangential_springs` and `rolling_springs` are `std::mem::take`n at the start of
`wall_contact_force`, rebuilt into new maps during the loop, then swapped back
(lib.rs:1076–1081, lib.rs:1627–1629). Contacts that ended are automatically
pruned (no explicit deletion needed). The key is `(wall_kind: u8, wall_index:
usize, particle_tag: u32)` where `wall_kind` is 0=plane, 1=cylinder, 2=sphere,
3=region.

### 7. Overlap cap at 0.5 × radius

`delta` is capped at `0.5 * radius` in all four wall-type loops (e.g.
lib.rs:1138, lib.rs:1335, lib.rs:1443, lib.rs:1547). This prevents runaway
forces from particles tunneling through the wall, but should not occur in a
well-timestep'd simulation.

### 8. Plane bounding box (`bound_x/y/z_low/high`)

Plane walls can have an AABB active region (lib.rs:274–290). The `in_bounds`
check (lib.rs:396–403) skips particles outside it. This is how the hopper's
funnel walls (finite faces) are implemented without needing a `region` wall. The
bounds default to ±∞ (infinite plane).

### 9. Servo one-step lag

The servo reads `force_accumulator` that was accumulated the previous timestep
(lib.rs:858–869). `wall_zero_force_accumulators` runs in `PreForce` (before
contact), `wall_contact_force` accumulates in `Force`, and `wall_move` runs next
step's `PreInitialIntegration`. This means there is a one-timestep lag between
measured force and velocity adjustment — normal for a discrete controller.

### 10. Config errors are fatal (`process::exit(1)`)

All `[[wall]]` parse and validation errors print `ERROR:` to stderr and call
`std::process::exit(1)` (lib.rs:617–619, lib.rs:655–659, lib.rs:670–672,
lib.rs:738–740). This is consistent with the rest of DIRT's fail-fast convention.

### 11. Named wall usage pattern (hopper example)

`first-simulation.md` shows the canonical pattern: give a blocker wall
`name = "blocker"`, then in a system call `walls.deactivate_by_name("blocker")`
when the simulation reaches a threshold. The `Walls` resource is mutable
(`ResMut<Walls>`) and the system runs on every qualifying step.

### 12. `region` wall vs. dedicated geometry wall: force equality

The test `region_wall_force_matches_dedicated_sphere` (lib.rs:2745) confirms that
a `Region::Sphere` region wall and a dedicated `WallSphere` produce equal forces
to `1e-6` relative tolerance. Use dedicated types for clarity; use `region` for
complex shapes (cone, union, intersection).

---

## Tutorial Outline

The existing `docs/src/physics/walls.md` already covers most of this. A tutorial
section should walk through:

1. **Minimal floor** — one `[[wall]]` entry; run and observe no particles escape.
2. **Named blocker** — add `name = "blocker"`; write a system that calls
   `deactivate_by_name`; show the hopper discharge pattern.
3. **Cylindrical container** — `type = "cylinder"`, `inside = true`; contrast
   with four bounding plane walls; mention axial `lo`/`hi` gap.
4. **Funnel with bounded planes** — use `bound_z_low` / `bound_z_high` to clip
   inclined plane walls to a finite frustum.
5. **Servo compaction** — add a moving top wall (`servo = { target_force = …
   }`); read `walls.planes[i].force_accumulator` in a diagnostic system.
6. **Region wall for complex geometry** — show a `Region::Cone` region wall
   (same as in config-anatomy); note adhesion limitation.

---

## Doc Gaps

| Gap | Severity | Notes |
|---|---|---|
| `deactivate_by_name` is one-way; no `activate_by_name` | Medium | Not documented; users wanting to re-enable a wall must reach into `Walls::active` directly |
| Adhesion asymmetry by geometry (cylinder/sphere/region silently ignore `surface_energy`) | High | Mentioned in module-level docs and README, but not in the mdBook; easy to waste hours wondering why JKR has no effect |
| Twisting friction plane-only limitation | Low | Undocumented in mdBook; material with `twisting_friction > 0` assigned to cylinder wall silently has no twisting effect |
| Motion is plane-only | Low | Mentioned in `physics/walls.md` but not in `reference/config.md` schema table |
| Bounding box on plane walls (`bound_x/y/z_low/high`) | Medium | Not mentioned in `physics/walls.md` or `reference/config.md`; key for finite-face funnel walls |
| Servo one-step force lag | Low | Worth a callout in the servo section — the controller reads last step's force |
| `temperature` field is stored, never read | Low | Mentioned in module docs and `physics/walls.md`; should note which future crate is expected to read it |
| No `activate_by_name` | Medium | Only deactivation is supported; re-activation requires direct `Walls::active[i] = true` |
| `inside = false` (outside) mode for cylinder and sphere | Low | Exists in code but not demonstrated or explained in any example config |

---

## Suggested Placement

The `docs/src/physics/walls.md` file already exists and is well-structured.
It should be **extended**, not replaced, with:

- A **"Plane bounding box"** subsection under Wall Types explaining `bound_*`
  keys with a worked funnel example.
- An **"Adhesion by geometry"** callout box (currently in the code docs but
  absent from mdBook prose).
- A **"Twisting friction: plane only"** note in the friction table.
- A **"Runtime control"** subsection explaining that deactivation is one-way
  and how to re-enable via direct field access.
- A **"Servo timing"** note under the Servo motion entry.

The `docs/src/reference/config.md` `## Walls` section (currently 6 lines) should
expand to the full key table documented above, parallel to the detailed
`## DEM physics` section for materials.
