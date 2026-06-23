# dirt_wall

Wall contact forces for DIRT (DEM) simulations: Hertz normal contact with viscous damping, tangential/rolling/twisting friction, and optional adhesion (JKR, DMT, or SJKR cohesion).

## What it does

Parses `[[wall]]` entries from the TOML config, resolves each wall's material against the `[[dem.materials]]` table, and computes wall–particle contact forces every timestep. Walls have infinite mass and infinite radius, so the effective radius is the particle radius and the reduced mass is the particle mass. Plane walls can move (constant velocity, sinusoidal oscillation, or servo force-feedback); cylinder, sphere, and region walls are static. Walls carry an optional `name` for runtime enable/disable and an optional `temperature` field.

## Friction & adhesion

Wall contacts reuse the same per-pair `MaterialTable` tables as particle–particle contacts (`friction_ij`, `rolling_friction_ij`, `twisting_friction_ij`, `e_eff_ij`, `g_eff_ij`, `beta_ij`), with `R* = particle_radius`.

- **Tangential (Mindlin sliding) friction** and **rolling resistance** (`constant` or `sds`) are applied by **all** wall types. Frictionless walls (`friction = 0`) are byte-for-byte identical to a pure-normal contact.
- **Twisting friction** is applied by **plane walls only**.
- **Adhesion is asymmetric by geometry:** plane walls support JKR/DMT (`surface_energy`) *and* SJKR cohesion (`cohesion_energy`); cylinder/sphere/region walls support **SJKR cohesion only** — their `surface_energy` is ignored.
- The optional **`temperature`** field is *stored but never read* by this crate — it is a hook for an external heat-transfer system.

Malformed `[[wall]]` config (bad TOML, unknown cylinder axis, wrong-length `center`, missing region, zero normal) prints an `ERROR:` and calls `std::process::exit(1)` at setup.

## Wall types

| Type (`type =`) | Description |
|------|-------------|
| `"plane"` (default) | Arbitrary-orientation infinite plane with optional XYZ bounding box and motion |
| `"cylinder"` | Axis-aligned cylinder with finite axial bounds and inside/outside modes |
| `"sphere"` | Sphere with inside/outside modes |
| `"region"` | Any `soil_core` `Region` shape used as a wall surface |

Plane motion modes: static, constant velocity, sinusoidal oscillation along the normal, and proportional servo (drives velocity toward a target contact force).

## Key types

| Item | Role |
|------|------|
| `WallPlugin` | Registers the `Walls` resource and contact/motion systems; depends on `DemAtomPlugin` |
| `WallDef`, `OscillateDef`, `ServoDef` | TOML config structs |
| `WallPlane`, `WallCylinder`, `WallSphere`, `WallRegion` | Runtime wall representations |
| `WallMotion` | `Static`, `ConstantVelocity`, `Oscillate`, `Servo` |
| `Walls` | Resource holding all walls with per-wall active flags; `deactivate_by_name` |
| `wall_move`, `wall_zero_force_accumulators`, `wall_contact_force` | Systems (PreInitialIntegration / PreForce / Force) |

## TOML configuration

```toml
# Plane floor at z=0, normal pointing up
[[wall]]
type = "plane"
point_z = 0.0
normal_z = 1.0
material = "glass"     # must match a [[dem.materials]] name
name = "floor"         # optional, for runtime enable/disable

# Cylindrical container (particles confined inside)
[[wall]]
type = "cylinder"
axis = "z"
center = [0.005, 0.005]
radius = 0.004
lo = 0.0
hi = 0.01
inside = true
material = "glass"

# Servo-controlled wall (adjusts velocity to reach target force)
[[wall]]
type = "plane"
point_z = 0.1
normal_z = -1.0
material = "glass"
servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }
```

## Usage

```rust
use dirt_wall::WallPlugin;

app.add_plugins(WallPlugin);
```

`WallPlugin` requires `DemAtomPlugin` (for the `MaterialTable` and `DemAtom` data) to be added first. Walls are then defined entirely through `[[wall]]` TOML entries.

## License

MIT OR Apache-2.0
