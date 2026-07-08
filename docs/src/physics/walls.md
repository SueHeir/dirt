# Walls

Walls are boundary surfaces particles collide with: planes, cylinders, spheres,
and arbitrary region surfaces. They reuse the same Hertz–Mindlin contact
machinery as particle–particle contact, plus runtime control (rename, deactivate,
move). Walls live in the `dirt_wall` crate; add `WallPlugin` to enable them.

## Contact mechanics

Wall contacts reuse the same per-pair mixing tables as particle–particle
contacts: the particle's material and the wall's `material` index a row of the
[`MaterialTable`](../reference/materials.md) (`e_eff_ij`, `g_eff_ij`, `beta_ij`,
`friction_ij`, `rolling_friction_ij`, `twisting_friction_ij`, …). Because a wall
is treated as **infinitely massive and infinitely flat**, the **effective radius
is the particle radius** (`R* = particle_radius`) and the reduced mass is the
particle mass.

Beyond the normal Hertz force + damping, walls apply:

- **Tangential (Mindlin sliding) friction** — incremental spring-history model
  with a per-contact tangential spring, Coulomb-capped at `μ |F_n|`. Supported
  by **all** wall types (plane, cylinder, sphere, region).
- **Rolling resistance** — `constant` (default) or `sds`, mirroring the
  particle–particle rolling model with the wall as a zero-spin second body.
  Supported by **all** wall types.
- **Twisting friction** — constant torque opposing particle spin about the local
  contact normal. Supported by **all** wall types.

Frictionless walls (`friction = 0`) are byte-for-byte unchanged from a
pure-normal contact. The tangential and rolling spring histories are stored on
the `Walls` resource and keyed by `(wall_kind, wall_index, particle_tag)`;
they are rebuilt each step, so contacts that end are pruned automatically.

### Adhesion asymmetry by geometry

> **Adhesion support depends on wall geometry.** Only **plane walls** support
> JKR/DMT adhesion (`surface_energy`), including the JKR extended-range pull-off
> regime. **Cylinder, sphere, and region walls ignore `surface_energy`** and
> apply only SJKR cohesion (`cohesion_energy`). Setting `surface_energy` on a
> material used by a curved or region wall emits a setup warning explaining that
> JKR/DMT `surface_energy` is plane-wall-only. Use a plane wall if you need
> JKR/DMT against a wall.

- **Plane walls** support JKR and DMT (`surface_energy`) *and* SJKR cohesion
  (`cohesion_energy`), including the JKR extended-range pull-off regime.
- **Cylinder, sphere, and region walls** support **SJKR cohesion only**
  (`cohesion_energy`); their `surface_energy` is not consulted, so JKR/DMT
  pull-off is unavailable on curved/region walls.

### Wall temperature (a hook, not yet read)

Every wall config accepts an optional `temperature` (K). `dirt_wall` **stores**
it on the wall but never reads it — it is a hook for an external heat-transfer
system to consult a wall's temperature. It has no effect on the contact force.

## Wall types

| Type | Description | Config key |
|---|---|---|
| **Plane** | Infinite flat plane defined by a point and unit normal | `type = "plane"` |
| **Cylinder** | Infinite cylinder along X/Y/Z axis with finite axial bounds | `type = "cylinder"` |
| **Sphere** | Sphere defined by center and radius | `type = "sphere"` |
| **Region** | Any `Region` shape used as a wall surface | `type = "region"` |

For curved/region walls, `inside = true` means particles live inside the surface
(the normal points inward); `inside = false` confines them outside.

### Bounded plane faces

Plane walls can be clipped to a finite axis-aligned active box with six optional
keys:

| Key | Default | Meaning |
|---|---|---|
| `bound_x_low` | `-inf` | Lowest particle-center `x` coordinate where the plane can apply contact. |
| `bound_x_high` | `+inf` | Highest particle-center `x` coordinate where the plane can apply contact. |
| `bound_y_low` | `-inf` | Lowest particle-center `y` coordinate where the plane can apply contact. |
| `bound_y_high` | `+inf` | Highest particle-center `y` coordinate where the plane can apply contact. |
| `bound_z_low` | `-inf` | Lowest particle-center `z` coordinate where the plane can apply contact. |
| `bound_z_high` | `+inf` | Highest particle-center `z` coordinate where the plane can apply contact. |

The bounds are checked against the particle center before the plane distance is
computed. Particles outside any bound do not see that plane at all. This makes a
plane wall act like a finite rectangular face, or like one face of a funnel,
while keeping the contact normal and distance calculation exactly planar.

The bounds are an activation mask only. Motion remains plane-only: `wall_move`
updates the plane's point and velocity for `velocity`, `oscillate`, or `servo`,
but it does not move, rotate, or reshape the bound box itself. If a moving
finite face is needed, choose bounds that cover the whole swept region or update
the `Walls` resource from your own system.

Minimal four-face funnel:

```toml
# +X side face, tilted inward toward the outlet
[[wall]]
type = "plane"
point_x = 0.04
normal_x = -0.94
normal_z = 0.34
bound_y_low = 0.0
bound_y_high = 0.04
bound_z_low = 0.02
bound_z_high = 0.12
material = "glass"

# -X side face
[[wall]]
type = "plane"
point_x = 0.00
normal_x = 0.94
normal_z = 0.34
bound_y_low = 0.0
bound_y_high = 0.04
bound_z_low = 0.02
bound_z_high = 0.12
material = "glass"

# +Y side face
[[wall]]
type = "plane"
point_y = 0.04
normal_y = -0.94
normal_z = 0.34
bound_x_low = 0.0
bound_x_high = 0.04
bound_z_low = 0.02
bound_z_high = 0.12
material = "glass"

# -Y side face
[[wall]]
type = "plane"
point_y = 0.00
normal_y = 0.94
normal_z = 0.34
bound_x_low = 0.0
bound_x_high = 0.04
bound_z_low = 0.02
bound_z_high = 0.12
material = "glass"
```

## Wall motion

| Motion | Description | Config |
|---|---|---|
| **Static** | Wall does not move (default) | — |
| **Constant velocity** | Wall translates at a fixed velocity each step | `velocity = [vx, vy, vz]` |
| **Oscillating** | Sinusoidal displacement along the wall normal | `oscillate = { amplitude, frequency }` |
| **Servo** | Proportional controller driving velocity to a target contact force | `servo = { target_force, max_velocity, gain }` |

> Wall motion is currently supported only for **plane walls**. The servo lid in
> `examples/wall_servo_lid` drives a plane down onto a bed until the contact
> force reaches the target, then is removed at runtime.

Each step, `wall_move` updates positions from the motion mode (in
`PreInitialIntegration`), `wall_zero_force_accumulators` resets per-wall force
sums (`PreForce`), and `wall_contact_force` computes the contact and accumulates
the scalar force a servo controller reads next step (`Force`).

## Named walls and runtime control

Give a wall a `name` and you can flip it on and off at runtime via
`Walls::deactivate_by_name` and `Walls::activate_by_name`. This is how the
hopper example removes a blocker to start discharge, and how the servo-lid
example releases the bed:

```rust
walls.deactivate_by_name("blocker");
```

If a staged run needs the same wall again, reactivate it by the same name:

```rust
walls.activate_by_name("blocker");
```

Both methods search plane, cylinder, sphere, and region walls. Unknown names are
ignored, so optional stages can call them without first checking whether a wall
was configured.

## TOML configuration

Walls are `[[wall]]` array-of-tables entries; each requires a `material` field
matching a name in `[[dem.materials]]`.

```toml
# Plane wall (floor at z=0, normal pointing up)
[[wall]]
type = "plane"
point_z = 0.0
normal_z = 1.0
material = "glass"
name = "floor"                  # optional, for runtime enable/disable

# Finite plane face: active only for particle centers inside the box
[[wall]]
type = "plane"
point_x = 0.04
normal_x = -1.0
bound_y_low = 0.0
bound_y_high = 0.04
bound_z_low = 0.02
bound_z_high = 0.12
material = "glass"

# Cylinder wall (particles confined inside a z-aligned cylinder)
[[wall]]
type = "cylinder"
axis = "z"
center = [0.005, 0.005]         # center in the XY plane
radius = 0.004
lo = 0.0                        # axial lo bound (default: -inf)
hi = 0.01                       # axial hi bound (default: +inf)
inside = true                   # particles live inside the cylinder
material = "glass"

# Sphere wall (particles confined inside a sphere)
[[wall]]
type = "sphere"
center = [0.005, 0.005, 0.005]
radius = 0.004
inside = true
material = "glass"

# Region wall (any Region shape as a wall surface)
[[wall]]
type = "region"
inside = true
material = "glass"
region = { type = "cone", center = [0.005, 0.005], axis = "z",
           rad_lo = 0.004, rad_hi = 0.002, lo = 0.0, hi = 0.01 }

# Moving wall with constant velocity
[[wall]]
type = "plane"
normal_z = 1.0
material = "glass"
velocity = [0.0, 0.0, -0.01]    # [vx, vy, vz]

# Oscillating wall (sinusoidal along normal)
[[wall]]
type = "plane"
point_z = 0.1
normal_z = 1.0
material = "glass"
oscillate = { amplitude = 0.001, frequency = 50.0 }

# Servo-controlled wall (adjusts velocity to reach target force)
[[wall]]
type = "plane"
point_z = 0.1
normal_z = -1.0
material = "glass"
servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }
```

> **Config errors are fatal.** A malformed `[[wall]]` entry (bad TOML, an
> unknown cylinder axis, a wrong-length `center`, a missing region, a zero
> normal, …) prints an `ERROR:` line and exits the process at setup rather than
> returning a `Result` — the run stops immediately and identically on every MPI
> rank. This is the same fail-fast convention used for materials.

## See also

- [Contact: Hertz–Mindlin](contact.md) — the underlying force law.
- [Materials & the MaterialTable](../reference/materials.md) — how wall and
  particle materials mix into the per-pair tables.
