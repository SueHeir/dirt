# Anatomy of a Config File

Every DIRT simulation reads its parameters from a TOML file passed on the
command line. The `main.rs` chooses *which* physics plugins run; the config
supplies *every number* they need. Below is the hopper config, section by
section.

## Communication & domain

```toml
[comm]
processors_x = 1
processors_y = 1
processors_z = 1          # MPI rank grid; 1×1×1 = serial

[domain]
x_low = 0.0
x_high = 0.04             # box extents [m]
y_low = 0.0
y_high = 0.02
z_low = 0.0
z_high = 0.08
boundary_x = "fixed"
boundary_y = "periodic"   # periodic in y → a 2D-like slice
boundary_z = "fixed"
```

`[comm]` is the domain-decomposition grid: the product of the three numbers must
equal your MPI rank count. `[domain]` is the simulation box and its boundary
conditions — `fixed` walls off an axis, `periodic` wraps it.

## Neighbor lists

```toml
[neighbor]
skin_fraction = 1.1       # search-radius multiplier (1.0–1.5 typical)
bin_size = 0.005          # neighbor-bin width [m] ≥ largest particle diameter
```

The neighbor list is the substrate's job (SOIL), configured here. The `skin`
lets particles move a little between rebuilds; `bin_size` sizes the spatial
grid.

## Gravity

```toml
[gravity]
gz = -90.81               # m/s² — enhanced 10× here to settle the bed faster
```

A common DEM trick: gravity is exaggerated to shorten settling time. Validation
examples are careful about when this is and isn't legitimate.

## Materials

```toml
[[dem.materials]]
name = "glass"
youngs_mod = 8.7e9        # Young's modulus [Pa]
poisson_ratio = 0.3
restitution = 0.3         # coefficient of restitution [0–1]
friction = 0.5            # sliding friction coefficient
```

Materials are named; particles and walls reference them by name. The contact
stiffness comes from `youngs_mod` and `poisson_ratio`; `restitution` sets the
damping.

> **Calibration note:** the input `restitution` is the *target* coefficient of
> restitution — DIRT numerically inverts the exact Hertz `COR(β)` curve to find
> the damping that realizes it, so for a binary collision the measured COR
> matches the number you type. (An older Tsuji polynomial fit overshot — e.g.
> 0.95 realized ≈ 0.965 — but that bias is gone.) See the
> [Materials](../reference/materials.md) and [validation](../reference/validation.md)
> chapters for the inversion and the measured values.

## Particle insertion

```toml
[[particles.insert]]
material = "glass"        # must match a [[dem.materials]] name
count = 200
radius = 0.001            # [m]
density = 2500.0          # [kg/m³]
velocity_z = -1.0         # initial downward velocity [m/s]
region = { type = "block", min = [0.005, 0.0, 0.055], max = [0.035, 0.02, 0.075] }
```

Each `[[particles.insert]]` block drops `count` particles randomly into a
`region`. You can have several blocks to build mixtures or layers.

## Walls

```toml
[[wall]]
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 1.0
normal_y = 0.0
normal_z = 0.0
material = "glass"
```

Each `[[wall]]` is a plane defined by a point and an outward normal. The hopper
config has side walls, a ceiling, the funnel faces, and a named `blocker` that
the `main.rs` removes at runtime. Walls, cylinders, cones, and spheres are all
available — see the [Walls](../physics/walls.md) chapter.

## The pattern

Config is **declarative and additive**: arrays of tables (`[[dem.materials]]`,
`[[particles.insert]]`, `[[wall]]`) let you stack as many of each as you need.
Each plugin you added in `main.rs` reads its own section. The full schema is in
the [Configuration Reference](../reference/config.md).
