# dirt_fixes

Group-based atom manipulation fixes and gravity for [DIRT](https://github.com/SueHeir/dirt).

Per-atom fixes that modify forces and velocities each timestep. Each fix targets a named atom group and is configured via TOML arrays in the simulation config: add/set forces, freeze atoms, prescribe motion, damp velocity, cap displacement, and apply gravity.

> The hard *position* constraint `[[pin]]` is **not** part of this crate — it lives in SOIL's `soil_fixes` (`SoilFixesPlugin`, with `PinDef`/`PinRegistry`/`PinState`). See the `soil_fixes` README for its schema. `freeze` here is the full-immobilization counterpart; see the freeze-vs-pin contrast below.

## Key Types

| Type | TOML key | Description |
|------|----------|-------------|
| `AddForceDef` | `[[addforce]]` | Adds a constant force vector to group atoms |
| `SetForceDef` | `[[setforce]]` | Overwrites the force vector on group atoms |
| `MoveLinearDef` | `[[move_linear]]` | Moves atoms at a prescribed constant velocity |
| `FreezeDef` | `[[freeze]]` | Zeros velocity and force (immobilizes atoms; also angular vel/torque/momentum for DEM atoms) |
| `ViscousDef` | `[[viscous]]` | Velocity-proportional damping (F = −γv) |
| `NveLimitDef` | `[[nve_limit]]` | Caps per-step displacement by scaling velocity (`vmax = max_displacement/dt`), direction-preserving; writes `n_limited` to thermo |
| `GravityConfig` | `[gravity]` | Gravitational body force (F = mg) |

Plugins: `FixesPlugin` registers the group-based fixes; `GravityPlugin` registers the gravity body force. `FixesRegistry` holds all parsed fix definitions.

### freeze vs. pin

`freeze` (this crate) zeros an atom's *velocity and force* every step — for DEM atoms it also zeros angular velocity, torque, and angular momentum — so the atom stops accelerating but stays wherever the last integration step left it. `pin` (in `soil_fixes`) additionally captures each atom's position at setup and *restores* it every step, holding the atom at a fixed point regardless of accumulated drift. Use `freeze` to kill all motion in place; use `pin` to anchor atoms to an exact location (BPM anchors, reaction-force abutments).

## TOML Configuration

### addforce — Add constant force to group atoms
```toml
[[addforce]]
group = "particles"
fx = 0.1
fy = 0.0
fz = 0.0
```

### setforce — Overwrite force with a constant value
```toml
[[setforce]]
group = "wall"
fx = 0.0
fy = 0.0
fz = 0.0
```

### move_linear — Prescribe constant velocity
```toml
[[move_linear]]
group = "piston"
vx = 0.0
vy = 0.0
vz = -0.001
```

### freeze — Immobilize atoms
```toml
[[freeze]]
group = "frozen"
```

> Looking for `[[pin]]`? It lives in `soil_fixes`, not here — see that crate's README.

### viscous — Velocity-proportional damping (F = −γv)
```toml
[[viscous]]
group = "all"
gamma = 0.1
```

### nve_limit — Cap max displacement per timestep
```toml
[[nve_limit]]
group = "all"
max_displacement = 0.0001
```

### gravity — Body force (F = m**g**)
```toml
[gravity]
gx = 0.0
gy = 0.0
gz = -9.81
```

## Schedule Phases

- **move_linear**: PreInitialIntegration (set velocity), PostForce (zero force)
- **addforce, setforce, freeze, viscous**: PostForce
- **nve_limit**: PostFinalIntegration
- **gravity**: Force

All fixes validate their group names at setup and print a summary on rank 0. `FixesPlugin` only registers a fix's per-step system when at least one definition is present. Ghost atoms (index ≥ nlocal) are skipped.

## License

MIT OR Apache-2.0
