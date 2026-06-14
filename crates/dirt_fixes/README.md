# dirt_fixes

Group-based atom manipulation fixes and gravity for [DIRT](https://github.com/SueHeir/dirt).

Per-atom fixes that modify forces and velocities each timestep. Each fix targets a named atom group and is configured via TOML arrays in the simulation config: add/set forces, freeze or pin atoms, prescribe motion, damp velocity, cap displacement, and apply gravity.

## Key Types

| Type | TOML key | Description |
|------|----------|-------------|
| `AddForceDef` | `[[addforce]]` | Adds a constant force vector to group atoms |
| `SetForceDef` | `[[setforce]]` | Overwrites the force vector on group atoms |
| `MoveLinearDef` | `[[move_linear]]` | Moves atoms at a prescribed constant velocity |
| `FreezeDef` | `[[freeze]]` | Zeros velocity and force (immobilizes atoms) |
| `PinDef` | `[[pin]]` | Hard position constraint — captures pos at setup, restores it every step |
| `ViscousDef` | `[[viscous]]` | Velocity-proportional damping (F = −γv) |
| `NveLimitDef` | `[[nve_limit]]` | Caps max displacement per timestep by scaling velocity |
| `GravityConfig` | `[gravity]` | Gravitational body force (F = mg) |

Plugins: `FixesPlugin` registers the group-based fixes; `GravityPlugin` registers the gravity body force. `FixesRegistry` holds all parsed fix definitions; `PinState` stores captured pin positions keyed by global atom tag.

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

### pin — Hard position constraint
```toml
[[pin]]
group = "anchor"
```
Unlike `freeze`, `pin` also restores each atom to its captured initial position every step, before the Verlet drift and after force computation. Position is captured lazily on the first step the group mask is populated and is keyed by global tag, so pinned atoms survive MPI migration. If `DemAtom` is registered, angular velocity, torque, and angular momentum are zeroed too — important for bonded-particle (BPM) anchors.

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
- **pin**: PreInitialIntegration and PostForce (restore pos, zero vel/force; capture on first populated-mask step)
- **addforce, setforce, freeze, viscous**: PostForce
- **nve_limit**: PostFinalIntegration
- **gravity**: Force

All fixes validate their group names at setup and print a summary on rank 0. `FixesPlugin` only registers a fix's per-step system when at least one definition is present. Ghost atoms (index ≥ nlocal) are skipped.

## License

MIT OR Apache-2.0
