# Fixes, Gravity & Damping

A **fix** is a per-atom operation applied to a named group every timestep — add
a force, hold a velocity, freeze a boundary, damp toward equilibrium. DIRT's
rotational fixes (and gravity) live in the `dirt_fixes` crate; add `FixesPlugin`
for the group fixes and `GravityPlugin` for the body force.

## The fix pattern

Each fix targets atoms in a named [group](../getting-started/config-anatomy.md)
and is configured via a TOML array of tables. A group is selected by `group =
"name"`; the fix then modifies the forces or velocities of every atom in that
group at a specific point in the timestep.

> **Typos are hard errors.** Every fix section (and `[gravity]`) uses
> `deny_unknown_fields`, so a misspelled key — `gammaa`, `feeze`, a stray field
> — fails TOML deserialization at startup rather than being silently ignored.
> `gamma` (viscous) and `max_displacement` (nve_limit) are *required* and have no
> default; omitting them is also a deserialization error. The force/velocity
> components (`fx`, `vz`, …) all default to `0.0`.

## The available fixes

| Fix | TOML key | Description |
|---|---|---|
| `AddForceDef` | `[[addforce]]` | Adds a constant force vector to atoms |
| `SetForceDef` | `[[setforce]]` | Overwrites the force vector on atoms |
| `MoveLinearDef` | `[[move_linear]]` | Moves atoms at a constant velocity |
| `FreezeDef` | `[[freeze]]` | Full immobilization — zeros velocity, force, and (for DEM atoms) angular velocity and torque |
| `ViscousDef` | `[[viscous]]` | Velocity-proportional damping (`F = −γv`) |
| `NveLimitDef` | `[[nve_limit]]` | Caps per-step displacement (`vmax = max_displacement/dt`), direction-preserving |
| `GravityConfig` | `[gravity]` | Gravitational body force (`F = mg`) |

### Gravity

```toml
[gravity]
gx = 0.0      # acceleration in x (default 0.0)
gy = 0.0      # acceleration in y (default 0.0)
gz = -9.81    # acceleration in z (default -9.81)
```

`GravityPlugin` applies `F_i = m_i g` to every local atom in the `Force` phase.
Ghost atoms are not affected. (A common DEM trick is to *exaggerate* gravity to
shorten settling time; the validation suite is careful about when that is
legitimate.)

### Forces: add vs. set

`[[addforce]]` accumulates `(fx, fy, fz)` on top of the pair/bond forces;
`[[setforce]]` **replaces** the force, discarding anything computed for those
atoms — useful for driving boundary atoms.

```toml
[[addforce]]
group = "fluid"
fz = 0.1

[[setforce]]
group = "wall"
fx = 0.0
fy = 0.0
fz = 0.0
```

You may declare multiple `[[addforce]]` (and `[[setforce]]`) blocks; each is
applied independently. But both kinds run in `PostForce`, and the scheduler does
not order independently-registered systems within a set. If the **same atom** is
in a group targeted by both an `addforce` and a `setforce`, the result is
undefined (last writer wins, with no guaranteed order). Do not overlap groups
between the two unless you control the ordering yourself.

### Prescribed motion

`[[move_linear]]` sets a group's velocity to a constant before the position
update (`PreInitialIntegration`), then zeros the force after force computation
(`PostForce`) so the final integration cannot alter it. The result: the group
translates at the prescribed rate regardless of applied forces.

```toml
[[move_linear]]
group = "piston"
vz = -0.001
```

`move_linear` prescribes **translation only** — it does not touch angular
velocity or torque, so a moved sphere can still spin up under contact. If you are
driving a bonded cluster or a body that must not rotate, pair `move_linear` with
`freeze` on a separate anchor or zero the rotation explicitly.

### Equilibration: viscous damping and displacement caps

`[[viscous]]` adds `F = −γv` to dissipate kinetic energy toward static
equilibrium. `[[nve_limit]]` rescales any atom whose speed would carry it more
than `max_displacement` in one step so that `|v| ≤ max_displacement/dt`. The
rescale is **direction-preserving** (heading unchanged, only speed clamped) and
the count of limited atoms is written to thermo as `n_limited` (this column
appears only if a `Thermo` resource is registered — i.e. with `PrintPlugin`
present). Its typical use is stabilizing the first few steps of a packing seeded
with overlaps, where huge contact forces would otherwise launch atoms across the
box.

```toml
[[viscous]]
group = "all"
gamma = 0.1

[[nve_limit]]
group = "all"
max_displacement = 0.0001
```

## Freeze vs. pin: the DIRT/SOIL split

`[[freeze]]` is **full immobilization** — translation *and* rotation. It zeros
velocity and force, and (if `DemAtom` is registered) angular velocity, torque,
and angular momentum. Because velocity is held at zero, the position never
drifts: the Verlet update adds `dt · 0`, so a frozen atom stays exactly where it
started without an explicit position restore.

Zeroing the rotational state is what makes a frozen boundary particle a true
immovable contact partner: it cannot spin up under contact torque, which would
otherwise corrupt the relative surface velocity at the contact. (The
bond-cantilever benchmark relies on this — a free-spinning anchor would let the
whole chain swing.)

> **Freeze (DIRT) vs. pin (SOIL).** `[[freeze]]` here is full immobilization,
> rotation included. For a **translation-only** positional constraint that
> *restores* position from a captured value and leaves rotation free, use SOIL's
> `[[pin]]` fix (`SoilFixesPlugin`). The two are deliberately split across
> tiers: rotational state is a DEM concept, so the rotation-aware freeze lives in
> DIRT, while the purely positional pin lives in the physics-agnostic substrate.

```toml
[[freeze]]
group = "anchor"
```

## Schedule ordering

| Fix | Schedule phase |
|---|---|
| `gravity` | `Force` |
| `addforce`, `setforce`, `freeze`, `viscous` | `PostForce` |
| `move_linear` | `PreInitialIntegration` (set velocity) + `PostForce` (zero force) |
| `nve_limit` | `PostFinalIntegration` (clamp after the final integration) |
