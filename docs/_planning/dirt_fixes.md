# Planning: `dirt_fixes` documentation

Crate: `dirt_fixes`  
Source: `/Users/suehr/Documents/GitHub/dirt/crates/dirt_fixes/src/lib.rs`  
README: `/Users/suehr/Documents/GitHub/dirt/crates/dirt_fixes/README.md`  
Existing doc target: `/Users/suehr/Documents/GitHub/dirt/docs/src/physics/fixes.md`

---

## Purpose

`dirt_fixes` is the DEM-tier crate for per-atom fixes and gravity. It provides a collection of group-targeted operations that modify atom forces, velocities, and (for DEM atoms) rotational state at specific points in each timestep's schedule. Two plugins are offered:

- **`FixesPlugin`** — registers all group-based fixes (addforce, setforce, move_linear, freeze, viscous, nve_limit)
- **`GravityPlugin`** — registers the gravitational body force

The crate depends on `dirt_atom::DemAtom` for rotational state (omega, torque, ang_mom), which is what makes it a "DEM-tier" rather than a SOIL-tier crate. The translational positional constraint (`[[pin]]`) deliberately lives in SOIL's `soil_fixes` crate (`SoilFixesPlugin`); see the freeze-vs-pin split below.

---

## Public Surface to Document

### Config structs (all `pub`, `#[derive(Deserialize)]`, `#[serde(deny_unknown_fields)]`)

| Struct | TOML key | Required fields | Optional fields (default) |
|---|---|---|---|
| `AddForceDef` | `[[addforce]]` | `group` | `fx` (0.0), `fy` (0.0), `fz` (0.0) |
| `SetForceDef` | `[[setforce]]` | `group` | `fx` (0.0), `fy` (0.0), `fz` (0.0) |
| `MoveLinearDef` | `[[move_linear]]` | `group` | `vx` (0.0), `vy` (0.0), `vz` (0.0) |
| `FreezeDef` | `[[freeze]]` | `group` | — |
| `ViscousDef` | `[[viscous]]` | `group`, `gamma` | — |
| `NveLimitDef` | `[[nve_limit]]` | `group`, `max_displacement` | — |
| `GravityConfig` | `[gravity]` | — | `gx` (0.0), `gy` (0.0), `gz` (−9.81) |

### Resource

- **`FixesRegistry`** — holds all parsed `Vec<*Def>` for every fix type; stored as an app resource; read by per-timestep systems (lib.rs:233–246)

### Plugins

- **`FixesPlugin`** — `impl Plugin for FixesPlugin` (lib.rs:257). Parses all fix arrays from config, stores `FixesRegistry`, adds `setup_fixes` at `ScheduleSetupSet::PostSetup`, then conditionally adds per-type update systems only when at least one definition of that type is present (lib.rs:303–344).
- **`GravityPlugin`** — `impl Plugin for GravityPlugin` (lib.rs:635). Loads `GravityConfig` from the `[gravity]` TOML section; adds `apply_gravity` at `ParticleSimScheduleSet::Force` unconditionally (lib.rs:647–649).

### Prelude / exports

All struct and plugin names are `pub` at crate root; no dedicated `prelude` module. Downstream crates import `dirt_fixes::{FixesPlugin, GravityPlugin, GravityConfig, AddForceDef, ...}` directly.

---

## Config / TOML Schema

All fix sections use `#[serde(deny_unknown_fields)]` — an unrecognized key is a hard error, not a silent ignore (lib.rs:65, 95, 128, 163, 184, 217, 604).

### `[gravity]` (GravityPlugin)

```toml
[gravity]
gx = 0.0     # x-component of gravitational acceleration [m/s²] (default 0.0)
gy = 0.0     # y-component [m/s²] (default 0.0)
gz = -9.81   # z-component [m/s²] (default -9.81)
```

Body force: `F_i = m_i * (gx, gy, gz)`. Sign convention: negative gz → downward in +z-up frame (lib.rs:598–615, 655–660). Default is standard Earth gravity in −z.

### `[[addforce]]` (FixesPlugin)

```toml
[[addforce]]
group = "fluid"   # required: name of a declared atom group
fx = 0.1          # force x-component [N] (default 0.0)
fy = 0.0          # force y-component [N] (default 0.0)
fz = 0.0          # force z-component [N] (default 0.0)
```

Accumulates onto existing force: `force[i] += (fx, fy, fz)`. Multiple `[[addforce]]` blocks are allowed and each is applied independently (lib.rs:66–78, 435–451).

### `[[setforce]]` (FixesPlugin)

```toml
[[setforce]]
group = "wall"
fx = 0.0
fy = 0.0
fz = 0.0
```

Replaces force: `force[i] = (fx, fy, fz)`, discarding all previously computed pair/bond forces. Useful for driving boundary atoms (lib.rs:94–108, 455–471).

### `[[move_linear]]` (FixesPlugin)

```toml
[[move_linear]]
group = "piston"
vx = 0.0
vy = 0.0
vz = -0.001   # velocity [m/s] (default 0.0)
```

Prescribes constant velocity. Group, velocity components; all velocity fields are optional and default to 0.0 (lib.rs:127–141).

### `[[freeze]]` (FixesPlugin)

```toml
[[freeze]]
group = "anchor"   # required; no other fields
```

Only `group`. No force or velocity parameters (lib.rs:163–168).

### `[[viscous]]` (FixesPlugin)

```toml
[[viscous]]
group = "all"
gamma = 0.1   # damping coefficient [force/velocity = kg/s] (required, no default)
```

`gamma` is required — no default. Units are force/velocity, i.e., kg/s in SI (lib.rs:184–191).

### `[[nve_limit]]` (FixesPlugin)

```toml
[[nve_limit]]
group = "all"
max_displacement = 0.0001   # max distance per step [m] (required, no default)
```

`max_displacement` is required — no default. The effective velocity cap is `vmax = max_displacement / dt` recomputed from `atoms.dt` each step (lib.rs:217–225, 566).

### Group selection (all fixes)

Every fix targets a named group via `group = "<name>"`. Groups must be declared elsewhere in the config (e.g., `[[group]]` blocks in the `soil_core` config). Group names are validated at setup time (via `GroupRegistry::validate_name`) and an unknown group name **panics** with a descriptive message (lib.rs:355–375). The mask is a per-atom boolean array; only atoms where `group.mask[i] == true` are affected.

---

## Key Behaviors, Invariants & Gotchas

### Schedule phases (lib.rs:24–35, 326–344)

| Fix | Phase(s) | Why |
|---|---|---|
| `gravity` | `Force` | Applied before PostForce so fixes can override it |
| `addforce`, `setforce`, `freeze`, `viscous` | `PostForce` | Run after all pair/bond forces are accumulated |
| `move_linear` | `PreInitialIntegration` + `PostForce` | Pre: set velocity before Verlet position update; Post: zero force so FinalIntegration cannot alter velocity |
| `nve_limit` | `PostFinalIntegration` | Clamps velocity after the timestep's final integration step |

The two-phase design of `move_linear` is load-bearing: setting velocity in PreInitialIntegration ensures the position advances by exactly `v * dt`; zeroing force in PostForce ensures the final Verlet half-step cannot perturb that velocity (lib.rs:113–141, 415–431, 513–529).

### `freeze` zeros rotational state (lib.rs:479–509)

`freeze` is DEM-aware: if `DemAtom` is registered, it zeros `omega[i]`, `torque[i]`, and `ang_mom[i]` in addition to translational velocity and force. This is critical for frozen boundary/wall particles: without it, a frozen sphere would spin up under contact torques from rolling contacts, corrupting `v_rel` at the surface and therefore the friction force (lib.rs:143–168, 474–509).

### `freeze` vs. `pin` split (lib.rs:33–35, README)

- **`freeze` (this crate, DEM tier)**: zeros velocity + force (+ rotational state for DEM atoms); position stability is maintained by holding velocity at zero so Verlet integration moves the atom by `dt * 0`. Position is not explicitly restored — there is no captured "home" position.
- **`pin` (`soil_fixes`, SOIL tier)**: translation-only constraint; captures each atom's position at setup and restores it every step. Leaves rotation free. Used for BPM anchors and reaction-force abutments.
- The split is intentional: rotational state is a DEM concept (DemAtom), so the rotation-aware fix belongs in the DEM tier (lib.rs:33–35; docs/src/physics/fixes.md:106–112).

### Ghost atoms are skipped for gravity (lib.rs:655)

`apply_gravity` iterates only `0..atoms.nlocal`. Ghost atoms (index ≥ nlocal) are not given a gravitational force — gravity is applied on the owning rank (lib.rs:655–660).

### Group fixes iterate only `nlocal` (lib.rs:420, 442, 461, 486, 540, 564)

All group fix systems also iterate only `0..nlocal`, consistent with the ghost-atom convention. Ghost atom forces are not modified by fixes.

### `FixesPlugin` is zero-overhead when unused (lib.rs:303–344)

If a fix type has no definitions, its per-timestep system is never registered. The `has_any` guard returns early (storing only the empty registry) if no fixes at all are configured. This means adding `FixesPlugin` to a sim with no fix sections costs nothing at runtime.

### `gamma` and `max_displacement` have no defaults (lib.rs:184–225)

Both `ViscousDef::gamma` and `NveLimitDef::max_displacement` are required keys — omitting them is a TOML deserialization error. This is different from the force component fields (`fx`, `fy`, `fz`, `vx`, `vy`, `vz`) which all default to 0.0.

### Sign convention for gravity (lib.rs:598–615)

The default is `gz = -9.81`. Users must set negative values for downward acceleration in a +z-up frame. The TOML snippet in `GravityPlugin::default_config` makes this explicit with a comment (lib.rs:637–643). There is no gravity magnitude + direction decomposition — the vector `(gx, gy, gz)` is used directly.

### `nve_limit` is direction-preserving (lib.rs:554–585)

All three velocity components are multiplied by the same `vmax / |v|` scale factor. Division-by-zero on a zero-velocity atom is safe: the `if vmag > vmax` guard ensures the rescale branch is never entered for zero-speed atoms (lib.rs:573–579).

### `nve_limit` writes `n_limited` to thermo (lib.rs:583–585)

The count of atoms clamped on the current step is set via `thermo.set("n_limited", ...)`. This field appears in thermo output only if `Thermo` is registered as an app resource; the system uses `Option<ResMut<Thermo>>` so it silently skips if thermo is absent.

### `addforce` vs. `setforce` ordering within PostForce

Both run in `PostForce`. If the same atom is in a group targeted by both `[[addforce]]` and `[[setforce]]`, the last writer wins — and because GRASS's scheduler does not guarantee ordering between independently registered systems in the same set, the interaction is undefined. Do not overlap groups between addforce and setforce unless you control ordering explicitly.

### Group name validation panics at setup (lib.rs:355–375)

Unknown group names cause a panic via `GroupRegistry::validate_name` in `setup_fixes`, which runs in `ScheduleSetupSet::PostSetup`. Errors surface at startup, not at runtime.

---

## Tutorial Outline

A good tutorial for `physics/fixes.md` would progress:

1. **Intro: what a fix is** — per-atom, per-timestep, group-targeted; point to group declaration syntax.
2. **Gravity** — the simplest case; show the default `[gravity]` block; explain sign convention; note ghost-atom skip; mention the "exaggerate gravity for settling" DEM trick (already in existing doc at fixes.md:38–40).
3. **Adding and overriding forces** — `[[addforce]]` for body forces (pressure gradients, buoyancy proxies); `[[setforce]]` to drive boundary atoms with known force; show the key difference (accumulate vs. replace).
4. **Prescribed motion: move_linear** — typical piston/indenter use; explain the two-phase implementation (pre: set vel; post: zero force); warn that angular velocity is NOT prescribed (atom can still spin).
5. **Immobilization: freeze** — wall/anchor particles; explain the DemAtom rotational zeroing; contrast with `[[pin]]` for position restoration.
6. **Equilibration tools** — `[[viscous]]` for settling; `[[nve_limit]]` for explosive-overlap startup; show them together; note `n_limited` in thermo output.
7. **Freeze vs. pin callout box** — canonical place to explain the DIRT/SOIL split.
8. **Schedule ordering table** — reproduce the schedule table for reference.

---

## Doc Gaps

The following are absent from, or only partially covered in, the current `docs/src/physics/fixes.md`:

1. **`nve_limit` thermo integration** — the existing doc (fixes.md:76–90) mentions the direction-preserving rescale and `n_limited`, but does not say HOW to see `n_limited` in output or that it requires `Thermo` to be registered.
2. **`addforce` + `setforce` ordering hazard** — no warning about overlapping groups or ordering within PostForce; this is a real user footgun.
3. **`move_linear` does NOT prescribe angular velocity** — not mentioned anywhere; a user moving a bonded cluster might expect rotation to also be controlled.
4. **`gamma` and `max_displacement` are required (no defaults)** — the existing TOML examples do show them, but there is no explicit "this field is required" callout to explain the deserialization error users would see if they omit it.
5. **Ghost atom behavior** — the gravity section (fixes.md:36–40) mentions ghost atoms, but no equivalent note exists for the group fixes.
6. **Zero-overhead lazy system registration** — the `FixesPlugin` optimization (no system registered if no definitions) is nowhere documented; worth one sentence for users worried about overhead.
7. **`FixesRegistry` as a public resource** — users writing custom systems can query `FixesRegistry` directly; this is undocumented.
8. **Multiple fix instances of the same type** — the config allows multiple `[[addforce]]` blocks targeting different (or even the same) groups; the doc shows only one example each.

---

## Suggested Placement

**Primary page: `docs/src/physics/fixes.md`** — this already exists and has the right structure. It should be expanded in-place.

Specific additions:
- Add required-vs-optional callout to the TOML schema tables for `gamma` and `max_displacement`.
- Add a "note" block after the `[[nve_limit]]` section explaining `n_limited` thermo output.
- Add a warning callout after the `[[addforce]]`/`[[setforce]]` section about overlapping groups and undefined ordering.
- The freeze-vs-pin section (fixes.md:92–126) is already good; no structural change needed.

**Secondary: `docs/src/reference/config.md`** — the existing "Fixes & gravity" stub (config.md:128–145) is appropriate as a brief cross-reference; no expansion needed there.

No new pages are required. The schedule-ordering table already present at fixes.md:119–126 should have `nve_limit` added (it is currently missing from the table).

---

## Source File Reference

| Topic | File:line |
|---|---|
| Module-level doc (fix list, plugins, schedule summary) | lib.rs:1–35 |
| `AddForceDef` | lib.rs:50–78 |
| `SetForceDef` | lib.rs:80–108 |
| `MoveLinearDef` | lib.rs:110–141 |
| `FreezeDef` | lib.rs:143–168 |
| `ViscousDef` | lib.rs:170–191 |
| `NveLimitDef` | lib.rs:193–225 |
| `FixesRegistry` | lib.rs:228–246 |
| `FixesPlugin::build` (conditional system registration) | lib.rs:287–344 |
| `setup_fixes` (group name validation + rank-0 summary) | lib.rs:349–411 |
| `apply_move_linear_pre` | lib.rs:415–431 |
| `apply_add_force` | lib.rs:433–451 |
| `apply_set_force` | lib.rs:453–471 |
| `apply_freeze` (incl. DemAtom rotational zero) | lib.rs:473–509 |
| `apply_move_linear_post` | lib.rs:511–529 |
| `apply_viscous` | lib.rs:531–550 |
| `apply_nve_limit` (incl. thermo write) | lib.rs:552–586 |
| `GravityConfig` (default gz = −9.81) | lib.rs:590–629 |
| `GravityPlugin::build` | lib.rs:635–650 |
| `apply_gravity` (local-only, F = mg) | lib.rs:652–661 |
| Existing doc page | docs/src/physics/fixes.md |
| Existing config reference stub | docs/src/reference/config.md:128–145 |
