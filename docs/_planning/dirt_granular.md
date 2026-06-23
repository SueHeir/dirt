# Planning: `dirt_granular` crate documentation

## Purpose

`dirt_granular` is the core contact-physics crate for DIRT. It implements the
full Hertz–Mindlin DEM contact loop (normal, tangential, rolling, twisting,
adhesion), rotational dynamics for angular DOF, and granular-temperature output.
Everything composes through the GRASS plugin system; the public entry point for
standard use is `GranularDefaultPlugins`.

Source files:
- `crates/dirt_granular/src/lib.rs` — crate-level doc, `GranularDefaultPlugins`, public re-exports
- `crates/dirt_granular/src/contact.rs` — fused Hertz-Mindlin contact (`hertz_mindlin_contact_force`) and Hooke path (`hooke_contact_force`)
- `crates/dirt_granular/src/tangential.rs` — `ContactHistoryStore` (the 7-slot per-contact spring)
- `crates/dirt_granular/src/rotational.rs` — `RotationalDynamicsPlugin`, `initial_rotation`, `final_rotation`
- `crates/dirt_granular/src/granular_temp.rs` — `GranularTempPlugin`, `print_granular_temperature`

---

## Public surface to document

### Plugin group
- **`GranularDefaultPlugins`** (`lib.rs:179`) — the standard DEM bundle.
  Registration order (lib.rs:182–190):
  1. `DemAtomPlugin` (from `dirt_atom`)
  2. `DemAtomInsertPlugin` (from `dirt_atom`)
  3. `VelocityVerletPlugin` (from `soil_verlet`)
  4. `HertzMindlinContactPlugin` (contact.rs:50)
  5. `RotationalDynamicsPlugin` (rotational.rs:49)

### Individual plugins
- **`HertzMindlinContactPlugin`** (`contact.rs:50`) — registers `ContactHistoryStore`
  in the `AtomDataRegistry` and routes to either `hertz_mindlin_contact_force` or
  `hooke_contact_force` at `ParticleSimScheduleSet::Force`. Dispatch is determined
  at `build()` time by reading `MaterialTable::contact_model` (contact.rs:70–89).
- **`RotationalDynamicsPlugin`** (`rotational.rs:49`) — registers two systems:
  `initial_rotation` at `ParticleSimScheduleSet::InitialIntegration` and
  `final_rotation` at `ParticleSimScheduleSet::FinalIntegration` (rotational.rs:52–56).
- **`GranularTempPlugin`** (`granular_temp.rs:33`) — opt-in; registers
  `print_granular_temperature` at `ParticleSimScheduleSet::PreExchange`.
  NOT included in `GranularDefaultPlugins` (lib.rs:112–114).

### Force systems
- `hertz_mindlin_contact_force` — Hertz normal + Mindlin tangential + constant/SDS rolling and twisting + JKR/DMT/SJKR adhesion/cohesion (contact.rs:103)
- `hooke_contact_force` — Hooke normal + Mindlin tangential + same rolling/twisting/cohesion; JKR/DMT surface_energy silently ignored (contact.rs:622)
- `initial_rotation` / `final_rotation` — half-step angular velocity Verlet (rotational.rs:58, 101)
- `print_granular_temperature` — computes `T_g = Σ m_i |v_i - v̄|² / (3M)`, writes `data/GranularTemp.txt` (granular_temp.rs:44)

### Key types
- `ContactHistoryStore` (`tangential.rs:49`) — per-contact `Vec<(u32, [f64;7], bool)>` keyed by `partner_tag`, canonical-frame; implements `AtomData` (pack/unpack for MPI migration)
- `MaterialTable` (in `dirt_atom`) — per-material and per-pair `*_ij` tables; read by all force kernels
- `GranularDefaultPlugins` (lib.rs:179)

### Constants
- `TANGENTIAL_EPSILON = 1e-30` — guards division when normalizing spring displacement (lib.rs:137)
- `LARGE_OVERLAP_WARN_THRESHOLD = 0.90` — ratio `dist / (r1+r2)` below which warns (lib.rs:145)
- `MAX_OVERLAP_WARNINGS = 500` — warnings before panic (lib.rs:152)
- `SQRT_5_6` (re-exported from `dirt_atom`) — damping factor in Tsuji formula

---

## Config / TOML schema

### `[dem]` section (global model selectors)

| Key | Type | Default | Meaning |
|---|---|---|---|
| `contact_model` | `"hertz"` \| `"hooke"` | `"hertz"` | Normal force law. Hertz: nonlinear `F_n = (4/3)E*√(R*δ)·δ`; Hooke: linear `F_n = kn·δ`. |
| `adhesion_model` | `"jkr"` \| `"dmt"` | `"jkr"` | Adhesion law when `surface_energy > 0`. Only consulted on Hertz path. |
| `rolling_model` | `"constant"` \| `"sds"` | `"constant"` | Rolling resistance type. |
| `twisting_model` | `"constant"` \| `"sds"` | `"constant"` | Twisting resistance type. |
| `track_orientation` | `bool` | `false` | Integrate per-sphere quaternion. For pure spheres this is causally inert; skip to save cost. (rotational.rs:63–65) |

Deserialized into `DemConfig` in `dirt_atom` (`dirt_atom/src/lib.rs:206`), `deny_unknown_fields`.

### `[[dem.materials]]` section (per-material)

> Note: The TOML section is `[[dem.materials]]` (nested under `[dem]`), not
> `[[materials]]`. The lib.rs crate-level doc example (lib.rs:46–55) uses
> `[[materials]]` — this is WRONG and must be corrected. Actual configs and
> `config.md` use `[[dem.materials]]`.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | `String` | (required) | Material name; referenced by `[[particles.insert]]` |
| `youngs_mod` | `f64` | (required) | Young's modulus E [Pa] |
| `poisson_ratio` | `f64` | (required) | Poisson's ratio ν (0–0.5) |
| `restitution` | `f64` | (required) | Target COR (0–1); converted to damping β via exact Hertz inversion |
| `friction` | `f64` | `0.4` | Sliding (Coulomb) friction coefficient μ |
| `rolling_friction` | `f64` | `0.0` | Rolling friction coefficient μ_r (0 = disabled) |
| `twisting_friction` | `f64` | `0.0` | Twisting friction coefficient μ_tw |
| `cohesion_energy` | `f64` | `0.0` | SJKR cohesion energy density [J/m²] (0 = disabled) |
| `surface_energy` | `f64` | `0.0` | JKR/DMT surface energy γ [J/m²] (0 = disabled; Hertz path only) |
| `kn` | `f64` | `0.0` | Linear normal stiffness [N/m] for `contact_model = "hooke"` |
| `kt` | `f64` | `0.0` | Linear tangential stiffness [N/m] for Hooke path |
| `rolling_stiffness` | `f64` | `0.0` | SDS rolling spring stiffness [N·m/rad] |
| `rolling_damping` | `f64` | `0.0` | SDS rolling dashpot coefficient |
| `twisting_stiffness` | `f64` | `0.0` | SDS twisting spring stiffness [N·m/rad] |
| `twisting_damping` | `f64` | `0.0` | SDS twisting dashpot coefficient |

Deserialized into `MaterialConfig` in `dirt_atom/src/lib.rs:144`, `deny_unknown_fields`.
`cohesion_energy` and `surface_energy` cannot both be non-zero — the code calls
`process::exit` (reference/config.md lines 107–109).

---

## Key behaviors, invariants & gotchas

### 1. Restitution → damping inversion (VERIFIED)

`restitution` is the **target COR**, not a damping ratio. In
`MaterialTable::build_pair_tables()` (dirt_atom/src/lib.rs:609–621) the geometric
mean `e_ij` is fed to `hertz_beta_for_cor(e_ij)`, which bisects the monotone
`COR(β)` curve of an RK4-integrated head-on Hertz collision (dirt_atom/src/lib.rs:78–131).
This makes the **input restitution equal the realized COR**.

For Hooke the inversion uses the analytic formula `β = −ln(e)/√(π² + ln²e)`
(a distinct code path at build_pair_tables; see dirt_atom/src/lib.rs:620–621 and
comment lines 617–619). The older Tsuji polynomial fit (no longer in use) was
documented to overshoot, e.g. 0.95 → 0.965; the exact inversion removes that
bias. The memory file `dirt-restitution-nominal-vs-realized.md` notes this affects
all calibration.

### 2. The 7-slot canonical contact history (VERIFIED)

`ContactHistoryStore` stores per-contact `[f64; 7]` (tangential.rs:49–55):

| Indices | Model | Content |
|---|---|---|
| `[0..3]` | Mindlin tangential | Tangential spring vector |
| `[3..6]` | SDS rolling | Rolling angular displacement vector |
| `[6]` | SDS twisting | Twisting scalar displacement |

Rolling/twisting slots are zero under the default constant-torque models (lib.rs:96–98).
Storage is **canonical**: lower-tag particle's perspective; a `sign = ±1` factor
converts to/from the local `(i, j)` frame each step (contact.rs:355–368).
The active-flag (`bool` in the tuple) resets to `false` before each pair loop;
stale entries are pruned by `retain` after the loop (contact.rs:125–129, 586–589).

### 3. Two-stage Coulomb cap (VERIFIED)

Applied in contact.rs:378–400:
1. The **stored spring** is capped: `|k_t s| ≤ μ|F_n|` — this truncates the
   history that carries to the next step.
2. The **assembled force** `F_t = k_t s − γ_t v_t` is capped again at `μ|F_n|`.

The same two-stage pattern is applied to SDS rolling (contact.rs:463–476) and SDS
twisting (contact.rs:522–526) via spring rescaling after the slider cap.

### 4. Hooke ignores JKR/DMT (VERIFIED)

In `hooke_contact_force` the `surface_energy` field is **read** (`mat_i`/`mat_j`
table lookups) but the Hooke normal force branch only handles `cohesion_energy`
(contact.rs:744–749). There is no JKR or DMT code path in the Hooke function.
The crate-level doc (lib.rs:38–39) and physics/contact.md (lines 83–88) both warn
about this asymmetry.

### 5. GranularTempPlugin is opt-in (VERIFIED)

`GranularTempPlugin` is explicitly excluded from `GranularDefaultPlugins`
(lib.rs:112–114, 164–166). It writes `data/GranularTemp.txt` with columns:
`step  time  T_granular  KE_total  |p_total|` (granular_temp.rs:14–16). The file
is truncated at step 0 and appended thereafter (granular_temp.rs:113–126).
Output fires at `ParticleSimScheduleSet::PreExchange` every `thermo` steps.

### 6. Schedule set phases

| System | Phase |
|---|---|
| `initial_rotation` | `ParticleSimScheduleSet::InitialIntegration` |
| `hertz_mindlin_contact_force` / `hooke_contact_force` | `ParticleSimScheduleSet::Force` |
| `final_rotation` | `ParticleSimScheduleSet::FinalIntegration` |
| `print_granular_temperature` | `ParticleSimScheduleSet::PreExchange` |

### 7. JKR extended interaction range

JKR extends the neighbor interaction range by a pull-off distance
`δ_pulloff = (π²γ²R* / (4E*²))^{1/3}` (contact.rs:166–173). Particles interact
even when separated by a gap up to `δ_pulloff`. In the gap regime (JKR adhesion
only, `jkr_adhesion_only = true`) the force is `−F_adhesion = −(3/2)π γ R*` and
no tangential, rolling, or history update occurs (contact.rs:244–347).
DMT has no extended range (contact.rs:246).

### 8. Clump sub-sphere inv_mass handling

`hooke_contact_force` and `hertz_mindlin_contact_force` both guard for clump
sub-spheres (`inv_mass = 0`) by falling back to `1.0 / atoms.mass[i]`
(contact.rs:234–236). Same-body pairs (sub-spheres of the same rigid body) are
skipped entirely via `dirt_atom::same_body` (contact.rs:140).

### 9. Large-overlap guard

When `distance / (r1 + r2) < LARGE_OVERLAP_WARN_THRESHOLD (0.90)`, the code
warns but still computes the repulsive force (overlap capped at `0.5 * r_min` to
keep Hertz valid). After `MAX_OVERLAP_WARNINGS (500)` such pairs per step, the
simulation panics with an actionable message (contact.rs:197–215).

### 10. quaternion update is opt-in within rotational plugin

The `initial_rotation` system checks `MaterialTable::track_orientation`
(rotational.rs:74, 84). For pure-sphere runs this is `false` (default), which
skips the quaternion half — the half-step ω kick still happens regardless
(rotational.rs:80–82). Only enables for clumps or cases requiring orientation output.

---

## Tutorial outline

1. **Minimal granular run** — `CorePlugins + GranularDefaultPlugins`, one
   `[[dem.materials]]`, one `[[particles.insert]]`. Annotate how `restitution`
   becomes realized COR.
2. **Adding rolling and twisting** — enable non-zero `rolling_friction`; show
   difference between `constant` and `sds` rolling models.
3. **JKR adhesion** — set `surface_energy`, explain pull-off distance, extended
   range, why you need `contact_model = "hertz"`.
4. **Granular temperature output** — add `GranularTempPlugin` explicitly, explain
   why it is opt-in, show output columns.
5. **Contact model swap: Hooke** — when to use linear stiffness; warn about
   JKR/DMT silently ignored; show required `kn`/`kt` values.
6. **Clumps with rotational dynamics** — explain `track_orientation = true`,
   sub-sphere `inv_mass = 0` skip, same-body pair exclusion.

---

## Doc gaps — accuracy check against physics/contact.md

### Verified accurate in physics/contact.md
- Normal force table (lines 11–15): formulas match code (contact.rs:254–255, 298–319).
- Tangential force formula (line 27): `γ_t = 2β√(5/6)√(k_t m_r)` matches contact.rs:388.
- 7-slot history layout (lines 42–47): correct (verified against tangential.rs:43–48).
- Two-stage Coulomb cap (lines 54–58): correct (verified against contact.rs:378–400).
- Rolling/twisting table (lines 62–66): formulas match contact.rs:489–503.
- JKR/DMT adhesion table (lines 75–79): formulas match contact.rs:296–319.
- Hooke/Hertz adhesion asymmetry warning (lines 83–87): correct.
- GranularDefaultPlugins registration order (lines 107–113): matches lib.rs:182–190.
- GranularTempPlugin opt-in statement (lines 121–123): correct.

### Inaccuracy found: TOML section names in lib.rs crate-level doc (lib.rs:43–65)

The crate-level doc example uses `[[materials]]` and `[materials]`:
```toml
[[materials]]
name = "glass"
...
[materials]
contact_model = "hertz"
```
The actual TOML sections are `[[dem.materials]]` and `[dem]`, as shown in
`docs/src/reference/config.md` (lines 76–105) and every real config file
(e.g. `examples/bench_hertz_rebound/config.toml:29,32`). The README.md has the
same error (`[[materials]]` / `[materials]`).

**Action:** Correct the TOML examples in lib.rs (lines 43–65) and README.md to
use `[[dem.materials]]` / `[dem]`.

### Stub in reference/validation.md
`docs/src/reference/validation.md` line 9 is a stub placeholder: "Stub. Summarize
the validation discipline and link each benchmark." The nominal-vs-realized
restitution caveat on lines 13–17 is now out of date — it describes the old
polynomial-fit overshoot, but as of the damping bug fix commit (c3ecd67) the exact
inversion is in place, so realized COR = nominal COR. The caveat should be updated
to describe the current behavior.

### `GranularTemp` not MPI-complete at step 0 in some edge cases
`print_granular_temperature` calls `all_reduce_sum_f64` unconditionally but only
rank 0 writes the file. This is correct, but the doc (granular_temp.rs:1) does not
note that step-0 truncate happens on rank 0 only. Minor, but worth noting for
multi-rank users who read output on other ranks.

---

## Suggested placement in docs/src

| Document | Placement |
|---|---|
| `physics/contact.md` | Already exists; **accurate** — no structural changes needed; fix validation.md stub |
| `reference/materials.md` | Already exists; accurate |
| `reference/config.md` | Already exists; accurate — `[dem]` / `[[dem.materials]]` naming is correct here |
| API reference for `ContactHistoryStore` | Could be added as `reference/contact-history.md` or folded into materials.md |
| `GranularTempPlugin` how-to | Short section in `physics/diagnostics.md` under a "Granular temperature" heading |
| Tutorial pages | `getting-started/` — two new pages: "Minimal DEM run" and "Adhesion & cohesion" |

---

## Summary of required fixes (actionable)

1. **lib.rs:43–65** and **README.md**: Change `[[materials]]` → `[[dem.materials]]`
   and `[materials]` → `[dem]` in the TOML examples.
2. **reference/validation.md:9**: Expand the stub; update the restitution caveat
   to reflect the current exact-inversion behavior (no longer overshoots).
3. **physics/contact.md** or **reference/materials.md**: Note that the Hooke β
   inversion uses the analytic formula (not bisection), distinct from Hertz.
   Currently neither doc distinguishes the two paths.
