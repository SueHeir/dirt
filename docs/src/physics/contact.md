# Contact: Hertz–Mindlin

DIRT's default contact law resolves every inter-particle overlap with a
nonlinear-elastic normal force, a history-dependent tangential (friction)
force, rolling resistance, twisting friction, and optional adhesion. The whole
family lives in the `dirt_granular` crate and is bundled into the
`GranularDefaultPlugins` group.

## The two normal models

| Model | Normal force | Damping |
|---|---|---|
| **Hertz** (default) | `F_n = (4/3) E* √(R* δ) · δ` | `∝ β √(S_n m_r)` (Tsuji) |
| **Hooke** | `F_n = k_n δ` (linear spring) | linear, `γ_n v_n` |

`E*` is the effective modulus mixed from each material's `youngs_mod` and
`poisson_ratio`; `R*` is the effective radius; `δ` is the overlap. Select the
model globally with `contact_model = "hertz"` or `"hooke"` in the `[dem]`
section.

## Tangential force: Mindlin with history

The tangential force is an **incremental spring-history** model with a Coulomb
friction cap `μ |F_n|`:

```text
F_t = k_t s − γ_t v_t        γ_t = 2 β √(5/6) √(k_t m_r)
```

The spring displacement `s` is stored per-contact and re-projected into the
current tangent plane each step (a contact's tangent plane rotates as the
particles move). This is what gives a frictional pile its memory — a grain that
has begun to slip remembers how far.

### The canonical-frame, 7-slot history

The Mindlin tangential force and the SDS rolling/twisting variants are all
history-dependent springs: a displacement is integrated across timesteps,
rotated to stay in the current tangent plane, and capped at a Coulomb limit.
That history lives in a per-contact store holding **7 `f64` per contact**:

- slots `[0..3]` — the tangential spring vector,
- slots `[3..6]` — the rolling spring vector,
- slot `[6]` — the twisting scalar.

(The rolling/twisting slots stay zero under the default constant-torque models;
they only carry state under the SDS models.)

Each entry is kept in **canonical form** — from the lower-tag particle's
perspective — so the spring is frame-consistent no matter which particle is `i`
versus `j` in the neighbor list. A `sign` factor of `±1` flips the canonical
spring into the local `(i, j)` frame each step.

The Coulomb friction limit is applied in **two stages**:

1. the **stored** spring is first capped so `|k_t s| ≤ μ|F_n|`, truncating the
   history that survives to the next step;
2. the **assembled** force `F_t = k_t s − γ_t v_t` is capped again at `μ|F_n|`.

## Rolling and twisting resistance

| Resistance | Default (`constant`) | `sds` (spring–dashpot–slider) |
|---|---|---|
| **Rolling** | `τ_r = μ_r \|F_n\| R*` opposing relative rolling | incremental rolling displacement with spring stiffness, viscous damping, and a Coulomb-style slider cap |
| **Twisting** | `τ_tw = μ_tw \|F_n\| R*` opposing relative twisting | incremental twist angle with spring, damping, and cap |

Select per-resistance with `rolling_model` and `twisting_model` in `[dem]`. The
SDS rolling model is matched 1:1 to LAMMPS' `pair_style granular ... rolling
sds`; see the [angle-of-repose benchmark](../reference/validation.md).

## Adhesion and cohesion

| Model | Force | Available under |
|---|---|---|
| **JKR** (Johnson–Kendall–Roberts) | pull-off `F = (3/2) π γ R*`, extended interaction range beyond contact | Hertz only |
| **DMT** (Derjaguin–Muller–Toporov) | constant attractive `F = 2π γ R*` during contact | Hertz only |
| **SJKR** (simplified) | `F = k_coh π δ R*`, proportional to contact area | both |

JKR and DMT are driven by `surface_energy`; SJKR by `cohesion_energy`. Select
between JKR and DMT with `adhesion_model` in `[dem]` (only consulted when
`surface_energy > 0`).

> **Hooke/Hertz adhesion asymmetry.** JKR and DMT adhesion are implemented
> **only on the Hertz contact path**. Under `contact_model = "hooke"` the
> `surface_energy` term is *silently ignored* — the linear-spring path applies
> SJKR cohesion (`cohesion_energy`) only. If you need JKR/DMT pull-off, use the
> default Hertz model.

## Nominal vs. realized restitution

The `restitution` you type for a material is the **target coefficient of
restitution (COR)**, not a damping ratio. DIRT inverts the exact head-on Hertz
collision to find the damping `β` that *realizes* that COR — so the input
restitution becomes the measured COR of a binary collision. See
[Materials & the MaterialTable](../reference/materials.md#restitution--damping)
for how the inversion works, and the
[validation chapter](../reference/validation.md) for why an older polynomial
fit overshot the nominal value.

## The `GranularDefaultPlugins` group

`GranularDefaultPlugins` is the DEM-physics half of a standard run (pair it with
`CorePlugins` for infrastructure). In registration order it adds:

- **`DemAtomPlugin`** — per-atom material properties (radius, density) and the
  `MaterialTable`.
- **`DemAtomInsertPlugin`** — random particle insertion from
  `[[particles.insert]]`.
- **`VelocityVerletPlugin`** — translational integration.
- **`HertzMindlinContactPlugin`** — the normal + tangential contact described
  above (labeled `"hertz_mindlin_contact"` in the scheduler).
- **`RotationalDynamicsPlugin`** — quaternion velocity-Verlet for angular
  degrees of freedom (`I = 2/5 m r²` for solid spheres).

> **Neither plugin group alone integrates motion.** Velocity Verlet lives in
> `GranularDefaultPlugins`, not `CorePlugins`. `CorePlugins` alone reads config,
> builds neighbor lists, and prints output but never moves a particle. For a
> non-granular method you would pair `CorePlugins` with your own integrator and
> force plugins instead.

Granular-temperature output (`GranularTempPlugin`, which writes
`data/GranularTemp.txt`) is **opt-in** — it is *not* part of
`GranularDefaultPlugins`. Add it explicitly when you want that file.

## See also

- [Materials & the MaterialTable](../reference/materials.md) — every parameter
  above, where it is stored, and the mixing rules.
- [Walls](walls.md) — the same Hertz–Mindlin machinery against wall surfaces.
- [Configuration Reference](../reference/config.md) — the `[dem]` and
  `[[dem.materials]]` schema.
