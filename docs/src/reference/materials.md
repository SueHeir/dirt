# Materials & the MaterialTable

Every contact force in DIRT reads its parameters from the **`MaterialTable`** —
the named materials you declare under `[[dem.materials]]` plus the precomputed
per-pair mixing tables the force kernels actually index. This page is the
conceptual home for the material system: the two-phase build, the mixing rules,
the exact restitution→damping inversion, and which parameter feeds which model.

## Materials are named; pairs are mixed

You declare materials by name; particles and walls reference them by name. The
contact between two materials `i` and `j` is governed not by the per-material
values directly but by a **mixed** per-pair value (`e_eff_ij`, `friction_ij`,
…). This is what lets glass-on-steel differ from glass-on-glass without
declaring a third material.

## The two-phase build contract

A `MaterialTable` is filled in **two phases**, and the contact force code reads
only the second-phase output:

1. **Register materials.** Each `add_material*` call appends one row to the
   per-material vectors (`youngs_mod`, `restitution`, `friction`, …) and returns
   its integer index. During this phase **every `*_ij` pair table is empty.**
2. **Build pair tables.** `build_pair_tables()` allocates the `N×N` `*_ij`
   tables and fills them from the registered per-material values using the
   mixing rules below. It **must** be called once, after the last material is
   added and before any contact force is evaluated. Indexing a `*_ij` table
   before this is an out-of-bounds panic.

When you drive DIRT from a TOML config you never see this — `DemAtomPlugin`
registers the materials and calls `build_pair_tables()` for you. The two-phase
pattern only surfaces if you build a table by hand in a test or tool:

```rust
use dirt_atom::MaterialTable;

let mut mat = MaterialTable::new();

// Phase 1 — register. Returns the material index.
let glass = mat.add_material(
    "glass",
    8.7e9, // Young's modulus E [Pa]
    0.3,   // Poisson's ratio ν
    0.95,  // restitution (target COR)
    0.5,   // sliding friction
    0.0,   // rolling friction
    0.0,   // cohesion energy
);
assert!(mat.beta_ij.is_empty()); // pair tables still empty in phase 1

// Phase 2 — build the per-pair mixing tables. Required before contact eval.
mat.build_pair_tables();
let beta = mat.beta_ij[glass as usize][glass as usize];
```

### The `add_material*` ladder

Four constructors form a wrapping ladder from fewest to most arguments; each
delegates to the next with sensible zero defaults. Use the shortest one that
covers the parameters you need:

| Constructor | Adds |
|---|---|
| `add_material` | basics (E, ν, restitution, friction, rolling friction, cohesion); `surface_energy = 0` |
| `add_material_full` | adds `surface_energy` |
| `add_material_extended` | adds twisting friction and Hooke linear stiffnesses `kn`/`kt` |
| `add_material_with_sds` | adds the SDS rolling/twisting spring–dashpot parameters (the full constructor) |

## Mixing rules

Phase 2 mixes per-material values into per-pair `*_ij` tables. The rule depends
on what the quantity is:

| Mixed quantity (per-pair) | Rule |
|---|---|
| `friction_ij`, `rolling_friction_ij`, `twisting_friction_ij` | **Geometric mean** of the per-material values |
| `cohesion_energy_ij`, `surface_energy_ij` | Geometric mean |
| `rolling_damping_ij`, `twisting_damping_ij` (SDS) | Geometric mean |
| `kn_ij`, `kt_ij` (Hooke) | **Harmonic mean** `2kᵢkⱼ/(kᵢ+kⱼ)` |
| `rolling_stiffness_ij`, `twisting_stiffness_ij` (SDS) | Harmonic mean (or one-sided if only one material sets it) |
| `e_eff_ij` (Hertz E*) | `1 / ((1−νᵢ²)/Eᵢ + (1−νⱼ²)/Eⱼ)` |
| `g_eff_ij` (Mindlin G*) | Effective shear modulus from per-material E, ν |
| `beta_ij` (damping) | COR inversion — see below |

The intuition: **friction-like** and **energy-like** quantities mix
geometrically (a smooth material against a rough one lands in between), while
**stiffnesses** mix harmonically (a soft material against a stiff one is governed
by the soft one, like springs in series). The effective moduli come from the
closed-form Hertz/Mindlin contact-of-two-elastic-bodies result.

## Restitution → damping

`restitution` is stored as the **target coefficient of restitution (COR)**, not
a damping ratio. In phase 2, `beta_ij` is computed by **inverting the exact
head-on Hertz collision** — a bisection on the monotone `COR(β)` curve of the
head-on Hertz model. This makes the **input restitution the realized COR** of a
binary collision, so DIRT's shear and cooling results land on the same
kinetic-theory line as LAMMPS/LIGGGHTS (`damping coeff_restitution`) rather than
on a shifted "realized-e" line.

The two models invert differently:

- **Hooke (linear):** `β = −ln(e) / √(π² + ln²e)` — exact for a
  constant-stiffness spring-dashpot, velocity-independent.
- **Hertz (nonlinear):** a numerically exact inversion of the Hertz-collision
  `COR(β)` curve. An older Tsuji *polynomial* fit realized a COR *above* nominal
  (e.g. 0.95 → 0.965); the exact inversion removes that bias so input `e` equals
  realized COR.

This matters because near the elastic limit the granular temperature scales as
`T* ∝ 1/(1 − e²)`, so any input-vs-realized gap throws the stress off. See the
[validation chapter](validation.md) for the measured COR caveat.

## Which model reads which field

`build_pair_tables()` produces these per-pair tables; each model branch reads a
subset:

| Model branch | `MaterialTable` inputs (per-pair) |
|---|---|
| Hertz normal | `e_eff_ij` (E*), `beta_ij` (from `restitution`) |
| Hooke normal | `kn_ij` (harmonic mean of per-material `kn`), `beta_ij` |
| Mindlin tangential | `g_eff_ij` (G*), `friction_ij` (μ), `beta_ij` |
| Hooke tangential | `kt_ij`, `friction_ij` |
| Rolling (constant / SDS) | `rolling_friction_ij`, `rolling_stiffness_ij`, `rolling_damping_ij` |
| Twisting (constant / SDS) | `twisting_friction_ij`, `twisting_stiffness_ij`, `twisting_damping_ij` |
| JKR / DMT adhesion (Hertz only) | `surface_energy_ij`, `adhesion_model` |
| SJKR cohesion | `cohesion_energy_ij` |

## Config-error convention

`add_material*` validates physically inconsistent input — for example, setting
both `cohesion_energy` and `surface_energy` on one material — by printing an
`ERROR:` line to stderr and **exiting the process** (it does not return a
`Result`). This is deliberate: a malformed material table is a config bug that
should stop the run immediately and identically on every MPI rank, rather than
propagate a half-built table. Walls use the same fail-fast convention.

## See also

- [Configuration Reference](config.md) — the full `[dem]` / `[[dem.materials]]`
  TOML schema.
- [Contact: Hertz–Mindlin](../physics/contact.md) — the force laws these
  parameters feed.
- [Particle Insertion](insertion.md) — how particles acquire a material.
