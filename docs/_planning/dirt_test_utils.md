# Planning doc: `dirt_test_utils`

_Target chapter: `docs/src/reference/writing-tests.md`_

---

## Purpose

`dirt_test_utils` is a narrow, deliberately dependency-light crate (only `soil_core`
+ `dirt_atom`) that builds minimal, valid test fixtures so that unit tests in every
DIRT consumer crate can focus on the logic under test rather than boilerplate setup.
It deliberately does **not** depend on `grass_app` or `grass_scheduler`, keeping
compilation fast and the crate boundary clean.

Source: `crates/dirt_test_utils/src/lib.rs:1–112` (module-level doc comment).

---

## Public surface to document

Five public free functions; no public types or traits.

### `make_atoms(n: usize) -> Atom`
`lib.rs:131–140`

Constructs an `Atom` with `n` test particles arranged along the x-axis at
`(i as f64, 0.0, 0.0)`, each with radius `0.5`, mass `1.0`, and tag `i as u32`.
Sets `atom.nlocal = n`, `atom.natoms = n`, and `atom.dt = 0.001`.

Used by: `dirt_fixes` tests (`lib.rs:668`) for non-DEM (fix/force-accumulation) tests
where particle geometry and DEM fields are irrelevant.

### `make_group_registry(name: &str, mask: Vec<bool>) -> GroupRegistry`
`lib.rs:156–171`

Builds a `GroupRegistry` containing exactly one `Group` whose membership is given by
`mask` (true = member). The count is derived automatically via `mask.iter().filter(…)`.
Only one group is created; tests that need multiple groups must call this helper
multiple times and merge, or build `GroupRegistry` directly.

Used by: `dirt_fixes` tests (`lib.rs:668`) to target a subset of atoms with a fix.

### `make_single_comm() -> CommResource`
`lib.rs:185–187`

Wraps `SingleProcessComm::new()` in a `CommResource`. No MPI initialization needed.
Suitable for all non-MPI unit tests.

### `push_dem_test_atom(atom: &mut Atom, dem: &mut DemAtom, tag: u32, pos: [f64; 3], radius: f64)`
`lib.rs:212–230`

Appends one solid sphere to the parallel `Atom` and `DemAtom` arrays:
- density: 2500 kg/m³ (hardcoded)
- mass: `ρ · (4/3)πr³`
- `inv_inertia`: `1 / (0.4 · m · r²)` (solid sphere)
- quaternion: `[1, 0, 0, 0]` (identity)
- omega, ang_mom, torque: zero
- body_id: `0.0`

Does **not** touch `atom.nlocal` or `atom.natoms`.

Used extensively in `dirt_granular/contact.rs` (via a local wrapper
`push_test_atom_with_history`, `contact.rs:1000–1010` that also pushes a per-atom
entry into `ContactHistoryStore`), in `dirt_wall/lib.rs:1695`, in
`dirt_bond/lib.rs:1377`, and in `dirt_granular/rotational.rs:132`.

### `make_material_table() -> MaterialTable`
`lib.rs:257–262`

Returns a single-material `MaterialTable` named `"glass"` with:

| Property        | Value  |
|-----------------|--------|
| Young's modulus | 8.7 GPa |
| Poisson ratio   | 0.3    |
| Restitution     | 0.95   |
| Friction        | 0.4    |
| Rolling friction| 0.0    |
| Cohesion energy | 0.0    |

Calls `mt.build_pair_tables()` before returning, so it is immediately usable as an
`App` resource in contact-law tests.

Tests that need non-default material parameters define private helpers locally — e.g.,
`make_material_table_cohesion`, `make_material_table_jkr`, `make_material_table_hooke`,
`make_material_table_rolling`, `make_material_table_sds_rolling`, etc. in
`dirt_granular/contact.rs:1149–1670`. This is intentional; the shared helper is the
neutral "glass" baseline.

---

## Config / TOML schema

None. This crate provides no plugin, no config key, and no TOML schema. There is
nothing to document here.

---

## Key behaviors, invariants, and gotchas

### 1. `push_dem_test_atom` does NOT set `nlocal` / `natoms`
`lib.rs:93–97`

The helper only appends to the backing `Vec`s. `atom.nlocal` and `atom.natoms` remain
at their default of `0`. Force systems iterate `0..nlocal`, so forgetting to set these
after the last `push_*` call produces a silent no-op: the system loops over an empty
range, no assertion fires, and the test passes with zero force applied. Set both
**after** the last `push_dem_test_atom` call:

```rust
atom.nlocal = 2;
atom.natoms = 2;
```

### 2. DEM contacts need `dt ≈ 1e-7`
`lib.rs:98–101`; confirmed in `contact.rs:1019`, `rotational.rs:130`

`make_atoms` sets `dt = 0.001`, which is appropriate for non-DEM fixture tests but
four orders of magnitude too large for stiff Hertz contacts. For any contact or bond
test, override before pushing atoms:

```rust
atom.dt = 1e-7;
```

### 3. `push_test_atom_with_history` is a local pattern, not in the crate
`contact.rs:1000–1010`

`dirt_granular` defines a private wrapper in its `#[cfg(test)]` module that calls
`push_dem_test_atom` and additionally pushes an empty per-atom slot into
`ContactHistoryStore`. This is not part of `dirt_test_utils` because
`ContactHistoryStore` is a `dirt_granular`-internal type — a deliberate layering
choice. Any crate with per-contact or per-bond history must define an equivalent
local wrapper.

### 4. `make_material_table` returns the neutral glass baseline; extend locally
`lib.rs:257–262`; see `contact.rs:1149–2030` for the pattern.

Cohesion, JKR adhesion, rolling-resistance, SDS, Hooke, and twisting variants are
built locally in each test module. The shared helper covers only the no-special-
behaviour case.

### 5. Density is hardcoded at 2500 kg/m³
`lib.rs:219`

There is no way to override density through `push_dem_test_atom`. If a test requires
a different density it must call `atom.push_test_atom` and populate `DemAtom` fields
directly.

### 6. Neighbor list must be built by hand (CSR)
`lib.rs:107` (what the crate does NOT provide)

There is no neighbor-list builder. The standard pattern is:

```rust
let mut neighbor = Neighbor::new();
neighbor.neighbor_offsets = vec![0, 1, 1]; // atom 0 has one neighbor; atom 1 has none
neighbor.neighbor_indices = vec![1];        // that neighbor is atom 1
```

Only one direction is needed (particle 0 → 1); `dirt_granular` contact systems use
Newton's-third-law pairing internally.

### 7. `AtomDataRegistry::register` order is significant
Inferred from `contact.rs:1036–1038`

`DemAtom` and `ContactHistoryStore` must both be registered into `AtomDataRegistry`
before it is added to the `App`. Systems retrieve them by type via the registry, so
omitting either registration produces a panic at system run time.

---

## Tutorial outline — writing a unit test for a new force term

Intended section in `docs/src/reference/writing-tests.md`.

1. **Add `dirt_test_utils` as a dev-dependency** in the consumer crate's `Cargo.toml`.
2. **Create an `Atom` + `DemAtom`**, set `atom.dt = 1e-7`.
3. **Push atoms** with `push_dem_test_atom`; choose positions so overlap is deliberate
   (e.g., two spheres of radius 0.001 m separated by 0.0019 m → 10 % overlap).
4. **Set `atom.nlocal` and `atom.natoms`** after the last push.
5. **Build `Neighbor` by hand** (CSR offsets and indices for the desired contact pair).
6. **Register per-atom data** into `AtomDataRegistry`. If your force law uses
   a history store, define a local `push_test_atom_with_history` wrapper that also
   pushes into that store.
7. **Assemble the App**: add `atom`, `neighbor`, `registry`, `make_material_table()`
   (or a local variant), schedule your system at `ParticleSimScheduleSet::Force`,
   call `organize_systems()`, then `run()`.
8. **Assert**: retrieve `Atom` via `app.get_resource_ref::<Atom>().unwrap()` and check
   `atom.force[i]`. For action–reaction tests, verify
   `(atom.force[0][0] + atom.force[1][0]).abs() < 1e-10`.

Reference: `contact.rs:1012–1051` (`fused_contact_repulsive_for_overlap`) is the
canonical minimal example.

---

## Doc gaps

1. **No API-level docs on `make_single_comm`** beyond a single-line comment — the
   chapter should state when tests need `CommResource` vs. when they can omit it
   entirely. (`lib.rs:185–187`)

2. **Density is not configurable.** The crate doc mentions "density 2500 kg/m³" but
   no doc explains what to do if a test needs a different density. The chapter should
   include a one-paragraph note pointing to `Atom::push_test_atom` + manual `DemAtom`
   field population as the escape hatch.

3. **`push_test_atom_with_history` pattern is undocumented.** It is reinvented in
   every crate that uses contact history. The chapter should describe the pattern and
   explain why it is local (layering constraint).

4. **No explanation of when `make_group_registry` is needed vs. irrelevant** for
   DEM tests. The `dirt_fixes` use is clear, but new contributors may add it
   unnecessarily to contact tests.

5. **The `AtomDataRegistry::register` ordering requirement** is not documented in the
   crate; it should appear as a footgun callout in the chapter.

6. **No cross-reference from material-variant helpers** — the chapter should direct
   readers to `contact.rs` for examples of `make_material_table_jkr` etc. as the
   pattern to follow.

---

## Suggested placement

This content belongs in `docs/src/reference/writing-tests.md`, which already exists
and has a skeleton. The chapter is already listed in `SUMMARY.md:30`. The existing
file (`writing-tests.md`) covers the basics; the gaps above (items 1–6) should be
added as additional sections or callout boxes, in roughly this order:

1. Quick start (already present — verify against `lib.rs:9–18`)
2. The full force-law test pattern (already present — verify against `contact.rs:1012`)
3. Two footguns (already present — expand with the registry ordering footgun)
4. When to use `make_group_registry` (new)
5. Customising materials beyond the glass baseline (new)
6. The `push_test_atom_with_history` local-wrapper pattern (new)
7. Density / geometry overrides (new)
8. What this crate does NOT provide (already present)
