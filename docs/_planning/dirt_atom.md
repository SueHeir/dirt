# dirt_atom — Documentation Planning

Crate: `dirt_atom`
Files: `crates/dirt_atom/src/lib.rs`, `crates/dirt_atom/src/insert.rs`, `crates/dirt_atom/src/radius.rs`

---

## Purpose

`dirt_atom` is the per-atom DEM data and material property layer. It provides:

- **`DemAtom`** — the SOIl `AtomData` extension that carries every per-sphere scalar and vector DEM needs (radius, density, inverse inertia, orientation quaternion, angular velocity, angular momentum, torque, rigid-body id). Registered via the `AtomData` derive macro; field attributes (`#[forward]`, `#[reverse]`, `#[zero]`) drive MPI pack/unpack and zeroing without any hand-written code.
- **`MaterialTable`** — a named-material store plus precomputed `N×N` per-pair mixing tables (`e_eff_ij`, `g_eff_ij`, `beta_ij`, `friction_ij`, …). The contact-force kernels in `dirt_granular` read only the `*_ij` tables, never the per-material rows.
- **`DemAtomPlugin`** — registers `DemAtom`, builds `MaterialTable` from `[dem]` / `[[dem.materials]]` TOML config, sets `Atom::ntypes`.
- **`DemAtomInsertPlugin`** — setup-time random/file insertion, runtime rate-based trickle insertion; all three modes use the **born-in-owner** parallel model.
- **`RadiusSpec`** — `Fixed(f64)` or `Distribution(RadiusDistribution)` (uniform, gaussian, lognormal, discrete).
- Free function **`same_body(dem, i, j) -> bool`** — clump contact-exclusion predicate.

---

## Public surface to document

### `DemAtom` fields and `AtomData` attributes

Defined at `lib.rs:719–746`. All fields are `Vec<T>`, indexed by local atom index.

| Field | Type | Attribute | Meaning |
|---|---|---|---|
| `radius` | `Vec<f64>` | _(none)_ | Particle radius [m]; set at insertion, used in contact detection |
| `density` | `Vec<f64>` | _(none)_ | Material density [kg/m³]; used to compute mass and inertia |
| `inv_inertia` | `Vec<f64>` | _(none)_ | Precomputed `1 / (2/5 · m · r²)` [1/(kg·m²)]; avoids divide-per-step |
| `quaternion` | `Vec<[f64;4]>` | _(none)_ | Orientation `[w,x,y,z]`; inert for spheres unless `track_orientation = true` |
| `omega` | `Vec<[f64;3]>` | **`#[forward]`** | Angular velocity [rad/s]; forwarded owner→ghost each step |
| `ang_mom` | `Vec<[f64;3]>` | _(none)_ | Angular momentum [kg·m²/s]; migrated on ownership transfer |
| `torque` | `Vec<[f64;3]>` | **`#[reverse]` `#[zero]`** | Torque [N·m]; ghost→owner accumulation; zeroed before each force computation |
| `body_id` | `Vec<f64>` | **`#[forward]`** | Clump body id; `0.0` = independent sphere; positive and equal → same rigid body |

`#[forward]` means forwarded owner-to-ghost each step. `#[reverse]` + `#[zero]` means ghost contributions are summed back to the owner, then zeroed before the next force step. Fields without attributes are communicated only on atom migration (ownership transfer). Source: `lib.rs:719–746`, attribute semantics in `lib.rs:36–43`.

### `MaterialTable` API

Source: `lib.rs:311–705`.

**Construction ladder** (shortest to fullest; each delegates to the next):

| Method | New vs. prior |
|---|---|
| `add_material(name, E, ν, e, μ, μ_roll, γ_coh)` | basics; `surface_energy = 0` |
| `add_material_full(…, γ_surf)` | adds JKR/DMT surface energy |
| `add_material_extended(…, μ_twist, kn, kt)` | adds twisting friction and Hooke stiffnesses |
| `add_material_with_sds(…, k_roll, γ_roll, k_twist, γ_twist)` | all SDS spring–dashpot params |

All return a `u32` material index. `find_material(name) -> Option<u32>` for name lookup.

**Commit phase:** `build_pair_tables()` must be called once after all materials are registered and before any contact evaluation; populates all `*_ij` fields. Called automatically by `DemAtomPlugin::build`.

**Per-pair tables** (all `Vec<Vec<f64>>`, indexed `[i][j]`):

| Field | Mixing rule |
|---|---|
| `beta_ij` | Exact COR inversion (Hertz) or analytic formula (Hooke) |
| `friction_ij` | Geometric mean of per-material `friction` |
| `rolling_friction_ij` | Geometric mean |
| `cohesion_energy_ij` | Geometric mean |
| `surface_energy_ij` | Geometric mean (JKR/DMT) |
| `twisting_friction_ij` | Geometric mean |
| `e_eff_ij` | Hertz effective Young's modulus: `1/((1−νᵢ²)/Eᵢ + (1−νⱼ²)/Eⱼ)` |
| `g_eff_ij` | Mindlin effective shear modulus |
| `kn_ij`, `kt_ij` | Harmonic mean of per-material `kn`/`kt` |
| `rolling_stiffness_ij`, `twisting_stiffness_ij` | Harmonic mean (or one-sided if only one is non-zero) |
| `rolling_damping_ij`, `twisting_damping_ij` | Geometric mean |

Source: `lib.rs:563–704`.

### Insertion plugin

`DemAtomInsertPlugin` (`insert.rs:308–360`): adds `dem_insert_atoms` setup system (after `domain_read_input`) and `dem_rate_insert` update system (`PreInitialIntegration`). Also triggers `calculate_delta_time` (Rayleigh-wave criterion, `insert.rs:1517–1561`).

`insert_single_particle` (`insert.rs:364–401`): canonical single-atom append; computes `mass = ρ · (4/3)π r³`, `inv_inertia = 1/(0.4 · m · r²)`, initializes quaternion to identity, all rotational fields to zero, `body_id` to 0.

`same_body(dem, i, j)` (`lib.rs:772–776`): returns `true` iff `body_id[i] > 0`, `body_id[j] > 0`, and `|body_id[i] − body_id[j]| < 0.5`. Used by clump contact exclusion.

---

## Config / TOML schema

### `[dem]` top-level

Source: `lib.rs:204–239`, `#[serde(deny_unknown_fields)]`.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `contact_model` | `String` | `"hertz"` | `"hertz"` or `"hooke"` |
| `adhesion_model` | `String` | `"jkr"` | `"jkr"` or `"dmt"`; used when any material has `surface_energy > 0` |
| `rolling_model` | `String` | `"constant"` | `"constant"` or `"sds"` |
| `twisting_model` | `String` | `"constant"` | `"constant"` or `"sds"` |
| `track_orientation` | `bool` | `false` | Integrate per-sphere quaternion; pure overhead for spheres (see `lib.rs:228–238`) |
| `materials` | `Vec<MaterialConfig>` | _(empty)_ | Via `[[dem.materials]]` sub-tables |

### `[[dem.materials]]`

Source: `lib.rs:144–190`, `#[serde(deny_unknown_fields)]`.

| Key | Required | Default | Meaning |
|---|---|---|---|
| `name` | yes | — | String name referenced by particles and walls |
| `youngs_mod` | yes | — | Young's modulus E [Pa] |
| `poisson_ratio` | yes | — | Poisson's ratio ν (0–0.5) |
| `restitution` | yes | — | Target coefficient of restitution (0–1) |
| `friction` | no | `0.4` | Sliding (Coulomb) friction coefficient |
| `rolling_friction` | no | `0.0` | Rolling resistance coefficient (0 = disabled) |
| `twisting_friction` | no | `0.0` | Twisting friction coefficient (0 = disabled) |
| `cohesion_energy` | no | `0.0` | SJKR cohesion energy density [J/m²] |
| `surface_energy` | no | `0.0` | JKR/DMT surface energy [J/m²]; mutually exclusive with `cohesion_energy` |
| `kn` | no | `0.0` | Linear normal stiffness for Hooke model [N/m]; 0 = use Hertz |
| `kt` | no | `0.0` | Linear tangential stiffness for Hooke model [N/m]; 0 = use Mindlin |
| `rolling_stiffness` | no | `0.0` | SDS rolling spring stiffness [N·m/rad] |
| `rolling_damping` | no | `0.0` | SDS rolling viscous damping |
| `twisting_stiffness` | no | `0.0` | SDS twisting spring stiffness [N·m/rad] |
| `twisting_damping` | no | `0.0` | SDS twisting viscous damping |

### `[[particles.insert]]`

Source: `insert.rs:96–166`, `InsertConfig` (no `deny_unknown_fields`).

Three modes; selected implicitly by field presence:

**Random mode** (default, `source = "random"` or omitted):

| Key | Required | Default | Meaning |
|---|---|---|---|
| `material` | yes | — | Name of a `[[dem.materials]]` entry |
| `count` | yes | — | Number of particles to insert at setup time |
| `radius` | yes | — | `RadiusSpec`: fixed `f64` or distribution table |
| `density` | yes | — | Particle density [kg/m³] |
| `velocity` | no | `0.0` | Gaussian random speed magnitude [m/s] |
| `velocity_x/y/z` | no | `0.0` | Directional velocity offset [m/s] (additive with random) |
| `region` | no | domain inset by max radius | Insertion region (see Region types below) |
| `seed` | no | `0` | RNG seed; determines the global packing |

**Rate-based mode** (presence of `rate` field triggers this mode):

| Key | Required | Default | Meaning |
|---|---|---|---|
| `rate` | yes | — | Particles to insert per interval |
| `rate_interval` | no | `1` | Insert every N timesteps |
| `rate_start` | no | `0` | First timestep for insertion |
| `rate_end` | no | _(never)_ | Last timestep for insertion |
| `rate_limit` | no | _(unlimited)_ | Maximum total particles inserted |

Rate mode also requires `material`, `radius`, `density`. The step-derived RNG seed (`config_seed ^ step_hash ^ entry_hash`) keeps candidates deterministic per-step per-entry across ranks. New atoms may be born overlapping existing particles (the overlap scratch for rate insertion covers only the new atoms, not existing local atoms — see `insert.rs:1431–1440` for the rationale).

**File-based mode** (`source = "file"`):

| Key | Required | Default | Meaning |
|---|---|---|---|
| `source` | yes | — | Must be `"file"` |
| `file` | yes | — | Path to input file |
| `format` | yes | — | `"csv"`, `"lammps_dump"`, or `"lammps_data"` |
| `material` | yes | — | Default material for atoms without a type_map entry |
| `density` | yes (csv/lammps_dump) | — | Particle density; embedded in sphere/bpm/sphere data files |
| `radius` | optional | — | Default radius; overridden by file column or file data |
| `columns` | no | `{x=0,y=1,z=2,radius=3}` | Column-index map for CSV |
| `type_map` | no | — | `{ 1 = "glass", 2 = "steel" }` — maps integer file types to materials |
| `atom_style` | no | auto-detected | LAMMPS data atom style: `"atomic"`, `"sphere"`, `"bpm/sphere"` |

LAMMPS data `sphere` and `bpm/sphere` styles encode `diameter density x y z` (7 columns); `atomic` style needs `radius` and `density` from config. Style is auto-detected from the `Atoms # style` comment in the data file (`insert.rs:1082–1100`).

**Region types** (from SOIL `Region`): `{ type = "block", min = [x,y,z], max = [x,y,z] }` and `{ type = "cylinder", center = [x,y], radius = r, axis = "z", lo = z1, hi = z2 }` (see `insert.rs:1803–1813` test for cylinder example).

**`RadiusSpec`** in `radius` field (`radius.rs:18–46`):

| Variant | TOML example |
|---|---|
| `Fixed` | `radius = 0.001` |
| `Uniform` | `radius = { distribution = "uniform", min = 0.0008, max = 0.0012 }` |
| `Gaussian` | `radius = { distribution = "gaussian", mean = 0.001, std = 0.0001 }` |
| `Lognormal` | `radius = { distribution = "lognormal", mean = 0.001, std = 0.0001 }` — `mean`/`std` are of the radius, not of the underlying normal |
| `Discrete` | `radius = { distribution = "discrete", values = [0.001, 0.0015], weights = [0.7, 0.3] }` |

Gaussian clamps samples to `max(value, 1e-15)` to prevent non-physical negative radii (`radius.rs:73`). Lognormal converts actual-distribution mean/std to the underlying normal parameters (`radius.rs:77–80`).

---

## Key behaviors, invariants, and gotchas

### 1. Two-phase MaterialTable build contract (`lib.rs:250–263`)

All `*_ij` pair tables are empty `Vec::new()` until `build_pair_tables()` is called. Indexing any `*_ij` table before phase 2 panics with an out-of-bounds error. `DemAtomPlugin::build` (`lib.rs:812–853`) always calls this after loading config; tests and standalone tools must call it manually.

### 2. Exact COR inversion for Hertz (`lib.rs:84–132`, `lib.rs:606–624`)

The input `restitution` value in `[[dem.materials]]` is the **realized COR**, not a damping ratio. For Hertz contact, `build_pair_tables` calls `hertz_beta_for_cor(e_ij)` — a 60-iteration bisection on the monotone `COR(β)` curve obtained by RK4 integration of the dimensionless Hertz collision (`lib.rs:84–109`). This is velocity-independent (Tsuji scaling) so one integration maps the full β↔COR curve. The old Tsuji polynomial fit realised COR above nominal (e.g. 0.95 → 0.965); the exact inversion removes this bias. For Hooke contact, the analytic formula `β = −ln(e)/√(π²+ln²e)` is used (`lib.rs:622–624`). Caveat: `config-anatomy.md` line 69–73 still mentions the old note that "restitution is not the realized COR" — this is now **stale documentation** for DIRT; the exact inversion makes it true.

### 3. `cohesion_energy` and `surface_energy` are mutually exclusive (`lib.rs:528–533`)

Setting both `> 0` on one material causes an immediate `eprintln! + process::exit(1)`. JKR/DMT adhesion (`surface_energy`) is only honored under Hertz contact; under Hooke contact it is silently ignored — only SJKR cohesion (`cohesion_energy`) applies (`lib.rs:307–311`).

### 4. `track_orientation = false` is correct for pure-sphere runs (`lib.rs:228–238`)

A sphere's orientation quaternion never enters any contact force law (forces depend on `ω`, not absolute orientation). Integrating it is pure overhead. Only needed if something downstream reads orientation (e.g. surface-marker viz). Non-spherical bodies store orientation in `BodyData`, not here.

### 5. Born-in-owner parallel insertion (`insert.rs:1–30`)

Every rank seeds its RNG from the same `seed` and generates the bit-identical candidate stream. Each rank materialises only atoms whose position satisfies `low ≤ pos < high` (half-open, matching `Domain::exchange()` ownership). Overlap detection uses a replicated global spatial hash (`SpatialHash`, `insert.rs:230–303`) — O(1) per candidate via a 3×3×3 cell neighborhood. This makes the global packing identical whether you run on 1 or 64 ranks.

### 6. Rate insertion overlap scratch is local-to-step-only (`insert.rs:1431–1440`)

For rate-based insertion the overlap scratch covers only the **new atoms this step**, not existing local atoms (which differ per rank and would desync accept/reject across ranks). New atoms may therefore be born overlapping existing particles. This is documented as intentional — rate-insert regions are normally placed in free space.

### 7. Spatial hash PBC correction for small boxes (`insert.rs:286–299`)

When any periodic axis is smaller than `3 × cell_size + 2.2 × max_radius`, the standard 3×3×3 neighborhood may miss periodic images; the code falls back to a brute-force scan of all atoms using minimum-image distances.

### 8. Auto-computed timestep (`insert.rs:1517–1561`)

`calculate_delta_time` runs at `PostSetup` (after insertion). It estimates the Rayleigh wave transit time per particle: `dt_R = π·r / α · √(ρ/G)`, `α ≈ 0.1631·ν + 0.8766`, `G = E/(2(1+ν))`. Takes the global minimum (MPI all-reduce) and uses `0.15 × dt_min`. This is overridden by an explicit `dt` in the `[[run]]` block.

### 9. Material config-error convention (`lib.rs:299–304`)

`add_material*` validators print `ERROR:` to stderr and call `std::process::exit(1)` — no `Result`, no panic. Deliberate: a malformed material table is a config bug that should stop all MPI ranks immediately rather than propagate a half-built table.

### 10. `same_body` uses floating-point equality via `< 0.5` tolerance (`lib.rs:772–776`)

`body_id` is stored as `f64` for MPI communication convenience. Two atoms are in the same body if `|bi − bj| < 0.5`, which is fine for integer body ids but would fail for fractional ids. Any code generating `body_id` values must ensure they are integers.

### 11. Lognormal distribution parameterisation (`radius.rs:76–80`)

The `mean` and `std` in `{ distribution = "lognormal", mean = …, std = … }` are the desired mean and standard deviation **of the radius distribution** (not the underlying normal's μ and σ). The code converts: `σ² = ln(1 + (std/mean)²)`, `μ = ln(mean) − σ²/2`.

### 12. `atom_type` encodes material index

`Atom::atom_type[i]` is set to the material index (returned by `resolve_material`) at insertion time. This is the coupling point to `MaterialTable`; contact code reads `atom.atom_type[i]` and `atom.atom_type[j]` to look up `*_ij[i][j]`.

---

## Tutorial outline

For `reference/materials.md` and `reference/insertion.md`:

**Materials page** (largely drafted, `docs/src/reference/materials.md`):
- Correct the stale "realized COR ≠ input COR" note inherited from `config-anatomy.md` (see gotcha 2 above).
- Add a Rust code block showing `MaterialTable::new()` + `add_material_with_sds` + `build_pair_tables()` with SDS parameters (currently only the short form is shown).
- Document JKR vs. DMT model selection and the Hertz-only restriction for JKR/DMT.
- Add the "which field is consulted only under Hertz, which under Hooke" table already in the file but verify against `dirt_granular` implementations.
- Add `track_orientation` guidance as a sidebar.

**Insertion page** (largely drafted, `docs/src/reference/insertion.md`):
- Expand the radius distribution section to a full table with all four variants and their parameter meanings (including the lognormal parameterisation note).
- Add rate-based insertion: full field table (`rate`, `rate_interval`, `rate_start`, `rate_end`, `rate_limit`), the step-derived seed derivation, and the overlap-scratch caveat.
- Expand file-based section: distinguish `lammps_dump` vs. `lammps_data`; document `atom_style` auto-detection; document `type_map` + fallback material.
- Add a note on the auto-computed timestep and when to override it with explicit `dt`.
- Region types: document `block` and `cylinder` variants with their fields.

---

## Doc gaps

1. **`config-anatomy.md` line 69–73** — the "input restitution is not the realized COR" note is now wrong for DIRT. It describes the old Tsuji behaviour. Should be corrected to reflect exact COR inversion.
2. **`reference/insertion.md`** — rate-based insertion mode is described only superficially (one sentence). Missing: full field table, overlap-scratch caveat, RNG seeding per step.
3. **`reference/insertion.md`** — `lammps_data` format is not mentioned (only `lammps_dump` and CSV). Auto-detected `atom_style`, `bpm/sphere` style, and the per-atom density/radius in sphere-style files are undocumented.
4. **`reference/insertion.md`** — Region types are not enumerated. The `cylinder` region variant appears in a test (`insert.rs:1803`) but is not in the docs.
5. **`reference/insertion.md`** — `type_map` and the fallback-material pattern for mixed-type files are undocumented.
6. **`reference/materials.md`** — JKR vs. DMT: no guidance on Tabor parameter or when to choose each (the code comment at `lib.rs:215–220` has the physics; it should migrate to docs).
7. **`reference/materials.md`** — No example with `rolling_model = "sds"` and the full SDS field set.
8. **`DemAtomPlugin::default_config`** (`lib.rs:790–808`) and **`DemAtomInsertPlugin::default_config`** (`insert.rs:313–344`) contain the canonical commented examples. These should stay synchronized with the reference docs.
9. **`same_body` floating-point gotcha** (gotcha 10) is not mentioned anywhere in the docs; clumps.md only describes the predicate's semantics.

---

## Suggested placement

- **`reference/materials.md`** — primary home for `MaterialTable`, `DemConfig`, `[[dem.materials]]` schema (all keys, mixing rules, restitution→damping, JKR vs. DMT, SDS). Page exists and has solid bones; needs the gaps above filled.
- **`reference/insertion.md`** — primary home for `[[particles.insert]]` full schema, `RadiusSpec` distributions, region types, rate-mode fields, file formats, born-in-owner model, auto-timestep. Page exists; needs rate and file sections expanded significantly.
- **`reference/config.md`** — aggregate reference that links to the above two; already correctly cross-references both. Should be kept as a summary, not duplicated detail.
- **`physics/clumps.md`** — already cross-references `same_body`; the floating-point tolerance gotcha could be added as a warning callout there.
- **Rust API docs (`cargo doc`)** — `MaterialTable`, `DemAtom`, `InsertConfig`, `RadiusSpec`, `RadiusDistribution` all have inline doc comments sufficient for `cargo doc`; no standalone API doc page needed in mdBook.

---

## File:line reference index

| Topic | Location |
|---|---|
| `DemAtom` struct definition | `lib.rs:719–746` |
| `#[forward]` / `#[reverse]` / `#[zero]` semantics | `lib.rs:36–43` |
| `same_body` predicate | `lib.rs:772–776` |
| `MaterialConfig` struct (all TOML keys) | `lib.rs:144–190` |
| `DemConfig` struct | `lib.rs:204–239` |
| `MaterialTable` struct fields | `lib.rs:311–367` |
| `add_material*` ladder | `lib.rs:447–551` |
| `build_pair_tables` mixing rules | `lib.rs:563–704` |
| Hertz COR inversion (RK4 + bisection) | `lib.rs:84–132` |
| `DemAtomPlugin::build` | `lib.rs:812–853` |
| `InsertConfig` struct (all TOML keys) | `insert.rs:96–166` |
| Born-in-owner model explanation | `insert.rs:1–30` |
| `insert_single_particle` (mass + inertia formula) | `insert.rs:364–401` |
| `owns_position` half-open interval | `insert.rs:410–412` |
| `SpatialHash` (3×3×3 + PBC fallback) | `insert.rs:230–303` |
| Rate-based insertion loop | `insert.rs:1341–1512` |
| Rate insertion overlap-scratch rationale | `insert.rs:1431–1440` |
| `calculate_delta_time` (Rayleigh criterion) | `insert.rs:1517–1561` |
| `RadiusSpec` / `RadiusDistribution` enum | `radius.rs:16–46` |
| Lognormal parameter conversion | `radius.rs:77–80` |
| Gaussian positive clamp | `radius.rs:73` |
| `cohesion_energy` + `surface_energy` mutual exclusion | `lib.rs:528–533` |
| JKR/DMT Hertz-only note | `lib.rs:307–311` |
| `track_orientation` rationale | `lib.rs:228–238` |
| Stale COR note in existing docs | `docs/src/getting-started/config-anatomy.md:69–73` |
