# Planning: dirt_clump documentation

**Crate:** `dirt_clump`
**Target page:** `docs/src/physics/clumps.md` (already exists; this plan informs a rewrite/expansion)

---

## Purpose

`dirt_clump` adds rigid multisphere/clump composites to DIRT, enabling
non-spherical DEM particles. A clump is a set of overlapping spheres that are
permanently fused into one rigid body: each sub-sphere participates in contact
detection as an ordinary atom, but forces are aggregated to the parent body,
which integrates full 6-DOF (translational + rotational Euler) dynamics. The
sub-sphere positions and velocities are derived back from body state after
each step. This is the standard "multisphere" approach, modeled after
LIGGGHTS `FixMultisphere`.

---

## Public surface to document

### Plugin
- **`ClumpPlugin`** (`lib.rs:308`) — the entry point. Registers `ClumpAtom`
  per-atom data, `ClumpRegistry`, `MultisphereBodyStore`, and all systems
  (setup, exchange, integration, force aggregation, position update, lost-atom
  diagnostic). Declares dependency on `DemAtomPlugin`. Adding this plugin to
  `App` is the only user-facing step to enable clumps.

### Configuration structs (TOML → Rust)
- **`ClumpSphereConfig`** (`lib.rs:143`) — one sphere: `offset: [f64;3]`,
  `radius: f64`.
- **`ClumpDef`** (`lib.rs:153`) — a named clump type: `name: String`,
  `spheres: Vec<ClumpSphereConfig>`. Loaded from `[[clump.definitions]]`.
- **`ClumpInsertConfig`** (`lib.rs:162`) — one insertion block:
  `definition`, `count`, `density`, `material`, optional `velocity`, optional
  `region`. Loaded from `[[clump.insert]]`.
- **`ClumpTopConfig`** (`lib.rs:183`) — top-level `[clump]` container
  (separate from `[dem]` which uses `deny_unknown_fields`).

### Per-atom data
- **`ClumpAtom`** (`lib.rs:198`) — `body_id: Vec<f64>` (0 = not in a clump;
  encoded as f64 for `AtomData` compatibility) and `body_offset: Vec<[f64;3]>`
  (body-frame offset of each sub-sphere, fixed at creation).

### Body types
- **`MultisphereBody`** (`body.rs:21`) — full rigid body state: `com_pos`,
  `com_vel`, `quaternion` [w,x,y,z] (body→space), `omega`, `angmom`,
  `principal_moments`, `principal_axes` (body→principal), `total_mass`,
  `inv_mass`, `force`, `torque`, `image`, `body_offsets`, `sub_sphere_radii`,
  `sub_sphere_tags`. Both `pack`/`unpack` (full MPI exchange) and
  `pack_forward`/`pack_reverse` (ghost comm) are implemented.
- **`MultisphereBodyStore`** (`body.rs:73`) — flat `Vec<MultisphereBody>` plus
  an O(1) ID→index flat map (`map: Vec<usize>`). The map must be regenerated
  with `generate_map()` after any add/remove; called automatically after
  insertion and exchange.
- **`ClumpRegistry`** (`lib.rs:228`) — runtime list of `ClumpDef`s; used by
  the insert system to look up definitions by name.

### Public functions
| Function | File:line | Role |
|---|---|---|
| `same_body(data, i, j)` | `lib.rs:1208` | Contact-exclusion predicate — must be called by every force plugin |
| `is_body_atom(data, i)` | `lib.rs:1219` | True if atom i is a clump sub-sphere |
| `insert_clump(...)` | `lib.rs:1097` | Programmatic clump insertion (used by the config-driven setup system; can be called directly) |
| `compute_inertia_tensor_analytical` | `body.rs:295` | Exact parallel-axis tensor for non-overlapping spheres |
| `compute_inertia_tensor_montecarlo` | `body.rs:330` | Monte Carlo tensor for overlapping spheres (hardcoded 100 000 samples) |
| `has_overlap` | `body.rs:400` | Detect whether any sub-sphere pair overlaps (drives path selection) |
| `diagonalize_inertia` | `body.rs:529` | Jacobi eigdecomposition → principal moments + axes quaternion |
| `jacobi_eigendecomposition` | `body.rs:421` | Raw 3×3 symmetric Jacobi solver (exported) |
| `rotation_matrix_to_quaternion` | `body.rs:493` | Shoemake method: 3×3 rotation → unit quaternion |
| `angmom_to_omega` | `body.rs:191` | L (space) → ω (space), via principal frame rotation |
| `quat_rotate` | `lib.rs:252` | Rotate vector by quaternion [w,x,y,z] |
| `cross` | `lib.rs:271` | 3-vector cross product |
| `compute_clump_inertia` | `lib.rs:289` | **Legacy** — trace/3 scalar moment; kept for backward compat only |

### Integration entry points (public, called by systems)
- `body::integrate_body_initial(body, dt)` — half-kick vel + drift pos + half-kick angmom + Richardson quaternion update + recompute omega. (`body.rs:547`)
- `body::integrate_body_final(body, dt)` — half-kick vel + half-kick angmom + recompute omega; no quaternion update. (`body.rs:593`)

---

## Config / TOML schema

```toml
# Top-level section — separate from [dem] (which denies unknown fields)
[clump]

# Named clump type definitions — at least one sphere required
[[clump.definitions]]
name = "dimer"              # string; referenced by [[clump.insert]].definition
spheres = [
    # offset: body-frame displacement of sub-sphere center from COM [x, y, z] (m)
    # radius: sub-sphere radius (m)
    { offset = [-0.0003, 0.0, 0.0], radius = 0.001 },
    { offset = [ 0.0003, 0.0, 0.0], radius = 0.001 },
]

# Insertion blocks — one block per clump type to insert
[[clump.insert]]
definition = "dimer"        # must match a [[clump.definitions]].name
count      = 100            # number of clumps to insert (u32)
density    = 2500.0         # kg/m³; used for both mass and inertia
material   = "glass"        # must match a [[dem.materials]].name entry
velocity   = 0.5            # optional f64 (m/s); each component drawn uniform in [-v,+v]
                            # omit or null → zero initial velocity
region     = { type = "block", min = [0.001, 0.001, 0.001],
                               max = [0.019, 0.019, 0.019] }
             # optional Region; defaults to domain inset by effective clump radius
```

**Effective clump radius** for overlap checks and default region inset is
`max over spheres of (|offset| + radius)` (`lib.rs:977`). Insertion uses a
5% margin: `min_sep = 2 * eff_radius * 1.05` (`lib.rs:1017`).

Insertion runs only on rank 0, then atoms and bodies are distributed by the
normal MPI exchange (`lib.rs:944`).

---

## Key behaviors, invariants, and gotchas

### 1. Contact-exclusion contract — `same_body` (`lib.rs:1208`)

Sub-spheres of the same body must **never** generate forces against each other.
The crate does **not** filter the neighbor list — it is the responsibility of
every force plugin to call `same_body(&clump_data, i, j)` and skip on `true`.
`dirt_granular` already does this. A custom force plugin that forgets this will
cause self-repulsion and explosive body separation. `same_body` compares
`body_id` as f64 using `(ci - cj).abs() < 0.5`, so the encoding as f64 is
intentional (IDs are small integers, exact in f64).

### 2. Sub-sphere `inv_mass` is zero (`lib.rs:1171`)

Sub-sphere atoms have `atoms.inv_mass[i] = 0.0`. This is intentional: they
are not integrated by the standard Verlet path. All translational dynamics go
through the body. Any custom integrator that tries to advance sub-spheres
independently using `inv_mass` will silently freeze them.

### 3. Angular momentum as the integrated state (`body.rs:547`, `body.rs:593`)

`angmom` (space-frame) is the primary rotational variable, half-kicked by
torque each step. `omega` is derived on demand via `angmom_to_omega`, never
integrated directly. This avoids numerical instability for asymmetric inertia
tensors. The ordering is: half-kick angmom → derive omega → Richardson
quaternion update → re-derive omega. In `integrate_body_final` there is **no**
quaternion update (`body.rs:593`) — only a half-kick of angmom and a final
omega derivation.

### 4. Two quaternions — body→space and body→principal (`lib.rs:58–72`)

`body.quaternion` is the live orientation (body→space), updated each step.
`body.principal_axes` is a fixed rotation (body→principal frame), set at
creation via `diagonalize_inertia`. Wherever a principal-frame operation is
needed, they are composed as `q_ps = quat_mul(body.quaternion, body.principal_axes)`.
Confusing these two is a common error: passing `body.quaternion` alone to
`angmom_to_omega` yields wrong omega if `principal_axes` is not identity.

### 5. Inertia path auto-selection and Monte Carlo noise (`body.rs:400`, `lib.rs:1110`)

`has_overlap` checks all pairs; if any overlap exists, the Monte Carlo path
fires with **exactly 100 000 samples** (`lib.rs:1111`). This is hardcoded.
At 100k samples, the test `montecarlo_single_sphere` (`body.rs:777`) uses
500k samples and still only asserts 5% tolerance. Production runs with the
MC path (any overlapping geometry, including the `sphere7` benchmark) carry
~5% stochastic noise in the principal moments. This is per-body and seeded
fresh each insertion (uses `rand::rng()`), so there is no reproducibility
unless a seed is fixed.

### 6. Richardson quaternion extrapolation (`body.rs:237`)

The orientation update uses `q ← 2*q_half - q_full` with renormalization
(`body.rs:278`). This is second-order accurate but costs ~2× the work of a
naive Euler step. The `dtq` passed to `richardson` is `half_dt = 0.5*dt`
(`body.rs:571`), matching the LIGGGHTS convention. Docs should note that the
Richardson routine is `pub(super)` (private to the `body` module) — users
cannot invoke it directly.

### 7. MPI exchange ordering (`lib.rs:350–427`)

The pre-exchange system `snap_subspheres_to_body_com` (`lib.rs:467`) moves
all sub-sphere atom positions to the body COM before atom exchange, ensuring
all sub-spheres migrate to the same rank as their body. `exchange_bodies`
(`lib.rs:613`) then migrates bodies whose COM left the local subdomain.
`restore_subsphere_positions` (`lib.rs:496`) undoes the snap and rebuilds
actual offsets. If the ordering is violated — e.g., a force plugin writes
sub-sphere positions after snap but before restore — positions will be wrong.
The labels enforce ordering: `snap_subspheres_to_body_com.before("exchange")`,
`exchange_bodies.before("exchange")`, `restore_subsphere_positions.after("exchange")`.

### 8. Ghost cutoff extension (`lib.rs:437`)

`extend_ghost_cutoff_for_clumps` adds `2 * max_R_bound` to `domain.ghost_cutoff`
at setup (`lib.rs:451`). This ensures sub-spheres of bodies near subdomain
boundaries are visible as ghosts on the body-owning rank. Forgetting this
(e.g., adding a clump definition after plugin build) will silently lose contacts.

### 9. Insertion only on rank 0 (`lib.rs:944`)

The `clump_insert_atoms` setup system returns early for ranks != 0. Bodies
and atoms are then distributed via the normal exchange path at step 0.
This means any programmatic clump insertion in user code (calling `insert_clump`
directly) must follow the same rank-0-only discipline or all ranks will create
duplicate atoms.

### 10. Lost-atom diagnostic (`lib.rs:880`)

`check_lost_clump_atoms` runs every 1000 steps and warns (stderr) if a body's
expected sub-sphere count (`body.sub_sphere_tags.len()`) does not match the
count of local atoms with that `body_id`. It does **not** delete atoms or
abort — it is diagnostic only. This can silently mask MPI exchange bugs.

---

## Tutorial outline

A tutorial for `physics/clumps.md` could proceed:

1. **Concept** — what a clump is (rigid multisphere), when to use it vs.
   bonds (permanent rigidity vs. flexible/breakable), and the "no phantom
   parent atom" design.
2. **Minimal config** — dimer definition + one insert block, pointing at an
   existing material. Annotate every field with type and unit.
3. **Adding the plugin** — `app.add_plugin(ClumpPlugin)` (depends on
   `DemAtomPlugin`; order matters).
4. **Contact-exclusion contract** — why `same_body` exists, how to call it,
   what happens if you forget. Include the code snippet from the module doc.
5. **Inertia** — which path fires (non-overlapping vs. overlapping), what
   the 5% noise means practically, and when to care.
6. **Integration mechanics** — one paragraph on angular momentum as state,
   two quaternions, Richardson extrapolation. Readers who extend the integrator
   need this; most readers can skim.
7. **MPI** — ghost cutoff extension (logged at startup), snap/exchange/restore
   ordering, body exchange. Mention that insertion is rank-0-only.
8. **Validation** — `bench_clump_haff_cooling` (7-sphere "sphere7" clump,
   500 bodies, Haff cooling). Config at
   `examples/bench_clump_haff_cooling/config.toml`.

---

## Doc gaps

1. **No documented API for programmatic insertion.** `insert_clump` is `pub`
   (`lib.rs:1097`) but undocumented in the book. Users who want to insert
   clumps from code (not TOML) have no guidance. Should add at least a
   code example showing the call signature and resource access pattern.

2. **Monte Carlo noise is not quantified for users.** The 5% figure is buried
   in inline comments and test tolerances (`body.rs:783`). The config reference
   should warn that overlapping-sphere clumps have stochastic inertia and
   document the hardcoded sample count (100 000).

3. **`compute_clump_inertia` legacy warning is missing from the book.**
   The README and `lib.rs` doc warn about this, but `physics/clumps.md` only
   has a blockquote note that could be missed. The reference/config page has no
   mention.

4. **No documentation on output / recording clump state.** Sub-sphere atoms
   appear in VTP output as ordinary atoms; there is no body-level output
   (COM position, orientation, angmom). Users wanting to track rigid-body
   kinematics have no documented path.

5. **Triclinic/LEBC PBC for body COM is implemented (`lib.rs:559`) but
   undocumented.** The `pbc_multisphere_bodies` function handles LEBC
   streaming-velocity remapping when a body crosses a y-boundary, but this
   is not mentioned anywhere in the docs.

6. **`region` field semantics in `[[clump.insert]]` are undocumented.** The
   config reference (`docs/src/reference/config.md`) should note that `region`
   defaults to the full domain inset by `eff_radius` and that overlap checks
   use a 5% margin.

7. **No guidance on `neighbor.bin_size` for clumps.** The `sphere7` benchmark
   config uses `bin_size = 0.004` with `eff_radius ≈ 0.0011 m`, commenting
   "~4x sub-sphere diameter". Users should know that bin size should be
   significantly larger than the clump bounding radius to avoid neighbor-list
   misses when sub-spheres are spread across bins.

---

## Suggested placement

The existing `docs/src/physics/clumps.md` already covers the architecture
well. The main work is:

- **`physics/clumps.md`** — expand with: programmatic insertion API, Monte
  Carlo noise warning (quantified), LEBC/triclinic PBC note, bin-size
  guidance, and a worked config annotated field-by-field.
- **`reference/config.md`** — add `[clump]` / `[[clump.definitions]]` /
  `[[clump.insert]]` schema table (types, units, defaults, constraints).
- **`reference/validation.md`** — confirm `bench_clump_haff_cooling` entry
  explains what the benchmark tests (rotational energy equipartition and Haff
  decay, not just translational).

The `physics/clumps.md` page is already in `SUMMARY.md` at the correct
position (`docs/src/SUMMARY.md:21`).
