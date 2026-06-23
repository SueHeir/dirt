# Documentation Plan: `dirt_contact_analysis`

Sources read: `crates/dirt_contact_analysis/src/lib.rs`, `Cargo.toml`, `README.md`,
`docs/src/physics/diagnostics.md`, `docs/src/reference/config.md`.

---

## Purpose

`dirt_contact_analysis` characterizes the contact network of a DEM packing in a
single neighbor-list traversal each step. It is an **opt-in diagnostic** — it
contributes no forces and imposes no model coupling — but it does depend on the
neighbor list and the `DemAtom` radius store to detect contacts (overlap > 0).
The four analyses it offers are independent flags in a single `[contact_analysis]`
TOML section: coordination number, rattler detection, fabric tensor, and
per-contact CSV dump.

---

## Public Surface to Document

| Symbol | Kind | Description |
|--------|------|-------------|
| `ContactAnalysisPlugin` | `Plugin` | Entry point; reads config, registers resources and systems. `lib.rs:275` |
| `ContactAnalysisConfig` | `struct` (public) | Serde-deserialized TOML config; all five fields are `pub`. `lib.rs:127` |
| `ContactAnalysis` | `AtomData` struct | Per-atom `coordination: Vec<f64>`, tagged `#[zero]` (reset each force pass) and `#[forward]` (ghost comm). `lib.rs:165` |
| `ContactRecord` | `struct` (public) | One contact's geometry: `i_tag`, `j_tag`, `overlap`, `cx/cy/cz`, `nx/ny/nz`. `lib.rs:197` |
| `ContactOutput` | Resource (public) | `records: Vec<ContactRecord>` for the current step; cleared at the start of each analysis pass. `lib.rs:220` |
| `FabricTensorAccum` | Resource (private) | Running sum for the six tensor components plus `nc`; private but documented. `lib.rs:246` |

Systems (all private, but their scheduling is a public contract):

| System | Schedule set | Ordering |
|--------|-------------|---------|
| `compute_contact_analysis` | `PostForce` | `.label("contact_analysis").after("hertz_mindlin_contact")` — `lib.rs:328` |
| `push_coordination_to_thermo` | `PostForce` | `.after("contact_analysis")` — `lib.rs:319` |
| `push_fabric_tensor_to_thermo` | `PostForce` | `.after("contact_analysis")` — `lib.rs:343` |
| `dump_contact_records` | `PostFinalIntegration` | (unconditional ordering, added only when `interval > 0`) — `lib.rs:336` |

---

## Config / TOML Schema

Section key: `[contact_analysis]`.  Uses `#[serde(deny_unknown_fields)]`, so
any unknown key is a hard error (`lib.rs:126`).

| Key | Rust type | Default | Meaning |
|-----|-----------|---------|---------|
| `interval` | `usize` | `0` | Dump per-contact CSV every N steps. `0` disables CSV output entirely. The gate check is `step % interval == 0`. `lib.rs:133` |
| `coordination` | `bool` | `false` | Compute per-atom contact count; pushes `coord_avg/max/min` to thermo and registers `coordination` as a per-atom dump scalar. `lib.rs:139` |
| `rattlers` | `bool` | `false` | Detect particles with < 4 contacts (requires `coordination = true` to have per-atom data). Pushes `n_rattlers` and `rattler_fraction` to thermo. `lib.rs:145` |
| `fabric_tensor` | `bool` | `false` | Compute F_ij = (1/Nc) Σ n_i n_j; pushes six tensor components and `contacts` to thermo. `lib.rs:151` |
| `file_prefix` | `String` | `"contact"` | Prefix for CSV filenames. Default provided by `default_file_prefix()` via `#[serde(default)]`, NOT by `Default::default()` which gives `""`. See doc-gap below. `lib.rs:117,153` |

Default-config snippet emitted by `Plugin::default_config` (`lib.rs:279`):

```toml
[contact_analysis]
interval = 0
coordination = true
rattlers = false
fabric_tensor = false
file_prefix = "contact"
```

Note: the plugin default has `coordination = true`; the struct field default is
`false`. The snippet wins when a user generates config via the framework hook,
but a hand-written TOML with no `coordination` key gets `false`.

---

## Key Behaviors, Invariants & Gotchas

### Contact detection threshold (`lib.rs:408–419`)

A pair is "in contact" when `delta = (r_i + r_j) - distance > 0`.  The guard
`distance == 0.0` prevents a divide-by-zero from identical-position atoms; those
pairs are silently skipped, which may mask initialization bugs.

### Newton's-third-law accounting (`lib.rs:62–79`, `lib.rs:424–455`)

All three accumulators are corrected for half vs. full neighbor lists so reported
quantities are list-independent:

- **Coordination** (`lib.rs:425–429`): Newton on → increment both `i` and `j`
  (but `j` only if `j < nlocal`, because ghost atoms are owned and counted by
  another rank). Newton off → increment only `i`; the pair is visited again as
  `(j, i)` where `j` becomes the new `i`.
- **Fabric tensor** (`lib.rs:441–451`): Newton on → full weight `vs = 1.0`;
  Newton off → half weight `vs = 0.5` per visit, so the sum and `nc` are
  identical either way.
- **Per-contact CSV** (`lib.rs:455`): Newton on → always record; Newton off →
  record only when `i < j` to suppress the duplicate visit.

### MPI reduction strategy (`lib.rs:520–526`)

Coordination stats use `all_reduce_sum_f64` (sum) for `coord_sum` and
`n_rattlers`, then divide by `atoms.natoms` (global atom count). For
`coord_max`, there is no `all_reduce_max` available, so it is computed as
`-all_reduce_min(-max_val)` (`lib.rs:523`). Fabric tensor components are each
individually all-reduced with sum (`lib.rs:571–577`), then normalized by global
`nc` after reduction.

### Thermo gating (`lib.rs:489–491`, `lib.rs:563–565`)

Coordination and fabric tensor pushes are guarded by `thermo.interval == 0 ||
step % thermo.interval != 0` — they fire only on thermo steps. CSV dumps are
guarded by `step % config.interval == 0` in both `compute_contact_analysis`
(record collection, `lib.rs:370`) and `dump_contact_records` (file write,
`lib.rs:611`), so records are only collected on dump steps.

### `ContactAnalysis` reset (`lib.rs:163`)

The `#[zero]` attribute resets `coordination` to zero at the start of each force
computation. This means per-atom coordination reflects contacts at the **current
step only** — there is no cumulative counter.

### Plugin-ordering panics (`lib.rs:306–311`)

When `coordination = true`, `build()` fetches `DumpRegistry` and will
**panic** with `"DumpRegistry not found — PrintPlugin must be added first"` if
`PrintPlugin` (part of `CorePlugins`) was not registered. The fix is to add
`CorePlugins` before `ContactAnalysisPlugin`. The `hertz_mindlin_contact` label
dependency is soft (scheduler ordering only) — missing it will not panic but
will produce incorrect results if the analysis runs before forces settle.

### CSV file naming convention (`lib.rs:638`)

Files are written to `<output_dir>/contact/<prefix>_<NNNNNN>_rank<rank>.csv`
where `NNNNNN` is the step number zero-padded to 6 digits. The subdirectory
`contact/` is created automatically (`lib.rs:637`). CSV errors print a warning
and do not abort the run (`lib.rs:622`).

### `rattlers` requires `coordination` (`lib.rs:536`)

`config.rattlers` is checked inside the `push_coordination_to_thermo` system,
which only runs when `config.coordination = true` (the system is only added in
that branch, `lib.rs:319`). Setting `rattlers = true` with `coordination =
false` silently produces no rattler output.

### CSV records contain geometry only, no force (`lib.rs:3–6`)

`ContactRecord` holds `overlap`, contact point, and contact normal. Tangential
force / history is stored in the Hertz-Mindlin plugin and is not accessible
here. Document explicitly so users do not expect force columns.

---

## Tutorial Outline

Suggested flow for the tutorial section in `physics/diagnostics.md`:

1. **"What does it measure?"** — one paragraph: coordination Z, rattlers,
   fabric tensor, CSV; distinguish geometry-only from force data.
2. **Minimal config** — enable coordination + rattlers only; show thermo output
   and what the numbers mean physically (Z ≈ 6 for random dense packing in 3D,
   rattler fraction as a packing quality check).
3. **Fabric tensor** — add `fabric_tensor = true`; explain F trace = 1,
   isotropic → diagonal ≈ 1/3, shear anisotropy; show how to plot eigenvalues
   from thermo columns.
4. **Per-contact CSV** — add `interval = 1000`; explain filename pattern; show a
   Python snippet loading a single file (`pd.read_csv`) and plotting contact-
   normal rose diagram.
5. **Plugin ordering** — brief callout box: add `CorePlugins` and
   `GranularDefaultPlugins` before `ContactAnalysisPlugin`; consequence of
   getting it wrong (panic vs. silent mis-order).
6. **MPI note** — CSV is rank-local; to reconstruct the global contact network,
   concatenate all `rank*` files for the same step.

---

## Doc Gaps

1. **`file_prefix` serde vs. `Default` mismatch** (`lib.rs:117,153`): The
   `Default` impl gives `file_prefix = ""` but the serde default function gives
   `"contact"`. A user instantiating `ContactAnalysisConfig::default()` in Rust
   code (rather than deserializing TOML) gets the empty string. The in-code doc
   comment and module-level docs say the default is `"contact"`, which is only
   true after TOML deserialization. Either the test at `lib.rs:919` documents
   the gap or it should be closed by implementing `Default` to call
   `default_file_prefix()`.

2. **`rattlers` implicitly requires `coordination = true`** but this is not
   validated at build time — no error or warning is emitted if `rattlers = true,
   coordination = false`. A validation step in `build()` with a clear error
   message would improve UX.

3. **No `interval > 0` guard in `default_config`** — the plugin default emits
   `interval = 0` (CSV disabled) but `coordination = true`. A new user copying
   the default gets coordination but no CSV. Good default, but could use a
   comment in the emitted snippet explaining the 0 = disabled semantic.

4. **`ContactOutput` is public** (`lib.rs:220`) but there is no documented use
   case for reading it downstream. If the intent is to allow other plugins to
   consume the contact network (e.g., for thermal conduction), that should be
   documented; if not, it could be private.

5. **No per-contact force** — the README and module doc both state this clearly,
   but the docs do not explain how a user *would* access tangential force if they
   needed it (pointer to Hertz-Mindlin tangential-history store would close this).

6. **`FabricTensorAccum` is private** but its fields are individually documented.
   The normalization-after-reduction pattern (accumulate locally, MPI reduce,
   then divide by global `nc`) is a subtle correctness point that deserves a
   callout in the user-facing docs, not just inline comments.

7. **No example in `examples/`** — there is no benchmark or tutorial example
   demonstrating `ContactAnalysisPlugin` in a real simulation. A minimal
   `examples/bench_contact_analysis/` or integration into an existing example
   (e.g., angle-of-repose at rest) would make the tutorial concrete.

---

## Suggested Placement

**Primary location**: `docs/src/physics/diagnostics.md` — the chapter already
exists and already has a `## Contact analysis` section with the correct physics
descriptions. That section should be expanded with the tutorial outline above,
the Newton accounting note (already present in abbreviated form), and the plugin
ordering callout. The doc gaps (especially `rattlers` requiring `coordination`,
and `file_prefix` default mismatch) belong in the caveats subsection there.

**Secondary location**: `docs/src/reference/config.md` — the `[contact_analysis]`
entry (`config.md:202`) is currently just a two-line snippet with a pointer to
`diagnostics.md`. Expand it to a full field table (matching the schema section
above) for reference completeness.

**No new chapter needed.** The crate is a single-plugin diagnostic that fits
cleanly under the existing diagnostics chapter structure.
