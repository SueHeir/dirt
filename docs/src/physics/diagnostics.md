# Diagnostics: Measurement Planes & Contact Analysis

Two opt-in plugins instrument a running simulation: **measurement planes** count
particles crossing a surface (throughput), and **contact analysis** characterizes
the packing (coordination number, rattlers, fabric tensor, per-contact CSV).
Both publish their results as thermo keys.

## Measurement planes

A measurement plane is an infinite plane defined by a point and a normal. Each
step the plugin computes the signed distance of every local particle from the
plane; when a particle's signed distance changes from `≤ 0` to `> 0` between
consecutive steps it is counted as having **crossed** in the positive-normal
direction. Add `MeasurePlanePlugin` and configure `[[measure_plane]]` blocks.

```toml
[[measure_plane]]
name = "outlet"           # unique name; used in thermo output keys
point = [0.1, 0.0, 0.0]   # any point on the plane
normal = [1.0, 0.0, 0.0]  # outward normal (automatically normalized)
report_interval = 1000    # averaging window in timesteps (default 1000)
```

### Output: thermo keys only

All results are exposed **only as thermo keys**, written every `report_interval`
steps. There is no public read API (the `MeasurePlanes` resource is opaque), so
downstream code reads the thermo columns. For a plane named `<name>`:

- `crossings_<name>` — cumulative crossing count (positive direction, global
  all-reduced, never reset).
- `flow_rate_<name>` — mass flow rate, averaged over the window.
- `cross_rate_<name>` — particle crossing rate (1/time), averaged over the
  window.

### Caveats — it is a directional gate, not a flux meter

Read these before trusting the numbers:

- **Directional, not net flux.** Only `≤ 0 → > 0` transitions are counted.
  Reverse crossings are ignored entirely — neither counted nor subtracted. A
  particle that oscillates across the plane is **recounted** on every forward
  pass, so the totals are *gross* positive crossings, not net throughput. Place
  planes where flow is essentially one-way (e.g. below a hopper outlet).
- **`prev_signed_dist` grows without bound.** The per-plane state stores one map
  entry per atom tag it has *ever* seen and never evicts them. In a long run
  with continuous insertion this is a slow memory leak proportional to the
  number of distinct tags seen near the plane.
- **MPI migration can mis/double-count.** Crossing detection runs over `nlocal`
  only, and the previous-distance map lives independently on each rank. When a
  particle migrates between subdomains its previous distance does not follow it,
  so a crossing straddling a migration step can be missed or counted on the
  wrong rank. Summing across ranks at report time does not repair this.
- **Variable `dt` makes the window time approximate.** `window_time` uses the
  *current* timestep; if `dt` changes within a window (e.g. across run stages)
  the reported rates are only approximate for that window.
- **Degenerate normal silently falls back to `[1, 0, 0]`.** A normal with
  magnitude `< 1e-30` is replaced by the +x direction without warning — a
  mis-specified plane silently measures the wrong cross-section.

## Contact analysis

`ContactAnalysisPlugin` characterizes a settled packing. Everything is computed
in a single pass over the neighbor list and toggled from `[contact_analysis]`:

```toml
[contact_analysis]
interval = 1000         # dump per-contact CSV every N steps (0 = disabled)
coordination = true     # per-atom coordination number
rattlers = true         # detect rattlers (< 4 contacts)
fabric_tensor = true    # fabric tensor to thermo
file_prefix = "contact" # prefix for contact CSV filenames
```

- **Coordination number** — active contacts per particle, exposed as thermo
  (`coord_avg`, `coord_max`, `coord_min`) and as a per-atom dump scalar
  (`coordination`).
- **Rattler detection** — particles with fewer than 4 contacts are mechanically
  unstable in 3D (they lack the *d* + 1 = 4 constraints for static equilibrium).
  Reports `n_rattlers` and `rattler_fraction`.
- **Fabric tensor** — `F_ij = (1/Nc) Σ nᵢ nⱼ` measures the directional
  distribution of contact normals. Its trace is always 1 (unit normals); an
  isotropic packing gives `F ≈ (1/3) I`, and shear-induced anisotropy shows in
  the diagonal. Six components go to thermo (`fabric_xx`, …, `fabric_yz`) with a
  `contacts` count.
- **Per-contact CSV** — geometry only (tags, overlap, contact point, contact
  normal) — **no force**. Force data lives in the contact plugin's
  tangential-history store and is not coupled here.

### Newton's-third-law accounting

Every metric is computed in one neighbor-list pass, so each depends on whether
the list uses **Newton's third law** — whether a pair `(i, j)` is visited once
(half list, `newton == true`) or twice (full list, `newton == false`). The three
accumulators correct for this so the reported quantities are list-independent:

- **Coordination.** When `newton`, the single visit increments both endpoints
  (`j` only when local). When `!newton`, each visit increments only `i`; across
  the two visits each particle ends with the right count.
- **Fabric tensor.** Each contact's `n ⊗ n` must contribute once: weight `1.0`
  under `newton`, weight `0.5` per visit under `!newton`. The contact count is
  incremented by the same weight, so the normalized tensor is identical either
  way.
- **Per-contact CSV.** Each physical contact appears once: the single visit
  under `newton`; the `i < j` ordering only under `!newton`.

### Plugin-ordering contract

`ContactAnalysisPlugin` does not own a force or output pipeline — it hooks into
existing ones, so **registration order matters**:

- A **Hertz–Mindlin contact plugin must be registered first.** The analysis runs
  `.after("hertz_mindlin_contact")` in `PostForce`; with no system carrying that
  label the scheduler has no ordering anchor.
- **`PrintPlugin` must be registered first** when `coordination = true` — the
  plugin registers the `coordination` dump scalar against the `DumpRegistry` at
  build time, and **panics** with `"DumpRegistry not found — PrintPlugin must be
  added first"` if it is absent.

In practice both are satisfied by adding `GranularDefaultPlugins` (contact) and
`CorePlugins` (which includes `PrintPlugin`) **before** `ContactAnalysisPlugin`:

```rust
app.add_plugins(CorePlugins)            // includes PrintPlugin
   .add_plugins(GranularDefaultPlugins) // labels "hertz_mindlin_contact"
   .add_plugins(GravityPlugin)
   .add_plugins(WallPlugin)
   .add_plugins(ContactAnalysisPlugin); // must come after the two above
```
