# Planning: `dirt_measure_plane` documentation

**Crate:** `dirt_measure_plane`  
**Source:** `crates/dirt_measure_plane/src/lib.rs`  
**README:** `crates/dirt_measure_plane/README.md`  
**Cargo.toml:** `crates/dirt_measure_plane/Cargo.toml`

---

## Purpose

`dirt_measure_plane` is an opt-in ECS plugin that places one or more infinite
measurement planes into a simulation and counts how many particles cross each plane
in the outward-normal direction each timestep. It is the primary tool for measuring
particle throughput (hopper discharge rate, chute flux, conveyor yield) without
requiring a custom recorder. Results are written to thermo columns only; there is no
programmatic read API.

---

## Public surface to document

| Item | Kind | Notes |
|---|---|---|
| `MeasurePlanePlugin` | `Plugin` | Add to `App`; registers all systems. Entry point for users. |
| `MeasurePlaneDef` | `struct` (public, `Deserialize`) | TOML config for one `[[measure_plane]]` block. All fields public. |
| `MeasurePlanes` | `struct` (public resource) | Opaque — `planes` field is private; no accessor. Document that output is thermo-only. |

Systems registered (internal, not public but worth documenting for ordering):
- `measure_plane_detect_crossings` — every timestep, `PostFinalIntegration`
- `measure_plane_report` — every timestep, `PostFinalIntegration` (no-ops unless `step % report_interval == 0`)

Both run at `ParticleSimScheduleSet::PostFinalIntegration`
(`lib.rs:245-253`).

---

## Config / TOML schema

```toml
[[measure_plane]]
name = "outlet"           # String — unique; used as suffix in thermo keys
point = [0.1, 0.0, 0.0]  # [f64; 3] — any point on the plane (length units)
normal = [1.0, 0.0, 0.0] # [f64; 3] — outward normal direction (auto-normalized)
report_interval = 1000   # usize — averaging window in timesteps; default 1000
```

Schema is `#[serde(deny_unknown_fields)]` (`lib.rs:107`), so typos in field names
are caught at startup. Multiple `[[measure_plane]]` blocks are supported; each is
independent.

**Field semantics:**

| Field | Type | Default | Meaning |
|---|---|---|---|
| `name` | `String` | required | Suffix for thermo keys: `crossings_<name>`, `flow_rate_<name>`, `cross_rate_<name>` |
| `point` | `[f64; 3]` | required | One point on the plane in simulation length units |
| `normal` | `[f64; 3]` | required | Direction toward which crossings are counted; need not be unit-length |
| `report_interval` | `usize` | `1000` | Window length in timesteps; rates are time-averaged over this window then reset |

**Thermo output keys** (for plane named `<name>`):

| Key | Meaning |
|---|---|
| `crossings_<name>` | Cumulative count of positive-direction crossings since simulation start; never reset; global (all-reduced) |
| `flow_rate_<name>` | Mass crossing rate (mass/time) averaged over the last `report_interval` steps |
| `cross_rate_<name>` | Particle crossing rate (count/time) averaged over the last `report_interval` steps |

All three are written at every `report_interval` step, only on rank 0 for the
`println!` side-channel; thermo itself is available on all ranks (`lib.rs:342-344`).

---

## Key behaviors, invariants, and gotchas

### Crossing detection algorithm (`lib.rs:265-291`)

Per-timestep, for each local particle:

1. Compute signed distance `d = (pos - point) · normal` (`lib.rs:187-192`).
2. If the stored `prev_dist ≤ 0` and the new `d > 0`: crossing detected; increment
   `crossings_window`, `total_crossings`, and `mass_window`.
3. Update `prev_signed_dist[tag] = d` unconditionally.

First-seen particles (no entry in `prev_signed_dist`) are never counted as a
crossing on their debut step — only recorded.

### Normal normalization (`lib.rs:163-168`)

The normal is normalized to unit length in `MeasurePlaneState::new`. A near-zero
magnitude (`< 1e-30`) silently falls back to `[1, 0, 0]` with no warning or panic.
A mis-specified normal produces a silently wrong measurement.

### Gross positive crossings, not net flux (`lib.rs:54-58`, README)

The plugin is a **directional gate**: only `≤ 0 → > 0` transitions are counted.
Reverse crossings (`> 0 → ≤ 0`) are neither counted nor subtracted. A particle
oscillating back and forth accumulates a crossing on every forward pass. Totals are
gross positive crossings. Place planes where the flow is reliably one-way (e.g.,
below a hopper outlet, not at a stagnation zone).

**Normal direction matters.** For a plane with normal pointing down (`[0, 0, -1]`),
a particle falling under gravity crosses in the *positive-normal* direction and is
counted exactly once if it never bounces back. See
`examples/measure_plane_throughput/config.toml:103-107` for the canonical pattern.

### MPI handling (`lib.rs:319-324`)

At each `report_interval`, `crossings_window` and `mass_window` are all-reduced
across MPI ranks via `comm.all_reduce_sum_f64`. `total_crossings` is also
all-reduced for the cumulative count. However:

- Detection (`measure_plane_detect_crossings`) runs over `nlocal` atoms only.
- `prev_signed_dist` is a per-rank `HashMap` keyed by atom tag (`lib.rs:144`).
- **When a particle migrates between subdomains, its previous distance entry does
  not follow it.** The receiving rank has no entry, so the first step on the new
  rank records the distance without checking for a crossing. A crossing that
  straddles a migration step is silently missed or double-counted on the origin rank.
  The all-reduce sums whatever each rank accumulated but does not repair the
  migration gap (`lib.rs:66-69`).

This is acceptable for steady, one-way flows but can introduce bias in highly mobile
or multi-stage systems with frequent domain decomposition updates.

### Memory growth (`lib.rs:59-63`)

`prev_signed_dist` accumulates one entry per distinct atom tag the plane has ever
seen. Entries are never evicted. In runs with rate-based insertion (non-recycled
tags) this is a slow monotonic memory leak proportional to the total number of tags
that passed near the plane over the whole run. Not a concern for fixed-population
simulations.

### Window time accuracy (`lib.rs:328-329`)

`window_time = window_steps × dt` uses the *current* `dt`. If `dt` changes within a
reporting window (e.g., between `[[run]]` stages), the denominator is wrong for that
mixed window and rates are only approximate. Rates in single-stage runs with fixed
`dt` are exact.

### Plugin is a no-op if no planes are configured (`lib.rs:237-240`)

If `config.parse_array::<MeasurePlaneDef>("measure_plane")` returns empty, the
plugin inserts an empty `MeasurePlanes` resource and returns without registering any
systems. The resource is still present (safe to query), but neither detection nor
reporting runs.

### Hopper example uses a custom recorder, not `MeasurePlanePlugin`

`bench_hopper_beverloo` implements its own cumulative-discharge CSV tracker
(`examples/bench_hopper_beverloo/main.rs`) rather than relying on
`MeasurePlanePlugin`. It tracks particles that have crossed a z-threshold using a
`HashSet` of discharged tags, not the sign-change algorithm. This is intentional:
the Beverloo benchmark needs a cumulative discharge curve written to CSV for
curve-fitting in Python, which the thermo-only output of `MeasurePlanePlugin` does
not provide. The standalone example (`examples/measure_plane_throughput/`) is the
canonical `MeasurePlanePlugin` usage.

---

## Tutorial outline: measuring hopper discharge rate

Goal: add a measurement plane below a hopper orifice and read `flow_rate_outlet`
from thermo to validate the Beverloo scaling in a simplified simulation.

1. **Add the plugin** — `app.add_plugins(MeasurePlanePlugin)` after
   `GranularDefaultPlugins` (no ordering constraint on `MeasurePlanePlugin`; it only
   depends on `Atom`, `CommResource`, `RunState`, and `Thermo`).
2. **Configure the plane** — add one `[[measure_plane]]` block to the TOML:
   ```toml
   [[measure_plane]]
   name = "outlet"
   point = [0.0, 0.0, 0.03]    # z just below orifice
   normal = [0.0, 0.0, -1.0]   # downward: falling particles cross positive-normal
   report_interval = 5000       # ~0.1 s at dt = 2e-5
   ```
   See `examples/measure_plane_throughput/config.toml:101-107` for a working example.
3. **Run** and observe thermo columns `crossings_outlet`, `flow_rate_outlet`,
   `cross_rate_outlet` in the CSV/stdout output.
4. **Read in Python** — parse the thermo CSV and take the mean of `flow_rate_outlet`
   over the steady discharge window to get W (mass/time). Compare to the Beverloo
   prediction `W = C ρ_b √g (D - kd)^(3/2)`.
5. **Multiple planes** — stack a second plane higher in the funnel to check that W
   is the same at both planes (conservation of mass in steady flow is a good
   sanity check).

**Pitfall to highlight:** if the normal points upward (`[0, 0, 1]`) and the orifice
is at the bottom, falling particles cross in the negative-normal direction and are
never counted. Always orient the normal *with* the flow direction.

---

## Doc gaps

1. **No example that combines `MeasurePlanePlugin` with the Beverloo benchmark.**
   `bench_hopper_beverloo` uses a custom recorder; the `measure_plane_throughput`
   example does not do a Beverloo sweep. A note in the diagnostics page linking
   the two would clarify the relationship.
2. **`MeasurePlanes` resource opaqueness** is documented in crate docs and in
   `docs/src/physics/diagnostics.md` but the *reason* (no programmatic read API
   is intentional so all I/O goes through the thermo pipeline) is not explained.
   Worth one sentence.
3. **Interaction with multi-stage `[[run]]` blocks** is not documented. The
   `report_interval` counts total timesteps (`run_state.total_cycle`), so a window
   that straddles a stage boundary mixes `dt` values and reports a blended rate.
   This is the variable-`dt` caveat but the stage-boundary case is more concrete.
4. **No guidance on `report_interval` sizing.** For Beverloo applications, the
   interval should span several grain transit times; a sentence on choosing the
   interval relative to `dt` and expected crossing timescale would help.
5. **Silent fallback for degenerate normals** — the `< 1e-30` threshold and the
   `[1, 0, 0]` fallback are noted in caveats but a user-facing warning or panic
   would be safer. Flag as a potential improvement.
6. **Memory leak not bounded.** Document the expected steady-state map size for
   fixed-population runs (bounded by N_particles) vs. continuous-insertion runs
   (unbounded). The distinction is not currently drawn.

---

## Suggested placement

The content belongs in **`docs/src/physics/diagnostics.md`**, which already has a
complete measurement-plane section (`diagnostics.md:9-59`) and a contact-analysis
section. The existing diagnostics page covers the caveats well; the gaps above
(multi-stage `dt`, choosing `report_interval`, Beverloo tutorial link) should be
added there as short subsections or callout boxes.

A cross-reference from the Beverloo example README (`examples/bench_hopper_beverloo/README.md`)
to the diagnostics page (explaining why the example uses a custom recorder instead of
`MeasurePlanePlugin`) would close the conceptual loop for users who arrive via the
example.

No new top-level chapter is warranted; `MeasurePlanePlugin` is a diagnostic
instrument, not a physics model, and the diagnostics chapter is the right home.
