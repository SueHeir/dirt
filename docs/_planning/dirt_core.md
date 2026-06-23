# Planning: `dirt_core` documentation

Source audited: `crates/dirt_core/src/lib.rs`, `crates/dirt_core/Cargo.toml`,
`crates/dirt_core/README.md`, and the full `docs/src/` tree.

---

## Purpose

`dirt_core` is the batteries-included umbrella crate. It has no physics logic of
its own. Its job is three things:

1. **Re-export** every sub-crate in the GRASS → SOIL → DIRT stack so users depend
   on exactly one crate.
2. **Define `CorePlugins`** — the plugin group for simulation infrastructure
   (input/config, communication, domain decomposition, neighbor lists, groups, run
   loop, output).
3. **Define `prelude`** — a flat star-import that surfaces every type needed for a
   typical simulation, including types from `grass_*`, `soil_*`, and `dirt_*`.

`GranularDefaultPlugins` (the DEM-physics half) lives in `dirt_granular`, but is
re-exported through the prelude so users never import it from two places.

---

## Public surface to document

### Prelude contents (`crates/dirt_core/src/lib.rs:246–278`)

| Category | Exported names |
|---|---|
| Plugin groups | `CorePlugins`, `GranularDefaultPlugins` |
| DEM atom + insertion | `DemAtomPlugin`, `DemConfig`, `MaterialTable`, `DemAtomInsertPlugin`, `ParticlesConfig` |
| Bonds | `DemBondPlugin` |
| Clumps | `ClumpPlugin`, `ClumpRegistry`, `ClumpAtom`, `ClumpDef`, `MultisphereBody`, `MultisphereBodyStore` |
| Granular contact | `HertzMindlinContactPlugin`, `RotationalDynamicsPlugin`, `GranularTempPlugin` |
| Contact analysis | `ContactAnalysisPlugin`, `ContactAnalysisConfig` |
| Measurement planes | `MeasurePlanePlugin`, `MeasurePlanes`, `MeasurePlaneDef` |
| Walls | `WallPlugin`, `Walls`, `WallDef`, `WallPlane`, `WallMotion` |
| Box deformation | `DeformPlugin`, `DeformConfig`, `DeformState` |
| Fixes (DEM) | `FixesPlugin`, `GravityPlugin`, `GravityConfig`, `AddForceDef`, `SetForceDef`, `FreezeDef`, `MoveLinearDef`, `ViscousDef`, `FixesRegistry` |
| Fixes (substrate) | `SoilFixesPlugin`, `PinDef`, `PinRegistry`, `PinState` |
| Derive macros | `ScheduleSet`, `StageEnum`, `AtomData` |
| Framework (glob) | `App`, `Plugin`, `PluginGroup`, `Res`, `ResMut`, `ScheduleSet` (trait) via `grass_app::prelude::*` + `grass_scheduler::prelude::*` |
| Core types (glob) | `Atom`, `Config`, `RunState`, `ParticleSimScheduleSet`, all of `soil_core::*` |
| Print (glob) | all of `soil_print::*` |
| Verlet (glob) | all of `soil_verlet::*` |

Note: `ParticleSimScheduleSet` is re-exported explicitly (line 275) to prevent
name collision with the `ScheduleSet` trait that arrives via the glob.

### `CorePlugins` membership (`lib.rs:195–211`)

Registration order matters for dependency injection — plugins registered earlier
can insert resources that later plugins consume.

| Order | Plugin | Source crate | Role |
|---|---|---|---|
| 0 | warning-fn setup (closure) | `soil_core::verlet_schedule_warnings` | Installs schedule-ordering warnings for common mistakes |
| 1 | `InputPlugin` | `soil_core` | CLI parsing, banner, TOML config load (skipped if `Config` already present) |
| 2 | `CommunicationPlugin` | `soil_core` | MPI or single-process comm backend (selected by `mpi_backend` feature) |
| 3 | `DomainPlugin` | `soil_core` | Domain decomposition, PBC, shrink-wrap |
| 4 | `NeighborPlugin` | `soil_core` | Bin-based neighbor list construction |
| 5 | `GroupPlugin` | `soil_core` | Atom group definitions and filtering |
| 6 | `RunPlugin` | `soil_core` | Run/cycle loop management |
| 7 | `PrintPlugin` | `soil_print` | Thermo output, dump files (CSV/binary/VTP), restart I/O |

**Not included:** `VelocityVerletPlugin` (translational integration),
`RotationalDynamicsPlugin`, any force/contact plugin.

### `GranularDefaultPlugins` membership (`crates/dirt_granular/src/lib.rs:179–`)

Defined in `dirt_granular`, re-exported via prelude. Registration order:

1. `DemAtomPlugin` — per-atom material (radius, density, `MaterialTable`)
2. `DemAtomInsertPlugin` — random/rate/file insertion from `[[particles.insert]]`
3. `VelocityVerletPlugin` — translational Velocity Verlet (from `soil_verlet`)
4. `HertzMindlinContactPlugin` — Hertz normal + Mindlin tangential contact
5. `RotationalDynamicsPlugin` — quaternion angular integration (I = 2/5 m r²)

**Opt-in, not bundled:** `GranularTempPlugin` (writes `data/GranularTemp.txt`).

### Feature flags (`Cargo.toml:12–14`)

| Feature | Default | Effect |
|---|---|---|
| `mpi_backend` | **yes** | Enables `soil_core/mpi_backend` + `dirt_clump/mpi_backend`; links an MPI library for multi-rank domain-decomposed runs. Disable with `--no-default-features` for a single-process serial binary. |

`soil_core` is depended on with `default-features = false`; MPI is only threaded
through when `dirt_core`'s own `mpi_backend` feature is active.

### Crate-level re-exports

`pub use dirt_atom;`, `pub use dirt_bond;`, `pub use dirt_clump;`, etc. — each
sub-crate is re-exported at the top level so users can do e.g.
`dirt_core::dirt_granular::…` or just `use dirt_core::prelude::*`. See
`lib.rs:120–159` for the full list with doc comments.

---

## Config / TOML schema consumed by `CorePlugins`

These are the sections `CorePlugins` plugins read. They are parsed by the SOIL /
GRASS substrate crates; the schema lives there, but users configure them here.

### `[comm]` (read by `CommunicationPlugin`)
```toml
[comm]
processors_x = 1
processors_y = 1
processors_z = 1      # product must equal MPI rank count
```

### `[domain]` (read by `DomainPlugin`)
```toml
[domain]
x_low = 0.0
x_high = 0.04
y_low = 0.0
y_high = 0.02
z_low = 0.0
z_high = 0.08
boundary_x = "fixed"     # "fixed" | "periodic"
boundary_y = "periodic"
boundary_z = "fixed"
```

### `[neighbor]` (read by `NeighborPlugin`)
```toml
[neighbor]
skin_fraction = 1.1   # search-radius multiplier; 1.0–1.5 typical
bin_size = 0.005      # bin width [m]; must be ≥ largest particle diameter
every = 1             # rebuild interval in steps
```

### `[[run]]` (read by `RunPlugin`)
```toml
[[run]]
name = "settle"
dt = 1.0e-5
steps = 200000
thermo = 2000
```
Multiple `[[run]]` blocks define a staged run.

### `[output]` / `[vtp]` (read by `PrintPlugin`)
```toml
[output]
dir = "examples/my_run"

[vtp]
interval = 2000
```

Note: sections consumed by `GranularDefaultPlugins` (`[dem]`,
`[[dem.materials]]`, `[[particles.insert]]`) are documented separately in the
contact and insertion pages and are NOT read by `CorePlugins`.

---

## Key behaviors, invariants, and gotchas

### `CorePlugins` does not move particles (`lib.rs:58–64`)

The warning comment at `lib.rs:59–64` is the single most important gotcha:
`CorePlugins` alone produces a run that reads config, builds neighbor lists, and
prints output — but never advances particle positions. Velocity Verlet lives in
`GranularDefaultPlugins`. A minimal run that prints thermo but does zero physics
is legal and intentional (useful for testing infrastructure), but easy to
introduce accidentally.

### `GranularTempPlugin` is opt-in (`lib.rs:66–68`)

`GranularTempPlugin` is not in `GranularDefaultPlugins`. Users who want
`data/GranularTemp.txt` must add it explicitly:
```rust
app.add_plugins(GranularDefaultPlugins)
   .add_plugins(GranularTempPlugin);
```

### `mpi_backend` vs `--no-default-features` (`Cargo.toml:12–14`)

`mpi_backend` is a default feature. Single-process builds — including all
single-machine dev work and examples — should pass `--no-default-features` to
avoid requiring an MPI toolchain:
```bash
cargo build --release --no-default-features
```
For multi-rank parallel runs, drop the flag and ensure OpenMPI or MPICH is
installed.

### Git dependencies — no crates.io (`Cargo.toml:17–21`, `README.md:77`)

All `grass_*` and `soil_*` dependencies are git-sourced from
`https://github.com/SueHeir/grass` and `https://github.com/SueHeir/soil`
(branch pinned to `main`/`master`). Because crates.io forbids git deps in
published crates, `dirt_core` cannot be published there. Users must depend on it
via git URL:
```toml
dirt_core = { git = "https://github.com/SueHeir/dirt" }
```

### `soil_core` dep uses `default-features = false` (`Cargo.toml:20`)

`soil_core` is pulled in with `default-features = false`; `dirt_core` threads
MPI through only when its own `mpi_backend` feature is active. This means
feature unification works correctly: enabling `mpi_backend` on `dirt_core`
propagates down to both `soil_core` and `dirt_clump`.

### `CommunicationPlugin` is backend-transparent

`CommunicationPlugin` resolves to `SingleProcessComm` or `MpiCommBackend`
depending on the feature flag — the application code sees the same `CommResource`
either way (as shown in `first-simulation.md:88`, `comm.all_reduce_sum_f64`).

### Warning function (`lib.rs:198–200`)

The first thing `CorePlugins::build` does is install
`soil_core::verlet_schedule_warnings` as the app-level warning function. This
catches common scheduling mistakes (e.g. adding a force system in the wrong
schedule set) and prints them to stderr at startup.

### `InputPlugin` skip-if-present (`lib.rs:167`)

`InputPlugin` is skipped if a `Config` resource is already registered. This
allows tests and integration harnesses to inject a pre-built `Config` and bypass
CLI / file I/O entirely.

---

## Tutorial outline — minimal `main.rs`

The canonical two-plugin pattern (cite: `lib.rs:14–23`, `README.md:31–39`,
`introduction.md:42–49`, `installation.md:65–73`):

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)            // I/O, comm, domain, neighbors, run, print
       .add_plugins(GranularDefaultPlugins); // Hertz–Mindlin contact + Velocity Verlet
    app.start();
}
```

The config file is passed on the command line:
```bash
cargo run --release --no-default-features -- config.toml
```

**Extending the assembly** — opt-in plugins bolted after `GranularDefaultPlugins`:

```rust
app.add_plugins(CorePlugins)
   .add_plugins(GranularDefaultPlugins)
   .add_plugins(GravityPlugin)      // reads [gravity] from TOML
   .add_plugins(WallPlugin)         // reads [[wall]] from TOML
   .add_plugins(DemBondPlugin)      // reads [bonds]
   .add_plugins(ClumpPlugin)        // reads [[clump.definitions]] + [[clump.insert]]
   .add_plugins(GranularTempPlugin) // writes data/GranularTemp.txt
   .add_plugins(ContactAnalysisPlugin);
```

**Finer control** — bypass plugin groups and add individual plugins:

```rust
app.add_plugins(CorePlugins)
   .add_plugins(soil_core::InputPlugin)       // already in CorePlugins; just for illustration
   .add_plugins(VelocityVerletPlugin)         // from soil_verlet, via prelude
   .add_plugins(HertzMindlinContactPlugin);   // from dirt_granular, via prelude
```

---

## Doc gaps

1. **No dedicated `dirt_core` page in the book.** The crate is described
   implicitly in `introduction.md` and the plugin-group table in `physics/contact.md`
   but has no standalone reference page. Users have no single place to see the
   full prelude contents or CorePlugins membership.

2. **Prelude contents not enumerated anywhere in the book.** `introduction.md`
   mentions "you get the whole stack re-exported" but never lists what that means.
   The only complete listing is in the Rust doc comment at `lib.rs:213–278`.

3. **`GranularTempPlugin` opt-in is mentioned in `contact.md` but not in any
   getting-started or config-reference page.** Easy to miss.

4. **`mpi_backend` / `--no-default-features` split is explained only in
   `installation.md`.** Config-reference and physics pages don't flag that serial
   builds need the flag; the hopper example command in `first-simulation.md` uses
   it but doesn't explain why.

5. **The `InputPlugin` skip-if-present behavior** (useful for testing) is
   documented only in the rustdoc, not in the book.

6. **`verlet_schedule_warnings` is not mentioned in user-facing docs.** A user
   who adds systems in the wrong set gets useful stderr output without knowing
   where it came from.

7. **`stack/overview.md` is a stub** (`stack/overview.md:36`): the promised
   "concrete trace of one timestep through all three tiers" is missing.

8. **Clump (`dirt_clump`) prelude exports** — `ClumpAtom`, `MultisphereBody`,
   `MultisphereBodyStore` — appear in the prelude (`lib.rs:254`) but are not
   mentioned in `physics/clumps.md` or anywhere in the book.

---

## Suggested placement

| Doc location | Content |
|---|---|
| **`introduction.md`** (already exists) | Keep the 3-line code snippet; add a one-paragraph expansion of what `CorePlugins` vs `GranularDefaultPlugins` each do; link to the new reference page. |
| **`getting-started/installation.md`** (already exists) | Already explains `--no-default-features`; add a sentence noting this is the `mpi_backend` feature toggle. |
| **`getting-started/first-simulation.md`** (already exists) | Already shows the two-plugin assembly; add the opt-in-plugin pattern (gravity, walls, bonds). |
| **NEW: `reference/dirt_core.md`** | Dedicated reference page: full `CorePlugins` membership table (ordered), full prelude contents table, feature flags, git-dep / crates.io caveat, `GranularTempPlugin` opt-in call-out, `InputPlugin` skip-if-present. |
| **`reference/config.md`** (already exists) | Already has the `CorePlugins` TOML section (`config.md:14`); consider adding a note next to `[neighbor]` about `bin_size ≥ largest diameter`. |
| **`stack/overview.md`** (stub) | Fill the promised timestep trace; this is a natural companion to the `CorePlugins` page once the architecture is written up. |
