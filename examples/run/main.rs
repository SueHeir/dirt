//! run — the generic, config-only DIRT driver.
//!
//! This is a *single* example binary that assembles the standard DIRT plugin
//! stack and then hands everything — geometry, materials, particle insertion,
//! walls, body forces, run stages, and output — to a TOML config supplied on
//! the command line. There is **no per-case Rust**: every case is a `config.toml`.
//!
//! ```bash
//! cargo run --release --example run -- examples/run/config.toml
//! ```
//!
//! Shipped configs (each a complete, self-contained scenario — pick one on argv):
//!
//! | Config | Scenario | Geometry | Insertion |
//! |--------|----------|----------|-----------|
//! | `config.toml`        | two-sphere settle       | floor plane           | 2 fixed spheres |
//! | `pour_settle.toml`   | granular pour & settle  | box (5 plane walls)   | polydisperse cloud |
//! | `pour_cylinder.toml` | bidisperse silo pour    | cylinder wall + floor | two materials |
//!
//! The plugin set below is deliberately a superset of what any one simple case
//! needs. Each plugin is a no-op when its config section is absent (gravity
//! defaults to zero body force only if you set it; walls/fixes register nothing
//! when no `[[wall]]` / `[[*]]` fixes are declared), so the same driver runs a
//! two-particle rebound, a settle-into-a-box, and a granular pour — chosen
//! entirely by the config.
//!
//! Standard dump output is produced via the config's `[dump]` (per-atom CSV /
//! binary snapshots) and/or `[vtp]` sections, handled by the core `PrintPlugin`.

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app
        // Infrastructure: CLI parse + TOML load, comm, domain, neighbor, run, print.
        .add_plugins(CorePlugins)
        // DEM: atom data + insertion, Velocity Verlet, Hertz-Mindlin contact, rotation.
        .add_plugins(GranularDefaultPlugins)
        // Optional, config-driven physics — each inert unless its section is present.
        .add_plugins(GravityPlugin) // [gravity] body force
        .add_plugins(WallPlugin) // [[wall]] container / boundary faces
        .add_plugins(FixesPlugin); // [[addforce]] / [[freeze]] / [[viscous]] / ...

    // Everything case-specific comes from the TOML config path on argv.
    app.start();
}
