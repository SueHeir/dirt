//! clump_dimer_drop — a minimal end-to-end rigid-body (clump) simulation.
//!
//! Drops a handful of two-sphere **dimer** clumps into a closed box under
//! gravity and lets them settle. Each dimer is a single rigid body: its two
//! sub-spheres detect contact like ordinary atoms, but their forces are
//! aggregated onto the body, which integrates translation + rotation and then
//! reconstructs the sub-sphere kinematics. This is the smallest assembly that
//! exercises `ClumpPlugin` end to end.
//!
//! Plugin stack:
//! - `CorePlugins` — input, comm, domain, neighbor lists, run loop, output.
//! - `GranularDefaultPlugins` — DEM atom data, particle insertion, velocity
//!   Verlet, Hertz–Mindlin contact, rotational dynamics.
//! - `GravityPlugin` — `[gravity]` body force.
//! - `WallPlugin` — the six `[[wall]]` container faces.
//! - `ClumpPlugin` — loads `[[clump.definitions]]`, inserts `[[clump.insert]]`
//!   bodies, aggregates sub-sphere forces, and integrates the rigid bodies.
//!
//! `ClumpPlugin` depends on `dirt_atom::DemAtomPlugin`, which
//! `GranularDefaultPlugins` already provides, so the order below is fine.
//!
//! Run (single process — MPI backend off):
//!
//! ```bash
//! cargo run --release --example clump_dimer_drop --no-default-features -- \
//!     examples/clump_dimer_drop/config.toml
//! ```

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(ClumpPlugin);

    // Clump definitions, insertion, walls, gravity, and run length are all
    // driven by the TOML config passed on the command line.
    app.start();
}
