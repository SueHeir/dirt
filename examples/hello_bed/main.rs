//! hello_bed — the smallest complete DIRT simulation.
//!
//! Drops a cloud of spheres into a closed box under gravity and lets them
//! settle into a bed. This is the minimal `CorePlugins + GranularDefaultPlugins`
//! assembly with gravity and six container walls — no stages, no runtime wall
//! control, no analysis. Start here, then graduate to the `hopper` example for
//! multi-stage runs and runtime wall removal.
//!
//! ```bash
//! cargo run --release --example hello_bed --no-default-features -- \
//!     examples/hello_bed/config.toml
//! ```
//!
//! The `--no-default-features` flag turns off the MPI backend for a simple
//! single-process run.

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins) // input, comm, domain, neighbor, run, print
        .add_plugins(GranularDefaultPlugins) // atom data, insertion, Verlet, Hertz-Mindlin, rotation
        .add_plugins(GravityPlugin) // [gravity] body force
        .add_plugins(WallPlugin); // [[wall]] container faces

    // Everything else — particle insertion, walls, gravity, run length — is
    // driven entirely by the TOML config. The app reads the config path from
    // the command line.
    app.start();
}
