//! measure_plane_throughput — count particles falling through a horizontal plane.
//!
//! A column of spheres is dropped under gravity inside a closed box. A single
//! horizontal [`MeasurePlanePlugin`] plane near the bottom counts how many
//! particles cross it (in the +z... actually −z, see below) and reports the
//! mass/particle throughput to thermo every `report_interval` steps.
//!
//! Because the measurement plane is a **directional gate** (it counts crossings
//! *with* its normal only), the plane normal is set to point **downward**
//! `[0, 0, -1]` so a particle falling under gravity crosses in the counted
//! direction exactly once. See the `dirt_measure_plane` crate docs for the full
//! caveat list (oscillating particles recount, etc.).
//!
//! Watch the thermo columns `crossings_outlet`, `flow_rate_outlet`, and
//! `cross_rate_outlet`.
//!
//! ```bash
//! cargo run --release --example measure_plane_throughput --no-default-features -- \
//!     examples/measure_plane_throughput/config.toml
//! ```

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins) // input, comm, domain, neighbor, run, print
        .add_plugins(GranularDefaultPlugins) // atom data, Verlet, Hertz-Mindlin contact
        .add_plugins(GravityPlugin) // [gravity] body force
        .add_plugins(WallPlugin) // [[wall]] container faces
        .add_plugins(MeasurePlanePlugin); // [[measure_plane]] throughput gate

    // All geometry, the material, insertion, the measurement plane, and the run
    // length are defined in the TOML config.
    app.start();
}
