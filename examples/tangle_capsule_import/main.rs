//! Loads a TANGLE adaptive-capsule DEM-BPM export into DIRT.
//!
//! TANGLE writes both the extended data file and a zero-step loading config:
//!
//! ```text
//! cargo run --release --example tangle_capsule_import -- \
//!   /path/to/fake_felted/output/dirt_capsules.toml
//! ```

use dirt_capsule::{CapsuleContactPlugin, CapsulePlugin};
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(dirt_core::soil_verlet::VelocityVerletPlugin::new())
        .add_plugins(CapsuleContactPlugin)
        .add_plugins(CapsulePlugin)
        .add_plugins(DemBondPlugin);
    app.start();
}
