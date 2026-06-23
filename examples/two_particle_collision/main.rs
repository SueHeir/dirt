//! two_particle_collision — the smallest contact-physics demo.
//!
//! Two spheres are placed on the x-axis moving toward each other; they collide
//! once and rebound. This is the minimal assembly that exercises the
//! Hertz-Mindlin contact force: `CorePlugins` (infrastructure) +
//! `GranularDefaultPlugins` (atom data, Verlet integration, contact) reading a
//! `MaterialTable` from the `[[dem.materials]]` config section. No gravity, no
//! walls — just two particles and one contact.
//!
//! ```bash
//! cargo run --release --example two_particle_collision --no-default-features -- \
//!     examples/two_particle_collision/config.toml
//! ```
//!
//! With `restitution = 0.9` in the config, the post-collision relative speed is
//! ~0.9× the pre-collision relative speed (the COR is the realized COR — see
//! `dirt_atom::MaterialTable`). Inspect the VTP output to watch the rebound.

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins) // input, comm, domain, neighbor, run, print
        .add_plugins(GranularDefaultPlugins); // atom data + Verlet + Hertz-Mindlin contact

    // Particle positions, velocities, material, and run length all come from the
    // TOML config. The two particles are inserted from a CSV so their initial
    // approach velocities are exact.
    app.start();
}
