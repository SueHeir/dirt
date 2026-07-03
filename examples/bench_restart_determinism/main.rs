//! bench_restart_determinism — restart continuity + run-to-run determinism.
//!
//! This example is deliberately identical in spirit to `two_particle_collision`:
//! it wires `CorePlugins` (input, comm, domain, neighbor, run, print — which
//! owns the `[restart]`/`[dump]` machinery) plus `GranularDefaultPlugins`
//! (atom data + Verlet + Hertz-Mindlin contact, including the per-contact
//! `ContactHistoryStore` tangential spring history). Everything that varies
//! between the runs of the benchmark — number of steps, the `[restart]`
//! interval/read flags, the `[dump]` interval, and the output directory — is
//! driven entirely from the TOML config, so this binary needs **no** bench
//! specific code. The driver (`sweep.py`) composes the per-run configs.
//!
//! The physics is a small frictional granular gas in a fully periodic box:
//! ~300 glass spheres given a seeded Gaussian velocity field collide and cool.
//! A moderately dense, frictional, chaotic system is the point — at any given
//! step there are live frictional contacts carrying tangential spring history,
//! so a restart that dropped that history (or any per-atom state) would visibly
//! diverge. Insertion is seeded (`[[particles.insert]] seed`), contact/Verlet
//! integration is deterministic, and single-rank output is written in atom
//! order, so two identical runs are bit-identical and a checkpoint→resume run
//! reproduces the uninterrupted trajectory.
//!
//! ```bash
//! # Standalone (uses the default control block at the end of config.toml):
//! cargo run --release --example bench_restart_determinism \
//!     --no-default-features --features precision-double -- \
//!     examples/bench_restart_determinism/config.toml
//!
//! # Full gated benchmark (composes the checkpoint / resume / twin runs):
//! python3 examples/bench_restart_determinism/sweep.py
//! ```

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins) // input, comm, domain, neighbor, run, print (restart + dump)
        .add_plugins(GranularDefaultPlugins); // atom data + Verlet + Hertz-Mindlin contact

    // Geometry, material, insertion (seeded), and run/restart/dump control all
    // come from the TOML config passed on the command line.
    app.start();
}
