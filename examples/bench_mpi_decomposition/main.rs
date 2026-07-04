//! bench_mpi_decomposition — MPI domain-decomposition correctness for DEM.
//!
//! Proves that running a contact-rich granular gas under a **>1 processor**
//! domain decomposition (`processors_x/y/z`) reproduces the single-rank
//! (`1×1×1`) trajectory to the floating-point associativity floor, that global
//! linear momentum and kinetic energy evolve identically regardless of the
//! decomposition, and that no atom is lost or duplicated as atoms migrate and
//! are ghost-exchanged across rank boundaries.
//!
//! Like `bench_restart_determinism`, this binary is deliberately just the stock
//! plugin wiring — `CorePlugins` (input, comm, domain, neighbor, run, print,
//! which owns the MPI-gathered `[dump]` writer) plus `GranularDefaultPlugins`
//! (atom data + Verlet + Hertz–Mindlin contact with the per-contact tangential
//! spring history that must survive migration). Everything that varies between
//! the runs of the benchmark — the `[comm]` processor grid, the number of
//! steps, the `[dump]` cadence, and the output directory — is driven entirely
//! from the TOML config, so this binary needs **no** bench-specific code. The
//! driver (`sweep.py`) composes the per-decomposition configs and gates the
//! cross-rank comparison.
//!
//! The physics is a small, dense, frictional granular gas in a fully periodic
//! box with gravity off: ~400 glass spheres given a seeded Gaussian velocity
//! field collide and cool (restitution < 1). At any step there are many live
//! frictional contacts carrying tangential spring history, and — with `2×1×1`
//! or `2×2×1` — a steady stream of atoms crossing rank boundaries, so ghost
//! exchange, migration, and per-contact history transfer are all exercised.
//!
//! ```bash
//! # Standalone single-rank smoke test (default control block in config.toml):
//! cargo run --release --example bench_mpi_decomposition -- \
//!     examples/bench_mpi_decomposition/config.toml
//!
//! # 2-rank MPI run along x (the decomposition under test):
//! cargo build --release --example bench_mpi_decomposition
//! mpiexec -n 2 target/release/examples/bench_mpi_decomposition \
//!     examples/bench_mpi_decomposition/config.toml
//!
//! # Full gated benchmark (composes 1×1×1 vs 2×1×1 vs 2×2×1 and compares):
//! python3 examples/bench_mpi_decomposition/sweep.py
//! ```

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins) // input, comm, domain, neighbor, run, print (MPI-gathered dump)
        .add_plugins(GranularDefaultPlugins); // atom data + Verlet + Hertz–Mindlin contact

    // Geometry, material, insertion (seeded), the [comm] processor grid, and the
    // run/dump control all come from the TOML config passed on the command line.
    app.start();
}
