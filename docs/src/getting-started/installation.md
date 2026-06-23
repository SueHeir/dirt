# Installation & Building

## Prerequisites

- **Rust** (stable, 2021 edition or newer) — install from
  [rustup.rs](https://rustup.rs).
- **A C compiler and an MPI library** *only if* you build with the default
  `mpi_backend` feature. For single-process runs you can skip MPI entirely with
  `--no-default-features` (see below).

You do **not** need to check out GRASS or SOIL yourself. DIRT pulls its
[`soil`](https://github.com/SueHeir/soil) and
[`grass`](https://github.com/SueHeir/grass) dependencies from GitHub
automatically during the build.

## Clone and build

```bash
git clone https://github.com/SueHeir/dirt
cd dirt
cargo build --release --no-default-features
```

`--no-default-features` disables the MPI backend and builds a single-process
binary — the quickest way to get running. Drop it once you have an MPI toolchain
and want multi-rank parallelism.

## Run an example

DIRT ships its examples as Cargo examples. Each one is a thin `main.rs` plus a
TOML config:

```bash
cargo run --release --example hopper --no-default-features -- examples/hopper/config.toml
```

The trailing argument is the config file; the engine reads all simulation
parameters from it. See [Your First Simulation](./first-simulation.md) for a
walk-through of what that command actually does.

## MPI builds

The default feature set enables `mpi_backend`, which links an MPI library for
multi-rank domain-decomposed runs:

```bash
cargo build --release           # mpi_backend is on by default
mpirun -np 4 ./target/release/examples/hopper examples/hopper/config.toml
```

If the build fails to find MPI, either install a system MPI (OpenMPI or MPICH)
or fall back to `--no-default-features` for serial runs.

## Using DIRT as a library

To build your own simulation binary, depend on `dirt_core` and pull the prelude:

```toml
# Cargo.toml
[dependencies]
dirt_core = { git = "https://github.com/SueHeir/dirt" }
```

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
       .add_plugins(GranularDefaultPlugins);
    app.start();
}
```

> **Note:** DIRT depends on `soil` and `grass` via git, so its crates cannot be
> published to crates.io (which forbids git dependencies). Depend on it by git
> URL, as above.
