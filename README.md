# DIRT: Discrete-element Interaction-Resolved Toolkit

**A modular, parallel discrete-element-method solver for granular materials,
rigid clumps, bonded particles, and fibers.**

<!-- disclaimer-banner -->
> **Research-software status:** This ecosystem is AI-authored and under active
> evaluation. DIRT's DEM claims are accompanied by reproducible analytical,
> cross-code, or empirical evidence; repositories prefixed `dev_` are
> experimental method demonstrations outside the author's domain expertise. See
> [DISCLAIMER.md](DISCLAIMER.md) and
> [examples/VALIDATION.md](examples/VALIDATION.md).
<!-- /disclaimer-banner -->

DIRT is a ground-up Rust DEM implementation built around replaceable physics
plugins. It provides particle and material data, contact mechanics, walls,
bonds, rigid multisphere bodies, loading, diagnostics, distributed execution,
and reproducible validation examples as one usable simulation code.

## What DIRT provides

### Particles and materials

- Spherical particles with per-type material properties and material-pair
  mixing.
- Fixed, uniform, and discrete particle-radius distributions.
- Random, lattice, rate-based, CSV, LAMMPS-data, and LAMMPS-dump insertion.
- Rigid multisphere clumps for representing nonspherical bodies.
- Translational and rotational dynamics.

### Contact mechanics

- Hertz–Tsuji nonlinear normal contact.
- Hooke linear spring-dashpot contact.
- MDR elastic-plastic normal contact.
- Mindlin tangential history, history-free tangential response, and
  LAMMPS-style Mindlin unloading-rescale variants.
- Coulomb sliding with rolling and twisting resistance.
- JKR and DMT adhesion, SJKR cohesion, and Willett liquid-bridge cohesion where
  supported by the selected contact and boundary path.

Not every model combination is meaningful or implemented. The contact
documentation records important boundaries, including adhesion support by
contact model and differences between plane and curved walls:
[Contact Models](docs/src/physics/contact.md).

### Boundaries, loading, and control

- Plane, cylinder, sphere, cone, and region-surface walls.
- Tangential wall friction and rolling/twisting resistance.
- Moving, named, activated, and servo-controlled boundaries.
- Gravity, Cundall damping, viscous damping, prescribed motion, freeze/pin,
  integration limiting, and applied-force fixes.
- Periodic, fixed, shrink-wrapped, deforming, and Lees–Edwards domains supplied
  by the underlying particle infrastructure.

### Bonded particles and fibers

- Beam-like normal, shear, bending, and twisting bond response.
- Elastic and elastoplastic axial and bending behavior.
- Stress- and strain-based breakage, including statistical strength models.
- Bond creation, sintering, and bonded-fiber configurations.

### Diagnostics and output

- Per-contact geometry and force output.
- Coordination number, fabric tensor, and rattler analysis.
- Measurement planes for particle counts, mass flow, flux, and profiles.
- Thermodynamic output, VTP visualization, text/binary dumps, and restart files.
- Examples for hopper flow, shear, impact, granular cooling, clumps, bonded
  specimens, and other DEM workflows.

## Two ways to run DIRT

### Assemble a simulation in Rust

A DIRT application is a set of plugins chosen for one problem:

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)            // domain, particles, comm, neighbors, I/O
       .add_plugins(GranularDefaultPlugins) // contact, rotation, Velocity Verlet
       .add_plugins(GravityPlugin)
       .add_plugins(WallPlugin);
    app.start();
}
```

A plugin registers state and scheduled systems. Custom physics or measurements
are ordinary Rust functions added at the appropriate typed DEM phase. The
[Your First Simulation](docs/src/getting-started/first-simulation.md) tutorial
walks through a complete example and adds a custom system.

### Run a declarative scenario

The prebuilt `run` example assembles the common plugin stack and reads geometry,
materials, insertion, walls, forces, and loading from TOML:

```console
cargo run --release --example run -- examples/run/pour_settle.toml
```

This is convenient for shipped scenarios and parameter sweeps because the same
binary can run many configurations without recompiling. TOML selects existing
capabilities; a new physical model still belongs in a Rust plugin.

See [Run from a Config](docs/src/getting-started/run-from-config.md).

## Parallel particle execution

DIRT uses SOIL for the spatial machinery beneath its DEM physics:

- MPI domain decomposition and local particle ownership;
- particle migration and ghost exchange;
- forward state replication and reverse force/torque accumulation;
- bin-based neighbor construction and rebuild decisions;
- restart-safe registered particle data;
- double, mixed, and single-precision modes.

The same scheduled DEM systems run in a single process or over distributed
subdomains. MPI is enabled by default; a no-MPI build is available for local use.

## Scientific evidence

DIRT keeps scientific validation separate from numerical and software
verification:

- [Scientific Validation](examples/VALIDATION.md) contains analytical,
  experimental/empirical, and independently executed cross-code comparisons.
- [Numerical and Software Verification](examples/VERIFICATION.md) contains
  convergence, reproducibility, MPI, restart, API, and build checks.

Cross-code agreement with LAMMPS tests consistency under a shared model; it is
not automatically evidence that the model represents a real material. Empirical
scaling checks are likewise identified separately from closed-form validation.
Known failures and withheld comparisons remain visible instead of being tuned
away.

The evidence covers important parts of the code, including elastic and damped
impact, tangential response, rotational resistance, adhesion/cohesion,
bonded-particle mechanics, granular flow, and distributed execution. It does not
make every implemented feature experimentally validated. Consult the ledgers
before relying on a particular model combination.

## Install

You need stable [Rust](https://rustup.rs/). DIRT pulls GRASS and SOIL during the
build.

For a single-process build with double precision:

```console
git clone https://github.com/SueHeir/dirt
cd dirt
cargo run --release --example hello_bed \
  --no-default-features --features precision-double \
  -- examples/hello_bed/config.toml
```

The default features enable MPI and double precision. With an MPI toolchain:

```console
cargo build --release
mpirun -np 4 ./target/release/examples/hopper examples/hopper/config.toml
```

See [Installation and Building](docs/src/getting-started/installation.md) for
library dependencies, feature selection, and MPI requirements.

## Crate map

`dirt_core` is the batteries-included umbrella crate. The other crates expose
individual pieces for applications that need direct control:

| Crate | Role |
|---|---|
| [`dirt_core`](crates/dirt_core/README.md) | prelude and core/default plugin groups |
| [`dirt_atom`](crates/dirt_atom/README.md) | DEM particle data, materials, radius distributions, and insertion |
| [`dirt_granular`](crates/dirt_granular/README.md) | normal/tangential contact, adhesion, rolling/twisting, and rotation |
| [`dirt_wall`](crates/dirt_wall/README.md) | wall geometry, contact response, motion, and servo control |
| [`dirt_bond`](crates/dirt_bond/README.md) | bonded-particle beams, plasticity, breakage, and sintering |
| [`dirt_clump`](crates/dirt_clump/README.md) | rigid multisphere composites |
| [`dirt_fixes`](crates/dirt_fixes/README.md) | gravity, damping, constraints, motion, and applied forces |
| [`dirt_contact_analysis`](crates/dirt_contact_analysis/README.md) | contact records, coordination number, fabric tensor, and rattlers |
| [`dirt_measure_plane`](crates/dirt_measure_plane/README.md) | measurement planes for counts, flow, flux, and profiles |
| [`dirt_schedule`](crates/dirt_schedule/README.md) | typed DEM scheduler labels shared by plugins |
| [`dirt_test_utils`](crates/dirt_test_utils/README.md) | shared test helpers |

The mdBook under `docs/` contains the complete user and physics documentation.

## Ecosystem

DIRT is the DEM tier of a one-way dependency stack:

```text
GRASS   scientific application framework
  └── SOIL   distributed particle infrastructure
        └── DIRT   discrete-element-method physics and applications
```

- [GRASS](https://github.com/SueHeir/grass) provides Apps, scheduling, plugins,
  lifecycle, configuration, and communication abstractions.
- [SOIL](https://github.com/SueHeir/soil) provides particle storage, domains,
  migration, ghost communication, and neighbor search.
- DIRT owns the DEM-specific materials, contact laws, boundaries, bonds, clumps,
  loading, diagnostics, and scientific evidence.

Because these layers expose typed state and scheduled behavior, DIRT can be used
as one component of a larger application. Such application-specific composition
is optional and is not part of DIRT's core scientific claim.


## License

MIT OR Apache-2.0
