# GRASS → SOIL → DIRT

DIRT is deliberately split across three repositories, one per tier. Lower tiers
never depend on higher ones.

```
GRASS    framework: App, Plugin, Scheduler, IO, coupling      (no particles)
  └─ SOIL   substrate: Atom, domain decomposition, comm, neighbor lists   (no physics)
       └─ DIRT   DEM physics: contact, bonds, walls, clumps   ← you are here
```

| tier | repo | owns | knows nothing about |
|---|---|---|---|
| **GRASS** | [grass](https://github.com/SueHeir/grass) | `App`, `Plugin`, scheduler, TOML config, MPI coupling | particles, physics |
| **SOIL** | [soil](https://github.com/SueHeir/soil) | base `Atom`, domain decomposition, ghost comm, migration, neighbor lists | contact forces, bonds, damage |
| **DIRT** | [dirt](https://github.com/SueHeir/dirt) | Hertz–Mindlin contact, bonds, walls, clumps, conduction | — |

## Why split it this way

The seam between SOIL and DIRT is a single contract: a physics tier registers
its per-particle state as an `AtomData` column, and the substrate then carries
that state through every migration, ghost exchange, and restart automatically.
Because that contract is physics-agnostic, the *same* substrate can carry a
completely different physics — SPH, peridynamics, your own force law — with no
change to SOIL.

That is the subject of the other two books:

- **[SOIL book](https://sueheir.github.io/soil)** — the `AtomData` contract and a
  tutorial for writing your own particle-physics tier on the substrate.
- **[GRASS book](https://sueheir.github.io/grass)** — the `App`/`Plugin`/scheduler
  model and how to build a solver from scratch.

## Two plugin groups

A DIRT binary is assembled from plugins. Two groups do the heavy lifting, and the
split between them is the single most important thing to understand about the
assembly:

```rust
use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)             // I/O, comm, domain, neighbors, run, print
       .add_plugins(GranularDefaultPlugins); // material, insertion, Verlet, contact, rotation
    app.start();
}
```

- **`CorePlugins`** (from `dirt_core`) is the infrastructure half: config input,
  the communication backend, domain decomposition, neighbor lists, atom groups,
  the run loop, and output. It is physics-agnostic — it carries `AtomData`
  columns through migration and ghost exchange but never decides a force.
- **`GranularDefaultPlugins`** (from `dirt_granular`) is the DEM-physics half:
  per-atom material properties, particle insertion, **Velocity Verlet**, the
  Hertz–Mindlin contact force, and rotational dynamics.

> **`CorePlugins` alone never moves a particle.** Velocity Verlet — the
> translational integrator — lives in `GranularDefaultPlugins`, *not* in
> `CorePlugins`. A binary that adds only `CorePlugins` will read its config,
> build neighbor lists, and print thermo, but every particle stays exactly where
> it started. This is intentional (it is useful for testing infrastructure in
> isolation), but it means forgetting `GranularDefaultPlugins` produces a run
> that looks alive in the logs yet does zero physics.

Everything else — gravity, walls, bonds, clumps, diagnostics — is an opt-in
plugin bolted on after these two. The prelude (`dirt_core::prelude::*`) re-exports
the whole stack, including `GranularDefaultPlugins`, so you import from one place.

## One timestep, end to end

The run loop advances each step through an ordered set of schedule phases. A
DEM step with Velocity Verlet, gravity, and contact runs roughly like this:

1. **Zero** — force and torque accumulators are cleared (the `#[zero]` `AtomData`
   columns). Each step starts from a clean slate.
2. **Initial integration** — the first Verlet half-kick: velocities advance by
   `½ dt · a`, then positions advance by `dt · v`. Clump bodies and per-sphere
   angular state take their half-kick here too.
3. **Forward (ghost exchange)** — owner ranks push the updated positions (and
   `#[forward]` columns such as `omega`, `body_id`) out to their ghost copies on
   neighboring ranks, so every rank sees a consistent picture of the particles
   straddling its boundary.
4. **Forces** — the neighbor list is traversed and the contact law writes forces
   and torques; gravity and other body forces add in; walls and bonds contribute.
   This is where the bulk of DIRT's physics happens.
5. **Reverse (ghost accumulation)** — the forces computed on ghost copies are
   summed back to the owner rank (the `#[reverse]` columns, `torque` among them),
   so each owner ends with the *total* force on its atoms.
6. **Final integration** — the second Verlet half-kick completes the velocity
   update with the new acceleration. Post-force fixes (freeze, viscous, setforce)
   apply here, and velocity caps (`nve_limit`) clamp after.
7. **Exchange / output** — atoms that left their subdomain migrate to their new
   owner (carrying every `AtomData` column), the neighbor list is rebuilt when
   due, and thermo/dump output is written on its interval.

The seam between SOIL and DIRT runs straight through this loop: steps 1, 3, 5,
and 7 are the substrate moving `AtomData` columns around; steps 2, 4, and 6 are
DIRT (and the Verlet integrator) deciding what those columns should be. Neither
tier needs to know the other's internals — only the `AtomData` contract.
