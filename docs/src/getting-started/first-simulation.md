# Your First Simulation

We'll run the **hopper** example: particles fall into a funnel under gravity,
settle, and then a blocker wall is removed so they discharge through the exit.
It exercises gravity, walls, contacts, and a runtime stage transition — a
representative slice of DIRT.

```bash
cargo run --release --example hopper --no-default-features -- examples/hopper/config.toml
```

You'll see the filling phase run, a line announcing that the particles have
settled and the blocker is being removed, then the flowing phase.

## Anatomy of the binary

The whole simulation program is a handful of plugins. Here is the example's
`main.rs`, lightly trimmed:

```rust
use dirt_core::prelude::*;

#[derive(Clone, Debug, PartialEq, Default, StageEnum)]
enum Phase {
    #[default]
    #[stage("filling")]
    Filling,
    #[stage("flowing")]
    Flowing,
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)            // config, comm, domain, neighbors, run loop, I/O
        .add_plugins(GranularDefaultPlugins) // Hertz–Mindlin contact + Velocity Verlet
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(StatesPlugin::new(Phase::Filling, ParticleSimScheduleSet::PostFinalIntegration))
        .add_plugins(StageAdvancePlugin::<Phase>::new(ParticleSimScheduleSet::PostFinalIntegration));

    // A custom system: watch the kinetic energy, and when the bed is quiet,
    // remove the blocker wall and advance to the flowing stage.
    app.add_update_system(
        check_settled.run_if(in_state(Phase::Filling)),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();
}
```

Three things to notice:

1. **Plugins compose the physics.** `CorePlugins` brings the infrastructure;
   `GranularDefaultPlugins` brings the contact force and integrator;
   `GravityPlugin` and `WallPlugin` add the rest. You assemble the simulation by
   choosing plugins.
2. **Stages are a state machine.** `Phase` is a two-state machine
   (`Filling` → `Flowing`). `StatesPlugin` and `StageAdvancePlugin` wire it into
   the scheduler. This is GRASS framework machinery, re-exported through DIRT.
3. **A system is just a function.** `check_settled` takes typed resources as
   arguments (`Res<Atom>`, `ResMut<Walls>`, …); the scheduler injects them. It
   runs only while `in_state(Phase::Filling)`.

Here is that system — your first piece of custom physics control:

```rust
fn check_settled(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    mut walls: ResMut<Walls>,
    mut next_state: ResMut<NextState<Phase>>,
) {
    let step = run_state.total_cycle;
    if step < 1000 || step % 100 != 0 {
        return; // give particles time to move, then check periodically
    }

    let nlocal = atoms.nlocal as usize;
    let local_ke: f64 = (0..nlocal)
        .map(|i| {
            let v = atoms.vel[i];
            0.5 * atoms.mass[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        })
        .sum();
    let global_ke = comm.all_reduce_sum_f64(local_ke); // sum across MPI ranks

    if global_ke < 1e-5 {
        walls.deactivate_by_name("blocker");
        next_state.set(Phase::Flowing);
    }
}
```

Note `comm.all_reduce_sum_f64` — the kinetic energy is summed across MPI ranks,
so the same code works serial or parallel. The substrate handles the
decomposition; you write ordinary physics.

## Where the physics lives

The `main.rs` says *what* to simulate. The TOML config says *with what numbers* —
domain size, materials, particle insertion, walls, gravity. The next chapter,
[Anatomy of a Config File](./config-anatomy.md), dissects the hopper config
section by section.
