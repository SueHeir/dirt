# Hopper

2D slot hopper with angled funnel walls, gravity, and named simulation stages.

```
Cross-section (x-z plane, periodic in y):

  |                     |    side walls (x=0, x=0.04)
  |   particles here    |
  |                     |
   \                   /     angled funnel walls
    \                 /
     \               /
      \             /
       \___     ___/         funnel exit (1 cm opening)
       |  blocker  |         removable blocker wall (z=0.015)
       |___________|         floor (z=0)
```

**Filling stage:** 200 particles are inserted in the upper region and settle under gravity onto the angled funnel walls and blocker.

**Flowing stage:** Once the total kinetic energy drops below 1e-5 J (particles nearly stationary), the blocker wall is automatically removed, the run advances to the next `[[run]]` stage, and particles flow through the funnel exit to the floor. The flowing stage has its own step count and thermo interval.

This example demonstrates DIRT's Tier 2 (Rust API) by adding a custom system alongside the standard TOML-configured plugins.

### How `main.rs` works

```rust
app.add_update_system(
    check_settled.run_if(in_stage("filling")),
    ParticleSimScheduleSet::PostFinalIntegration,
);
```

The `check_settled` system is registered with a **run condition**: it only executes while the named `"filling"` stage is active. Once the run advances to `"flowing"`, the system is skipped. `ParticleSimScheduleSet::PostFinalIntegration` places it after the Velocity Verlet update each timestep.

The `check_settled` function itself is a regular system that declares its dependencies as function arguments — the scheduler injects them automatically:

- `Res<Atom>` — read-only access to particle data (velocities, masses)
- `Res<RunState>` — current timestep
- `Res<CommResource>` — MPI communicator for global reductions
- `ResMut<Walls>` — mutable access to wall definitions
- `ResMut<SchedulerManager>` — mutable access to request the next run stage

Every 100 steps (after an initial 1000-step warmup), it computes the total kinetic energy across all MPI ranks via `comm.all_reduce_sum_f64()`. When KE drops below the threshold, it deactivates the named `"blocker"` wall and requests the next configured run stage.

This pattern — TOML config for standard physics, custom Rust systems for runtime logic — is the core design of DIRT.

### Config: Named `[[run]]` stages

```toml
[[run]]
name = "filling"
steps = 1000000
thermo = 2000

[[run]]
name = "flowing"
steps = 1000000
thermo = 2000
```

Each `[[run]]` stage has a name that code can reference with `in_stage("...")`. Each stage can have its own step count, thermo interval, and output settings.

## Run

```bash
# Single-process
cargo run --example hopper -- examples/hopper/config.toml

# With MPI
cargo build --release --example hopper
mpiexec -n 4 ./target/release/examples/hopper examples/hopper/config.toml
```

## Parameters

| Parameter | Value |
|-----------|-------|
| Particles | 200 |
| Radius | 0.001 m |
| Density | 2500 kg/m^3 |
| Young's modulus | 8.7 GPa |
| Poisson ratio | 0.3 |
| Restitution | 0.3 |
| Friction | 0.5 |
| Gravity | -90.81 m/s^2 (z) |
| Domain | 0.04 x 0.02 x 0.08 m |
| Boundaries | Non-periodic x/z, periodic y |
| Funnel angle | ~67 deg from horizontal |
| Funnel exit | 1 cm opening at z = 0.015 m |
| Blocker wall | z = 0.015 m (removed when KE < 1e-5 J) |
| Filling stage | 1,000,000 steps, thermo 2000 |
| Flowing stage | 1,000,000 steps, thermo 2000 |
| KE check | Every 100 steps after step 1000 |

## Validation

`validate.py` checks physics sanity: no NaN/Inf, non-negative temperature, bounded energy (no explosion). Run via `./validate.sh` or directly:

```bash
python3 examples/hopper/validate.py
```
