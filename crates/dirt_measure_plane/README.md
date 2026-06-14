# dirt_measure_plane

General-purpose measurement plane plugin for counting particle crossings and mass flow rates in DIRT simulations.

## What it does

A **measurement plane** is an infinite plane defined by a point and a normal vector. Each timestep, the plugin computes the signed distance of every local particle from the plane and detects when a particle crosses in the positive-normal direction. For each configured plane it tracks:

- **Cumulative crossing count** — total particles crossing since simulation start
- **Mass flow rate** — mass per unit time, averaged over `report_interval`
- **Crossing rate** — particles per unit time, averaged over `report_interval`

This is useful for measuring throughput in hoppers, chutes, conveyors, and other granular flows where you need the particle and mass throughput across a cross-section.

## Crossing detection algorithm

For each particle (tracked by tag):
1. Compute the signed distance `d = (pos - point) · normal`.
2. If the previous distance `d_prev ≤ 0` and the current `d > 0`, count a positive-direction crossing.
3. Reverse crossings (positive → negative) are ignored.

Counts and mass are summed across MPI ranks via an all-reduce at each plane's `report_interval`.

## Key types

| Item | Description |
|------|-------------|
| `MeasurePlanePlugin` | Plugin that registers crossing-detection and reporting systems |
| `MeasurePlaneDef` | TOML configuration for a single `[[measure_plane]]` block |
| `MeasurePlanes` | ECS resource holding the per-plane runtime state |

The plugin registers its systems at `ParticleSimScheduleSet::PostFinalIntegration`. If no `[[measure_plane]]` blocks are configured, only an empty `MeasurePlanes` resource is inserted.

## TOML configuration

```toml
[[measure_plane]]
name = "outlet"          # Unique name; used in thermo output keys
point = [0.1, 0.0, 0.0]  # Any point on the plane [length units]
normal = [1.0, 0.0, 0.0] # Outward normal (auto-normalized)
report_interval = 1000   # Averaging window in timesteps (default: 1000)
```

Multiple `[[measure_plane]]` blocks can be defined; each plane tracks crossings independently.

## Thermo output keys

For each plane named `<name>`:
- `crossings_<name>` — total cumulative crossing count
- `flow_rate_<name>` — mass flow rate (mass/time), averaged over `report_interval`
- `cross_rate_<name>` — particle crossing rate (1/time), averaged over `report_interval`

## Usage

```rust
use dirt_measure_plane::MeasurePlanePlugin;
use grass_app::prelude::*;

let mut app = App::new();
app.add_plugin(MeasurePlanePlugin);
// ... add other plugins, configure planes in TOML, and run
```

## License

MIT OR Apache-2.0
