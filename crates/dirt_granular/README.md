# dirt_granular

Granular contact physics for [DIRT](https://github.com/SueHeir/dirt): Hertz/Hooke normal contact, Mindlin tangential friction, rolling and twisting resistance, adhesion, rotational dynamics, and granular temperature output.

## What it does

Core DEM physics for spherical particle simulations:

- **Normal contact**: Hertz (nonlinear `F_n ∝ δ^1.5`, default) or Hooke (linear `F_n ∝ δ`), with viscoelastic damping
- **Tangential friction**: Mindlin spring-history with Coulomb cap and viscous damping; per-contact spring rotated into the tangent plane each step
- **Rolling resistance**: constant torque (default) or SDS spring-dashpot-slider
- **Twisting friction**: constant torque (default) or SDS
- **Adhesion / cohesion**: JKR (extended range), DMT (contact-only), or SJKR area-based cohesion
- **Rotational dynamics**: quaternion velocity Verlet (`I = 2/5 m r²` for solid spheres)
- **Output**: granular temperature (velocity fluctuation) time series

Overlap is guarded: pairs exceeding `LARGE_OVERLAP_WARN_THRESHOLD` warn (forces still computed, capped), and more than `MAX_OVERLAP_WARNINGS` per step aborts with an actionable message.

## Key types

| Item | Role |
| --- | --- |
| `GranularDefaultPlugins` | Default DEM granular plugin group (contact + rotation + insertion + temp) |
| `HertzMindlinContactPlugin` | Fused normal + tangential contact (primary code path) |
| `RotationalDynamicsPlugin` | Quaternion velocity Verlet for angular DOF |
| `GranularTempPlugin` | Granular temperature file output |

Modules: `contact` (fused Hertz–Mindlin + Hooke), `tangential` (`ContactHistoryStore`), `rotational`, `granular_temp`.

## TOML configuration

```toml
[[dem.materials]]
name = "glass"
youngs_modulus = 8.7e9      # Pa
poisson_ratio = 0.3
restitution = 0.95          # coeff. of restitution (0–1)
friction = 0.4              # sliding friction μ
rolling_friction = 0.1      # rolling friction μ_r
cohesion_energy = 0.0       # J/m² (SJKR, 0 = disabled)
surface_energy = 0.0        # J/m² (JKR/DMT, 0 = disabled)

[dem]
contact_model = "hertz"     # "hertz" (default) or "hooke"
adhesion_model = "jkr"      # "jkr" (default) or "dmt"
rolling_model = "constant"  # "constant" (default) or "sds"
twisting_model = "constant" # "constant" (default) or "sds"
```

## Usage

```rust
use dirt_core::prelude::*;

let mut app = App::new();
app.add_plugins(CorePlugins).add_plugins(GranularDefaultPlugins);
app.start();
```

`GranularDefaultPlugins` adds per-atom material properties, particle insertion, velocity Verlet, Hertz–Mindlin contact, rotational dynamics, and granular temperature output. It does not include infrastructure plugins; pair it with `CorePlugins` for input, comm, domain, neighbor, run, and print.

Part of the [DIRT](https://github.com/SueHeir/dirt) workspace.

## License

MIT OR Apache-2.0
