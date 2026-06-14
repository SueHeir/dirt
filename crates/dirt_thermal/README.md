# dirt_thermal

Contact-based heat conduction for DEM particles in DIRT.

## What it does

Transfers thermal energy between overlapping particles in a Discrete Element Method simulation. Heat flows from hotter to cooler particles through the contact area, and each particle's temperature is integrated forward every timestep. Energy is conserved (the per-pair flux is antisymmetric). Temperature-bearing walls (plane, cylinder, sphere, region) also conduct heat into contacting particles, acting as infinite thermal reservoirs.

## Physics model

Heat flux between contacting particles *i* and *j*:

```text
Q = k · 2a · (T_j − T_i)
```

where:
- `k` — thermal conductivity (W/(m·K))
- `a = √(r_eff · δ)` — Hertzian contact radius
- `r_eff = (r_i · r_j) / (r_i + r_j)` — effective radius
- `δ = (r_i + r_j) − d` — overlap depth (`d` = center-to-center distance)

Temperature is then integrated: `T += dt · Q / (m · c_p)`, where `m` is the particle mass and `c_p` the specific heat.

## Key types

| Item | Description |
|------|-------------|
| `ThermalPlugin` | Registers `ThermalAtom` data and, when `[thermal]` is configured, the heat-conduction systems. Requires `DemAtomPlugin` (radii) and `NeighborPlugin` (contact pairs). |
| `ThermalConfig` | Config from `[thermal]`: `conductivity` (1.0), `specific_heat` (500.0), `initial_temperature` (300.0). |
| `ThermalAtom` | Per-atom `temperature` (K, forward-communicated) and `heat_flux` (W, reverse-communicated and zeroed each step). |
| `compute_heat_conduction` | `Force` system: accumulates particle-particle contact heat flux. |
| `compute_wall_heat_conduction` | `Force` system added when a `Walls` resource exists: conducts heat from temperature-bearing walls. |
| `integrate_temperature` | `PostFinalIntegration` system: advances temperatures from accumulated flux. |

If the `[thermal]` section is omitted, only `ThermalAtom` is registered (default values) and no heat-transfer systems run.

## TOML configuration

```toml
[thermal]
conductivity = 1.0          # W/(m·K)
specific_heat = 500.0       # J/(kg·K)
initial_temperature = 300.0 # K (optional, default: 300.0)
```

## Usage

```rust
use grass_app::prelude::*;
use dirt_thermal::ThermalPlugin;

let mut app = App::new();
app.add_plugins(ThermalPlugin);
app.start();
```

Access temperatures via the `AtomDataRegistry`:

```rust
let thermal = registry.expect::<ThermalAtom>("my_system");
let temp_i = thermal.temperature[i];
```

## License

MIT OR Apache-2.0
