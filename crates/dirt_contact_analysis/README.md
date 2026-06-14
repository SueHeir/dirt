# dirt_contact_analysis

Contact network analysis for DEM simulations: per-atom coordination number, rattler detection, per-contact CSV dumps, and fabric tensor computation.

## What It Does

`dirt_contact_analysis` provides a `ContactAnalysisPlugin` that reads a `[contact_analysis]` TOML config section and registers post-force systems over the neighbor list. In a single neighbor traversal it computes any combination of:

- **Coordination number**: count of active contacts (overlap > 0) per particle, stored via the `ContactAnalysis` AtomData and exposed to thermo (`coord_avg`, `coord_max`, `coord_min`) and as the per-atom dump scalar `coordination`.
- **Rattler detection**: particles with fewer than 4 contacts (lacking the *d* + 1 = 4 constraints for static equilibrium in 3D), reported to thermo as `n_rattlers` and `rattler_fraction`.
- **Per-contact CSV output**: geometric data (atom tags, overlap, contact point, contact normal) for every contact pair, written every `interval` steps.
- **Fabric tensor**: the symmetric second-order tensor *F_ij = (1/Nc) Σ n_i n_j* describing the directional distribution of contact normals; its six independent components plus the total `contacts` count are pushed to thermo.

## Key Types

| Type | Description |
|------|-------------|
| `ContactAnalysisPlugin` | Registers the analysis systems based on config flags |
| `ContactAnalysisConfig` | `[contact_analysis]` config: `interval`, `coordination`, `rattlers`, `fabric_tensor`, `file_prefix` |
| `ContactAnalysis` | Per-atom `AtomData` holding the coordination number |
| `ContactRecord` | One contact's geometry (tags, overlap, contact point, normal) |
| `ContactOutput` | Resource holding the current step's contact records |

## Configuration

All fields are optional with sensible defaults:

```toml
[contact_analysis]
interval = 1000          # dump per-contact CSV every N steps (0 = disabled, default 0)
coordination = true      # compute per-atom coordination number (default false)
rattlers = true          # detect rattler particles with < 4 contacts (default false)
fabric_tensor = true     # output fabric tensor components to thermo (default false)
file_prefix = "contact"  # CSV filename prefix (default "contact")
```

## Usage

```rust
use dirt_contact_analysis::ContactAnalysisPlugin;

app.add_plugin(ContactAnalysisPlugin);
```

Thermo output (on thermo intervals):

- `coord_avg`, `coord_max`, `coord_min` — coordination statistics
- `n_rattlers`, `rattler_fraction` — when `rattlers` is enabled
- `fabric_xx`, `fabric_yy`, `fabric_zz`, `fabric_xy`, `fabric_xz`, `fabric_yz`, `contacts` — when `fabric_tensor` is enabled

When `interval > 0`, CSV records are written to `<output_dir>/contact/<file_prefix>_<step>_rank<rank>.csv` with columns `i_tag`, `j_tag`, `overlap`, `cx`, `cy`, `cz`, `nx`, `ny`, `nz`.

## License

MIT OR Apache-2.0
