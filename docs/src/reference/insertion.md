# Particle Insertion

DIRT seeds a simulation with `[[particles.insert]]` blocks — random packings,
rate-based trickle feeds, or particles loaded from a file. The mechanism behind
all three is the **born-in-owner** model, which is what makes a packing
identical whether you run on 1 rank or 64.

## The born-in-owner model

Insertion is **not** done on one rank and scattered. Instead it runs on *every*
MPI rank, and every rank seeds its RNG with the **same** `seed`, so all ranks
generate the bit-identical candidate stream — the same positions, radii,
velocities, and tags, in the same order. Each rank then keeps **only** the
candidates whose position falls inside its own subdomain, tested with a
half-open interval (`low ≤ pos < high`) that matches the substrate's ownership
rule exactly.

Three properties follow:

- **Born in its owner.** Each atom is materialized only on the rank that owns
  its position, so the first post-insertion exchange never has to migrate it —
  no insertion-time communication is needed.
- **Exactly-once.** The half-open convention guarantees every position is
  claimed by one rank, never two and never zero — no duplicated or dropped atoms
  at subdomain boundaries.
- **Deterministic across rank counts.** Because the candidate stream depends
  only on `seed` (not on the decomposition), the *global* packing is identical
  whether you run on 1 rank or 64 — only the partitioning of that packing
  changes. Overlap rejection uses a global spatial hash replicated on every
  rank, so packings are reproducible run-to-run and rank-count to rank-count.

File-based insertion follows the same rule: the file is parsed on every rank and
filtered to the local subdomain, with tags advanced identically so they stay
globally consistent.

## The three insertion modes

A `[[particles.insert]]` block selects a mode by its fields:

- **Random** (default, `source = "random"` or omitted) — drop `count` particles
  into a `region`. Requires `material`, `count`, `radius`, `density`.
- **Rate-based** — random insertion with a `rate` field present; registers for
  periodic insertion (trickle feed). The candidate stream is derived from a
  step-derived seed, so new atoms are still born in their owner.
- **File-based** (`source = "file"`) — load positions from CSV or a LAMMPS dump
  / data file. Requires `file` and `format`.

```toml
[[particles.insert]]
material = "glass"        # must match a [[dem.materials]] name
count = 200
radius = 0.001            # fixed value, or a distribution (see below)
density = 2500.0          # kg/m³
velocity_z = -1.0         # directional initial velocity
region = { type = "block", min = [0.005, 0.0, 0.055], max = [0.035, 0.02, 0.075] }
seed = 0                  # deterministic insertion RNG (default 0)
```

## The three modes in detail

**Random** drops `count` particles into `region` at setup, rejecting overlaps
against a global spatial hash. **Rate-based** is triggered by a `rate` field and
trickles particles in over time:

| Key | Default | Meaning |
|---|---|---|
| `rate` | — (triggers mode) | Particles inserted per interval. |
| `rate_interval` | `1` | Insert every N timesteps. |
| `rate_start` / `rate_end` | `0` / never | First / last step for insertion. |
| `rate_limit` | unlimited | Cap on total particles inserted. |

Rate mode still requires `material`, `radius`, and `density`. Its candidate
stream is seeded per-step-per-entry (`config_seed ^ step_hash ^ entry_hash`) so
candidates are deterministic across ranks. One caveat: the overlap scratch for
rate insertion covers only the *new* atoms of that step, not the existing local
population — new atoms can be born overlapping existing particles, so place
rate-insert regions in free space.

## Radius distributions

`radius` accepts a fixed value or one of four statistical distributions, to build
a polydisperse packing:

| `radius` value | Distribution |
|---|---|
| `0.001` | **Fixed** — every particle has this radius. |
| `{ distribution = "uniform", min = …, max = … }` | **Uniform** over `[min, max]`. |
| `{ distribution = "gaussian", mean = …, std = … }` | **Gaussian**, clamped to ≥ 1e-15 to avoid negative radii. |
| `{ distribution = "lognormal", mean = …, std = … }` | **Lognormal**. |
| `{ distribution = "discrete", values = […], weights = […] }` | **Discrete** — pick from `values` with the given `weights`. |

> **Lognormal `mean`/`std` describe the radius itself.** They are the desired
> mean and standard deviation *of the radius distribution*, not the μ and σ of
> the underlying normal. The code converts internally
> (`σ² = ln(1 + (std/mean)²)`, `μ = ln(mean) − σ²/2`), so you specify the
> statistics you actually want to see in the packing.

## Region types

`region` selects where particles may be placed (default: the domain inset by the
maximum radius). It accepts any `soil_core::Region`; the two most common:

```toml
region = { type = "block", min = [x0, y0, z0], max = [x1, y1, z1] }
region = { type = "cylinder", center = [cx, cy], radius = r, axis = "z", lo = z0, hi = z1 }
```

Region walls (`[[wall]]`) accept the same shapes plus cones, unions, and
intersections — see [Walls](../physics/walls.md).

## File-based insertion

For `source = "file"`, the `format` is `"csv"`, `"lammps_dump"`, or
`"lammps_data"`:

- **`csv`** — map columns with `columns = { x = 0, y = 1, z = 2, radius = 3 }`
  (the default). `radius` and `density` from the block fill in any column the
  file omits.
- **`lammps_dump`** — a LAMMPS dump snapshot; positions (and per-atom radius,
  where present) come from the file.
- **`lammps_data`** — a LAMMPS data file. The atom style is **auto-detected** from
  the `Atoms # <style>` comment: `sphere` and `bpm/sphere` encode
  `diameter density x y z` (per-atom radius and density read from the file), while
  `atomic` needs `radius` and `density` from the config. Override detection with
  an explicit `atom_style`.

When a file carries integer atom types, `type_map = { 1 = "glass", 2 = "steel" }`
maps each type to a named material; atoms whose type is absent from the map fall
back to the block's `material`.

```toml
[[particles.insert]]
source = "file"
file = "examples/bond_cantilever/chain.csv"
format = "csv"
material = "bpm"
radius = 0.001
density = 2500.0
columns = { x = 0, y = 1, z = 2 }
```

## A note on stable timesteps

The Rayleigh-wave criterion sets a stable DEM timestep: `dt_R = π r / α ·
√(ρ/G)` per particle (`α ≈ 0.1631 ν + 0.8766`, `G = E / (2(1+ν))`), with the
chosen step a fraction of the minimum across all particles. Stiff materials and
small particles force a tiny `dt` — typically `~1e-7 s` for the example configs.

## See also

- [Materials & the MaterialTable](materials.md) — how an inserted particle's
  `material` resolves into contact parameters.
- [Anatomy of a Config File](../getting-started/config-anatomy.md) — the
  `[[particles.insert]]` block in context.
