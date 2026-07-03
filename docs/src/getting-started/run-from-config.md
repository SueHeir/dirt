# Run a Simulation from a Config — Zero Rust

You do **not** need to write, edit, or compile any Rust to run a DIRT
simulation. DIRT ships a single generic driver — the `run` example — that
assembles the standard plugin stack and then takes *every* case-specific
detail (geometry, materials, particle insertion, walls, body forces, loading,
output, duration) from a **TOML config file** you pass on the command line.

The workflow is:

> **Describe** your simulation in a config, then **run the binary.**

You *describe* a scenario declaratively — what materials, what geometry, what
particles, how long. You do **not** *script* it: there is no user-authored
sequence of operations, no control flow, no per-case code. The driver owns the
step loop; your config only states *what* to simulate. (This is a deliberate
design choice — see [the note on declarative-vs-scripted](#declarative-not-scripted)
at the end.)

## 1. Install

Follow [Installation & Building](./installation.md) to get a serial build. In
short:

```bash
git clone https://github.com/SueHeir/dirt
cd dirt
cargo build --release --no-default-features
```

`--no-default-features` skips the MPI backend for the quickest single-process
build. That's all the setup you need — no code to touch.

## 2. Pick a config

The driver ships several complete, self-contained scenarios under
`examples/run/`. Each is a full simulation expressed **entirely** in TOML — pick
one and pass it to the driver:

| Config                       | Scenario                    | Geometry                | Insertion          |
|------------------------------|-----------------------------|-------------------------|--------------------|
| `config.toml`                | two-sphere settle           | floor plane             | 2 fixed spheres    |
| `pour_settle.toml`           | granular pour & settle      | box (5 plane walls)     | polydisperse cloud |
| `pour_cylinder.toml`         | bidisperse silo pour        | cylinder wall + floor   | two materials      |
| `shear_box.toml`             | Lees–Edwards simple shear   | triperiodic cell        | dense glass pack   |
| `compression_box.toml`       | uniaxial compression        | triperiodic cell        | loose glass pack   |

The **same** driver runs all of them; the config alone decides which physics is
active. A plugin whose config section is absent simply does nothing, so one
binary covers a two-particle rebound, a settle-into-a-box, and a sheared cell.

## 3. Run it

Run the granular pour-and-settle case — a cloud of glass beads dropped into an
open-topped box under gravity:

```bash
cargo run --release --no-default-features --example run -- examples/run/pour_settle.toml
```

The single trailing argument is the config path. There is **nothing else to
edit** — no `main.rs`, no recompiling per case. You'll see the driver banner,
an echo of the parsed config, and a running thermo line:

```text
Comm: processors 1 1 1
Domain: 0 0.06 0 0.06 0 0.16
Domain: boundary f f f
Neighbor: skin_fraction=1.2 bin_size=0.011 rebuild=displacement newton=on
Run: 80000 steps
DemAtomInsert: inserting 200 particles of material 'glass' (r=0.005, rho=2500, E=8700000000, nu=0.3)

Step         Atoms        Ke           Neighbors    Walltime     Stepps
-----------------------------------------------------------------------------
0            200          4.293040e-11 131          0.0006       0.0
```

When it finishes you'll have standard dump output next to the config
(`dump/dump_{step}.csv` per-atom snapshots and `vtp/*.vtp` visualization
frames), ready for post-processing or viewing in ParaView.

> To run any other scenario, swap the config path — e.g.
> `... --example run -- examples/run/shear_box.toml`. Same binary, different
> description.

## 4. Read the config you just ran

Everything the driver did came from `pour_settle.toml`. Here are the pieces,
each a plain declarative section:

```toml
[comm]                     # MPI rank grid; 1×1×1 = serial
processors_x = 1
processors_y = 1
processors_z = 1

[domain]                   # the simulation box + boundary conditions
x_low = 0.0
x_high = 0.06
y_low = 0.0
y_high = 0.06
z_low = 0.0
z_high = 0.16
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[gravity]                  # body force — activates GravityPlugin
gz = -9.81

[[dem.materials]]          # named material, referenced by particles & walls
name = "glass"
youngs_mod = 8.7e9
poisson_ratio = 0.3
restitution = 0.3
friction = 0.5

[[particles.insert]]       # ~200 polydisperse beads, seeded non-overlapping
material = "glass"
count = 200
radius = { distribution = "uniform", min = 0.0035, max = 0.0050 }
density = 2500.0
region = { type = "block", min = [0.006, 0.006, 0.07], max = [0.054, 0.054, 0.15] }
seed = 7

[[wall]]                   # floor plane; four more [[wall]] blocks make the box
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "glass"
name = "floor"

[dump]                     # standard per-atom output
interval = 5000
format = "text"

[[run]]                    # duration: 80 000 steps of 5 µs = 0.4 s
name = "pour_settle"
dt = 5.0e-6
steps = 80000
thermo = 5000
```

To make your *own* run, copy this file, change the numbers — bigger box, more
particles, a stiffer material, a longer run — and pass your copy to the same
command. Still no Rust. The full field-by-field breakdown is in
[Anatomy of a Config File](./config-anatomy.md) and the complete schema is the
[Configuration Reference](../reference/config.md).

## Declarative, not scripted

DIRT gives scientists a general runner binary on purpose — but it deliberately
**rejects the LAMMPS-style input-script model**. The config *describes* a
simulation; it is never a scripting language. There is no control flow, no
user-authored operation sequencing, no "do step A then step B" schedule. The
runner owns the schedule.

Even a *loaded* case — deforming a cell up to a target strain — stays
declarative. The `shear_box.toml` and `compression_box.toml` configs describe
the loading with a single `[loading]` block:

```toml
[loading]
type          = "shear"    # "shear" (Lees–Edwards xy) | "compression" (uniaxial)
rate          = 50.0       # engineering strain rate (1/s)
target_strain = 2.0        # dimensionless total strain to reach
dt            = 2.0e-7     # timestep (s)
```

You state *what* deformation to apply and *how far*; the runner derives the
step count (`target_strain / (|rate|·dt)`) and drives the box there itself. You
never write the loop. If a scenario can't be expressed declaratively today, the
answer is a new **declarative option**, not a scripting hook.

That's the whole idea: **describe your simulation in a config and run the
binary — never script it.**
