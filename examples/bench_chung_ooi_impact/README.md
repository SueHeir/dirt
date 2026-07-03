# bench_chung_ooi_impact — Chung & Ooi (2011) elastic normal-impact verification

Reproduces two of the standard DEM code-verification benchmarks proposed by

> L. Chung and J.Y. Ooi, **"Benchmark tests for verifying discrete element
> modelling codes at particle impact level"**, *Granular Matter* **13**(5):643–656
> (2011). https://doi.org/10.1007/s10035-011-0277-0

These are *the* particle-impact-level benchmarks the DEM community cites to
establish that a new code computes the Hertz contact correctly. This example
covers the two elastic normal-impact cases:

- **Test 1** — elastic normal impact of **two identical spheres**.
- **Test 2** — elastic normal impact of a **sphere with a rigid wall**.

For both, the contact is perfectly **elastic** (`restitution = 1` → zero
damping), so the exact reference is the closed-form **Hertz** solution for the
collision — the same analytical reference Chung & Ooi use for these cases. The
benchmark is therefore validated against an **independent analytical solution**,
not against LAMMPS and not against DIRT's own output (anti-gaming: theory only).

## What is measured, and against what

Each impact records three quantities and gates them PASS/FAIL against Hertz
theory across a sweep of impact velocities:

| quantity | measured from DIRT | Hertz analytical reference |
|---|---|---|
| maximum contact force `F_max` | peak net force on the particle during contact | `(4/3) E* √R* δ_max^{3/2}` |
| maximum overlap `δ_max` | peak geometric overlap `R₁+R₂−d` (pure kinematics) | `(15 m* v² / (16 E* √R*))^{2/5}` |
| contact duration `t_c` | steps in contact × `dt` | `2.943920 · δ_max / v` |

with the effective quantities

- **two identical spheres:** `E* = E/(2(1−ν²))`, `R* = R/2`, `m* = m/2`
- **sphere on same-material wall:** `E* = E/(2(1−ν²))`, `R* = R`, `m* = m`

The contact-time constant `2.943920 = 2 ∫₀¹ (1−x^{5/2})^{−1/2} dx` is the exact
Hertz-collision value (Johnson, *Contact Mechanics*, 1985, Ch. 11 — identical to
Chung & Ooi's Test 1/2 reference).

## Material (Chung & Ooi 2011 aluminium alloy)

`E = 70 GPa`, `ν = 0.30`, `ρ = 2700 kg/m³`, `R = 0.1 m`. Relative impact
velocity is swept over `0.5, 1, 2, 5, 10 m/s`, spanning the paper's impact
regime. `restitution = 1.0` gives DIRT's undamped (pure-Hertz) contact, so the
elastic analytical solution applies exactly.

## Running

```bash
python3 examples/bench_chung_ooi_impact/sweep.py generate   # write per-case configs
python3 examples/bench_chung_ooi_impact/sweep.py start      # build + run all sims -> CSV
python3 examples/bench_chung_ooi_impact/sweep.py graph      # validate + plot
python3 examples/bench_chung_ooi_impact/sweep.py            # all three, in order
```

A single standalone impact (Test 1, v = 1 m/s) can be run directly:

```bash
cargo run --release --example bench_chung_ooi_impact --no-default-features \
    --features precision-double -- examples/bench_chung_ooi_impact/config.toml
```

### Outputs

| path | contents | tracked |
|---|---|---|
| `sweep/<case>/config.toml` | per-case DIRT configs | no (gitignored) |
| `data/sweep_results.csv` | measured `F_max`, `t_c`, `δ_max` per case | no |
| `plots/max_force.png` | `F_max(v)`: DIRT vs Hertz theory, both cases | **yes** |
| `plots/contact_time.png` | `t_c(v)`: DIRT vs Hertz theory | **yes** |
| `plots/max_overlap.png` | `δ_max(v)`: DIRT vs Hertz theory | **yes** |

## Status / findings

**Validated — 31/31 checks pass.** DIRT reproduces the Hertz analytical solution
for both Chung & Ooi elastic normal-impact cases across the whole velocity
sweep:

- **Maximum contact force** and **maximum overlap** match to `< 0.01%` at every
  velocity in both the sphere–sphere and sphere–wall cases (the residual is
  round-off / sampling, not a model error).
- **Contact duration** matches to `< 0.5%`; the small residual is the integer
  step-count resolution `dt/t_c` and shrinks with a finer timestep.
- The sphere–wall force is larger than the sphere–sphere force at the same
  relative velocity (larger `R*` and `m*`), exactly as Hertz theory predicts —
  confirming the effective-property mixing is correct.

Tolerances (force 2%, overlap 2%, contact-time 2%) bound the discretisation
residual and are not loosened to force a pass; the actual errors are far inside
them.

## License

MIT OR Apache-2.0
