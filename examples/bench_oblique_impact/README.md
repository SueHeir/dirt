# bench_oblique_impact — tangential contact validation (Maw 1976)

Validates the **tangential** contact model (Hertz–Mindlin incremental spring +
Coulomb cap) against the oblique-impact theory of **Maw, Barber & Fawcett
(1976)** and the experiments of **Kharaz, Gorham & Salman (2001)**, and — when a
LAMMPS binary is available — against LAMMPS's equivalent `granular` Hertz–Mindlin
model. The normal model is validated separately by
[`bench_hertz_rebound`](../bench_hertz_rebound).

## Why a frozen sphere partner

This benchmark uses a frozen **sphere** as the contact partner to exercise the
**particle–particle** tangential model directly, and to mirror the LAMMPS
sphere–sphere setup for the cross-check. The target sphere is **frozen**
(`[[freeze]]`) — an immovable, non-spinning, fully frictional partner — and a
projectile strikes it obliquely. Both DIRT and LAMMPS aim the projectile
dead-center so the impact normal is exactly +z; DIRT additionally decomposes
velocities in the line-of-centers frame, so any residual tilt is handled exactly.
(`dirt_wall` now also carries a Mindlin tangential spring with a Coulomb cap, so a
flat wall could serve as the partner too; the sphere keeps this test focused on the
particle–particle path.)

## Reference

With contact-point tangential surface velocity `v_s` (incident, `= v_t` since
spin is initially zero) and `v_s' = v_t' − R·ω'` (rebound):

- non-dimensional incidence angle  `ψ₁ = [2(1−ν)/(μ(2−ν))] · (v_s/v_n)`
- tangential restitution            `β = −v_s'/v_s`
- exact gross-sliding (rigid) limit `β = −1 + 7(1+eₙ)(1−ν)/[(2−ν)·ψ₁]`

Maw's full solution interpolates a **sticking plateau** (small ψ₁, β → −1 at
near-normal incidence), a **microslip transition** through a positive peak, and
the **gross-slip** branch (large ψ₁) — the textbook S-curve.

## Running

Everything is driven by a single script, `sweep.py`:

```bash
python3 examples/bench_oblique_impact/sweep.py generate   # write per-case configs
python3 examples/bench_oblique_impact/sweep.py start      # build + run all sims -> CSV
python3 examples/bench_oblique_impact/sweep.py graph      # validate + plot
python3 examples/bench_oblique_impact/sweep.py            # all three, in order
```

If `lmp_serial` / `lmp` / `lmp_mpi` / `lammps` is on `PATH`, each case is also run
in LAMMPS and overlaid; otherwise only DIRT runs. A single standalone impact can
be run directly against `config.toml`:

```bash
cargo run --release --example bench_oblique_impact --no-default-features -- \
    examples/bench_oblique_impact/config.toml
```

### Outputs

| path | contents | tracked |
|---|---|---|
| `sweep/<case>/` | per-case DIRT configs + LAMMPS inputs | no (gitignored) |
| `data/sweep.csv`, `data/sweep_lammps.csv` | β(ψ₁) results | no |
| `data/trace_dirt.csv`, `data/trace_lammps.csv` | per-step contact-force trace | no |
| `plots/beta_vs_psi1.png` | β vs ψ₁: DIRT vs LAMMPS vs gross-slip theory | **yes** |
| `plots/contact_trace.png` | per-step normal & tangential force loops (ψ₁≈1.7) | **yes** |

## Status / findings

**Validated.** DIRT reproduces the full Maw S-curve and matches LAMMPS across the
whole range (max `|Δβ| ≈ 0.007`):

- **Normal restitution** is constant at `eₙ ≈ 0.985` independent of tangential
  velocity (spread < 0.002), confirming the normal model is decoupled from the
  tangential one.
- **Tangential restitution** follows Maw: the `β ≈ −1` sticking plateau at low
  ψ₁, the microslip rise through a `+0.32` peak near ψ₁≈3.3, and convergence onto
  the analytical gross-slip branch at high ψ₁.
- The per-step trace shows DIRT and LAMMPS tracing an identical normal curve and
  tangential loading/unloading hysteresis loop.

This validation also drove two fixes in the contact model: a tangential
damping-sign error (energy injection) and the requirement that a frozen contact
partner have its rotation frozen too (`[[freeze]]`, not a translation-only pin).

## License

MIT OR Apache-2.0
