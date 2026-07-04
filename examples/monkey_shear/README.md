# monkey_shear — LEBC volume-fraction shear campaign

Constant-volume Lees–Edwards simple-shear rheology of three particle shapes at a
**common equivalent-volume diameter** `D_eq = 0.1 m`, so each shape displaces the
identical solid volume `V_eq = π·D_eq³/6 = 5.236e-4 m³` and shares one `N(Φ)` table:

| type   | model                                             |
|--------|---------------------------------------------------|
| sphere | single sphere, `radius = 0.05`                    |
| rigid  | 44-sub-sphere "monkey" as a rigid `[[clump.insert]]` body (random orientation) |
| bpm    | same monkey as free sub-spheres welded intra-monkey by `auto_bond` |

Final box `2×2×1` (vol 4), fully periodic; flow=x, gradient=y, vorticity=z. Material =
reference glass (`E 7.0e7, ν 0.245, e 0.926, μ 0.16, ρ 2500`). Per (type, Φ) 3-stage
protocol: loose insert → settle → quasi-static isotropic `erate` compress to the target
Φ → constant-volume Lees–Edwards `xy erate` shear (γ̇ = 1.0 s⁻¹, strain ≈ 2).

## Files
- `main.rs` — the one binary (all shapes, chosen by config). Recorder reused verbatim
  from `bench_lebc_shear` (p, τ, N₁, N₂, granular T, Φ → `<output>/lebc_shear_results.csv`).
- `monkey_Deq0.1.toml` — the scaled monkey clump definition (Monte-Carlo D_eq-verified).
- `tools/gen_series.py` — generates the whole grid.
- `configs/<type>/phi_<Φ>.toml` — the 3×8 config grid.

## Regenerate the grid
```bash
python3 examples/monkey_shear/tools/gen_series.py all      # build-monkey + full grid
python3 examples/monkey_shear/tools/gen_series.py series --smoke   # Φ=0.1 short smoke
```

**bpm CSVs are generated, not committed.** Each `configs/bpm/phi_<Φ>.toml` reads a
companion `configs/bpm/phi_<Φ>.csv` (all sub-sphere centres + radii, placed with a
1.25·(Rᵢ+Rⱼ) inter-monkey gap so `auto_bond` welds only intra-monkey). Those CSVs total
~42 MB and are excluded by the repo's `*.csv` ignore rule; the placement is deterministic
(seeded), so `gen_series.py series` reproduces them exactly. **Run the generator before
launching any bpm config.**

## Run one case
```bash
source ~/projects/.build-env
cargo run --release --example monkey_shear --no-default-features --features precision-double \
    -- examples/monkey_shear/configs/bpm/phi_0.1.toml
```

## Smoke validation (Φ=0.1, this PR)
All three compile, run to completion, finite p/τ, no NaN, no overlap blow-up. bpm forms
`N×97 = 74108` bonds (no explosion). Monkeys tangle → p ≫ sphere p at equal Φ, so they are
near the dense regime already at Φ=0.1; the mid-series (Φ=0.3–0.4) inertial number falls to
the target `I ≈ 0.05–0.1`.
