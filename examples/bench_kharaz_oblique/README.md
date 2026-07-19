# bench_kharaz_oblique — replicate Kharaz, Gorham & Salman (2001)

Reproduces the oblique-impact experiment of **Kharaz, Gorham & Salman, "An
experimental study of the elastic rebound of spheres", Powder Technology 120
(2001) 281-291**: a 5 mm **aluminium-oxide** sphere striking a thick soda-lime
**glass** anvil (a fully elastic response) at a **fixed impact speed
Vᵢ = 3.85 m/s**, with the **angle of incidence swept** from near the surface
normal to near grazing. This is the classic elastic oblique-impact benchmark and
is complementary to [`bench_oblique_impact`](../bench_oblique_impact) (which
sweeps at fixed normal velocity and validates β(ψ₁) against Maw theory + LAMMPS):
here we reproduce **Kharaz's actual experimental protocol and figure format** —
the rebound/spin curves plotted against incidence angle.

## Geometry — a flat glass anvil (not a frozen sphere)

Kharaz used a thick glass block, i.e. a flat half-space. This example fires the
sphere at a real **`dirt_wall` z-plane** (normal +z) carrying the Mindlin
tangential spring + Coulomb cap. A flat wall keeps the contact normal exactly +z
**throughout** the collision at every incidence angle, so the normal restitution
stays angle-independent right up to grazing. (A frozen *sphere* partner — as used
by `bench_oblique_impact` to exercise the particle–particle path — curves the
projectile at grazing incidence and depresses the apparent eₙ; that is why this
replication uses a wall.) The wall is infinite, so no aiming is needed: the
projectile is launched with velocity `(v_t, 0, −v_n)`, `v_t = Vᵢ sinΘᵢ`,
`v_n = Vᵢ cosΘᵢ`.

## Conditions (Kharaz 2001, glass anvil)

| quantity | value |
|---|---|
| impact speed Vᵢ | 3.85 m/s (fixed) |
| incidence angle Θᵢ (from normal) | 5°–80° sweep |
| sphere | 5 mm alumina, ρ = 4000 kg/m³, E = 380 GPa, ν = 0.23 |
| normal restitution eₙ | **0.98** (measured; ~angle-independent) |
| sliding friction μ | **0.092** (measured) |

## Reference — independent, not self-consistent

In the **gross-sliding** regime the contact point slides throughout the collision
and the rebound is fixed by the exact rigid-body impulse relations (no free
parameters beyond Kharaz's measured μ and eₙ):

```
v_n' = eₙ v_n
v_t' = v_t − μ(1+eₙ) v_n          ⇒  e_t = 1 − μ(1+eₙ)/tan Θᵢ
ω'   = 5 μ(1+eₙ) v_n / (2R)        ⇒  Rω'/Vᵢ = 5 μ(1+eₙ) v_n /(2 Vᵢ)
tan Θ_r = v_t'/v_n'
```

valid while `tan Θᵢ > (7/2) μ(1+eₙ)` — here Θᵢ ≳ **32.5°**. This is exactly the
sliding-branch behaviour Kharaz's glass-anvil data confirmed. Below the boundary
the contact sticks / micro-slips (Maw's S-curve); there DIRT is checked
qualitatively (spin stays at or below the sliding line).

## Running

```bash
python3 examples/bench_kharaz_oblique/sweep.py generate   # per-angle configs
python3 examples/bench_kharaz_oblique/sweep.py start       # build + run -> CSV
python3 examples/bench_kharaz_oblique/sweep.py graph       # validate + plot
python3 examples/bench_kharaz_oblique/sweep.py             # all three
```

A single standalone impact (Θᵢ = 45°) can be run against `config.toml`:

```bash
cargo run --release --example bench_kharaz_oblique --no-default-features \
    --features precision-double -- examples/bench_kharaz_oblique/config.toml
```

### Outputs

| path | contents | tracked |
|---|---|---|
| `sweep/<case>/` | per-angle DIRT configs | no (gitignored) |
| `data/kharaz_sweep.csv` | reduced rebound/spin quantities | no |
| `data/kharaz_lammps.csv` | matched LAMMPS rebound/spin quantities | no |
| `plots/kharaz_rebound_spin.png` | the four Kharaz curves vs incidence angle | **yes** |
| `kharaz_experiment.csv` | *optional* digitised experimental points (overlaid if present) | yes |

## Status / findings

**PASS.** DIRT reproduces the Kharaz rebound/spin curves at the paper's exact
conditions and overlays the Maw/Hertz–Mindlin analytical micro-slip curves:

- **Normal restitution** is flat at **eₙ = 0.986** across the whole 5°–80° sweep
  (spread 0.0000), within 0.006 of Kharaz's measured glass value 0.98.
- **Sliding regime (Θᵢ ≳ 32.5°):** the rebound angle Θ_r, tangential restitution
  eₜ = v_t'/v_t, and non-dimensional rebound spin Rω'/Vᵢ match the exact
  rigid-body kinematics to three decimals (max |Δeₜ|, |ΔRω/Vᵢ| < 0.001,
  |ΔΘ_r| < 0.1°).
- **Sticking / micro-slip regime (Θᵢ < 32.5°):** DIRT traces the Maw S-curve —
  tangential restitution dips to a minimum eₜ ≈ 0.62 near Θᵢ ≈ 20°, and the
  rebound spin peaks at Rω'/Vᵢ ≈ 0.39 near Θᵢ ≈ 30° before both join the sliding
  branch. The contact-point restitution agrees with the analytical curve to
  max `|Δβ| = 0.0015`.
- **Cross-code:** all 16 matched LAMMPS wall impacts overlay DIRT. Maximum
  DIRT–LAMMPS differences are `0.00074` in eₜ, `0.00045` in Rω'/Vᵢ,
  `0.0105°` in rebound angle, and `0.00260` in contact-point β.

![Kharaz rebound/spin curves](plots/kharaz_rebound_spin.png)

*Rebound angle, tangential restitution, non-dimensional rebound spin, and normal
restitution vs incidence angle. Filled circles: DIRT; open squares: matched
LAMMPS wall impacts. Green: analytical Maw/Hertz–Mindlin curves. Dotted: exact
rigid-body sliding kinematics where valid. Dashed red: Kharaz's measured
eₙ = 0.98. The bottom row shows the contact-point β curve and the direct
DIRT–LAMMPS β residual; validation limits remain in the harness rather than as
shaded figure bands.*

### A note on the experimental points

Kharaz's paper reports its **glass-anvil per-point scatter only in figures** in
the (paywalled) publication; there is no open-access copy or supplementary data
file. This benchmark therefore validates against the exact rigid-body sliding
kinematics and Maw–Barber–Fawcett theory that those experimental points confirmed
(Kharaz's central result was that the elastic-rebound data follow this theory to
within experimental scatter), anchored to Kharaz's **measured** scalars
eₙ = 0.98 and μ = 0.092. If the digitised experimental points become available,
place them in `kharaz_experiment.csv` (columns any of `theta_deg,e_t,spin_nd,
theta_r_deg`) and `sweep.py graph` will overlay them automatically.

## License

MIT OR Apache-2.0
