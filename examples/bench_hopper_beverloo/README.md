# bench_hopper_beverloo — Beverloo discharge law (2D slot)

Validates the steady mass-flow rate of a hopper draining under gravity against
**Beverloo's correlation**. The geometry is a **quasi-2D slot** hopper: the
domain is periodic in `y` over a depth of a few grain diameters, so the orifice
is a long slot of opening width `D` (flow is per unit slot depth). The matching
Beverloo exponent for a 2D slot is **3/2** (the 3D circular-orifice exponent 5/2
does **not** apply here — see *Assumptions*).

## Physics

Granular discharge from an orifice is governed by a free-fall arch above the
opening, not by the bed height (no Torricelli `√h` dependence). Beverloo, Leniger
& van de Velde (1961) found the empirical law

- 3D circular orifice:  `W = C · ρ_b · √g · (D − k·d)^(5/2)`
- **2D slot (this example, per unit width):**  `W = C · ρ_b · √g · (D − k·d)^(3/2)`

with bulk density `ρ_b`, grain diameter `d`, `g` gravity, `k ≈ 1.4` (the "empty
annulus" correction so grains cannot pass within `k·d/2` of the edge), and `C`
an order-1 constant. We validate the **exponent** and the **`(D − k·d)` form**,
not the exact `C`.

The recorder (`main.rs`) fills the hopper, removes a blocker wall when the
flowing stage begins, then logs the **cumulative discharged mass** (grains whose
center has dropped below the orifice plane) vs time. The steady slope of that
curve is `W`; `sweep.py` fits it for each `D` and regresses `W` vs `(D − k·d)` on
log–log.

## Material Properties

| Property | Value | Notes |
|---|---|---|
| Grain diameter `d` | 4.0 mm (`r` = 2 mm) | monodisperse |
| Density | 2500 kg/m³ | glass-like |
| Young's modulus | 5×10⁷ Pa | softened — keeps the rigid-grain Beverloo regime while allowing `dt = 2×10⁻⁵ s` |
| Poisson ratio | 0.3 | |
| Restitution | 0.5 | |
| Friction `μ` | 0.5 | needed for a discharge arch |
| Gravity `g` | 9.81 m/s² | |
| Particles | 1400 | full bed, drains in ≈ 0.3–1.1 s |
| Slot depth (periodic `y`) | 12 mm = 3·d | quasi-2D |
| `k·d` | 1.4 · 4 mm = 5.6 mm | Beverloo edge correction |

The funnel is a symmetric wedge converging from the bin walls (`x` = 20→140 mm)
at `z` = 180 mm down to a central slot at `z` = 50 mm.

## Parameter Sweep

The slot opening `D` is swept over five values, all `> 2.5·d` above the Beverloo
cutoff so flow never jams:

| `D` (mm) | 16 | 20 | 24 | 28 | 32 |
|---|---|---|---|---|---|
| `D − k·d` (mm) | 10.4 | 14.4 | 18.4 | 22.4 | 26.4 |

`config.toml` is the single representative case (`D` = 24 mm); `sweep.py` templates
the same layout for every `D`.

## Validation Criteria

`sweep.py graph` fits `ln W = m·ln(D − k·d) + b` and **PASSES** when:

- fitted exponent `m` is within **±0.25** of **3/2**,
- log–log fit quality `R² ≥ 0.97`,
- `W` increases monotonically with `D` (so `W → 0` as `D → k·d`).

LAMMPS is **not** used for this benchmark — a multi-stage hopper with a
runtime-removed blocker is not a drop-in LAMMPS input — so validation is against
Beverloo theory only.

## How to Run

```bash
python3 examples/bench_hopper_beverloo/sweep.py generate   # write per-case configs
python3 examples/bench_hopper_beverloo/sweep.py start      # build + run all D cases -> CSV
python3 examples/bench_hopper_beverloo/sweep.py graph      # validate (PASS/FAIL) + plot
python3 examples/bench_hopper_beverloo/sweep.py            # all three, in order
```

A single standalone discharge can be run directly:

```bash
cargo run --release --example bench_hopper_beverloo --no-default-features -- \
    examples/bench_hopper_beverloo/config.toml
```

### Outputs

| path | contents | tracked |
|---|---|---|
| `sweep/<case>/config.toml` | per-`D` DIRT configs | no (gitignored) |
| `data/sweep.csv` | fitted `W` vs `D` | no |
| `data/curve_D<...>.csv` | per-`D` cumulative-discharge curves | no |
| `plots/beverloo_W_vs_D.png` | `W` vs `(D − k·d)` log–log with fitted power law + 3/2 reference | **yes** |
| `plots/discharge_curves.png` | cumulative discharged mass vs time, one curve per `D` | **yes** |

## Expected Plots

- **`beverloo_W_vs_D.png`** — the five DIRT points fall on a straight log–log line
  bracketing the Beverloo 3/2 reference slope; the title reports the fitted
  exponent and `R²`.
- **`discharge_curves.png`** — five cumulative-mass curves, each with a clean
  constant-slope steady region (that slope is `W`), steeper for larger `D`, all
  plateauing at the same total bed mass (≈ 0.117 kg).

## Status / findings

**Validated.** The discharge is steady and the exponent matches the 2D-slot
Beverloo form:

- fitted exponent ≈ **1.36** (target 3/2; well within tolerance), `R²` ≈ **1.00**,
- `W` rises monotonically with `D` and the curve points toward `W → 0` near
  `D ≈ k·d`, confirming the `(D − k·d)^{3/2}` form (not a bare `D^{3/2}`).

The fitted exponent sits slightly below 3/2 — expected for a finite hopper with a
modest range of `(D − k·d)` and a converging wedge feed; the law is exact only in
the large-silo / small-orifice asymptote. The fit is essentially a perfect power
law (`R²` ≈ 1.0).

## Assumptions

- **2D slot, not 3D circular.** DIRT's plane/cylinder/cone wall primitives make a
  clean *slot* orifice (two inclined plane walls leaving a gap) straightforward,
  but they cannot cleanly cut a *circular hole* in a flat floor (a region wall
  built from a cone closes its narrow end with a disk cap, sealing the orifice).
  So this benchmark uses the slot geometry with the matching **3/2** exponent. A
  3D-circular variant would need a disk-with-hole / annulus wall primitive.
- **Rigid-grain regime.** The softened `E` keeps overlaps small relative to `d`
  while permitting a practical timestep; Beverloo is a rigid-grain law, so this is
  consistent.
- **Free-fall arch, not bed-height driven.** Validity assumes `D` ≫ `d` (no
  jamming) and a bed tall enough that `W` is independent of fill height during the
  measured window — both hold for the swept range.
- **Steady-state window.** `W` is fit over the 10–90 % portion of the discharge
  (excludes the brief start-up transient and the final empty-out tail).

## References

- W. A. Beverloo, H. A. Leniger, J. van de Velde, "The flow of granular solids
  through orifices", *Chem. Eng. Sci.* **15** (1961) 260–269.
- R. M. Nedderman, *Statics and Kinematics of Granular Materials*, Cambridge
  Univ. Press (1992), ch. 10 (slot vs circular Beverloo exponents).

## License

MIT OR Apache-2.0
