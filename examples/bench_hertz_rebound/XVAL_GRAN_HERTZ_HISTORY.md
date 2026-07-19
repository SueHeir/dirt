# Independent LAMMPS cross-validation — `gran/hertz/history` vs DIRT

**Goal:** `dirt-lammps-xval-hertz` — validate DIRT's Hertz normal contact against
the **classic `pair_style gran/hertz/history`** implementation in the installed
LAMMPS binary, comparing **COR**, **contact time** (`t_c`), and **peak overlap**
(`δ_max`) within stated tolerances, for independent-code provenance.

Run it:

```bash
$BENCH_PYTHON examples/bench_hertz_rebound/sweep.py start   # produces DIRT reference (data/sweep_results.csv)
$BENCH_PYTHON examples/bench_hertz_rebound/xval_gran_hertz_history.py
```

Driver: [`xval_gran_hertz_history.py`](xval_gran_hertz_history.py). Generated
results are written to `data/xval_gran_hertz_history.csv`; per-case LAMMPS
inputs, logs, and traces are written under `xval/`.

## Why this is *independent* provenance

`sweep.py` already overlays a LAMMPS run, but it uses the **modern**
`pair_style granular` + `hertz/material` model (`src/GRANULAR/granular_model.cpp`,
`gsm_normal.cpp`). This cross-check instead drives the **classic**
`pair_style gran/hertz/history` / `fix wall/gran hertz/history` path
(`src/GRANULAR/pair_gran_hertz_history.cpp`, the pre-2019 granular kernel) — a
*separate implementation* with its own force assembly and a **different damping
law**. Agreement between DIRT and this path is therefore genuine cross-code
provenance for the Hertz normal contact, not self-consistency.

Installed binary: LAMMPS **22 Jul 2025 - Update 4** (`~/.local/bin/lmp_serial`,
`stable_22Jul2025_update4`).

## Physics mapping (derivation, not a fit)

The classic Hertzian style uses (LAMMPS `doc/src/pair_gran.rst`):

```
F_hz = √δ · √(Ri·Rj/(Ri+Rj)) · [ Kn·δ·n_ij − m_eff·γ_n·v_n − (tangential) ]
```

with the documented **material mapping** `Kn = 4G/(3(1−ν)) = 2E/(3(1−ν²))`,
`G = E/(2(1+ν))`. For a sphere on a flat wall (`Rj → ∞`) the geometric prefactor
is `√R`, so the elastic normal force is

```
F = Kn · √R · δ^{3/2} = (4/3)·E* · √R · δ^{3/2},    E* = E/(2(1−ν²)),
```

**identical** to true Hertz and to DIRT's spring (DIRT uses the same
`E* = E/(2(1−ν²))` convention). Hence `Kn = (4/3)·E*` — the script asserts this
equality — and `δ_max`, `t_c` are **stiffness+geometry predictions, not fits**.

Material: `E = 70 GPa`, `ν = 0.22`, `R = 5 mm`, `ρ = 2500 kg/m³` (same as the
bench). Computed `Kn = 4.904021×10¹⁰ Pa`. Frictionless (`xmu = 0`), no gravity,
`fix nve/sphere`, same timestep as DIRT.

## Two comparisons

### 1. Near-elastic anchor (all four impact speeds)

The classic parser sets the tangential coefficient to `gammat/gamman`
(`granular_model.cpp:53`), so **exactly-zero normal damping is a 0/0 → NaN** and
is not representable in `gran/hertz/history`. We therefore use the smallest
`γ_n` that realizes `COR ≥ 0.9999` (calibrated at the fastest speed, where a
constant-γ_n COR is lowest). Damping is then negligible and `δ_max`, `t_c` sit on
the undamped Hertz values.

| v₀ (m/s) | COR (ghh) | δ_max ghh | δ_max DIRT | δ_max Hertz | t_c ghh | t_c DIRT | t_c Hertz |
|---|---|---|---|---|---|---|---|
| 0.1 | 0.99993 | 1.8599e-6 | 1.8602e-6 | 1.8600e-6 | 5.4882e-5 | 5.4882e-5 | 5.4779e-5 |
| 0.5 | 0.99993 | 6.7382e-6 | 6.7393e-6 | 6.7406e-6 | 3.9637e-5 | 3.9637e-5 | 3.9703e-5 |
| 1.0 | 0.99989 | 1.1734e-5 | 1.1735e-5 | 1.1736e-5 | 3.5064e-5 | 3.5064e-5 | 3.4564e-5 |
| 2.0 | 0.99994 | 2.0428e-5 | 2.0427e-5 | 2.0434e-5 | 3.0490e-5 | 3.0490e-5 | 3.0089e-5 |

`gran/hertz/history` reproduces DIRT's `δ_max`/`t_c` to **< 0.02 %** and the
analytic Hertz closed form to **< 0.5 % (δ_max) / < 1.5 % (t_c, ≈1 timestep)**,
with COR = 1.000 to `< 5×10⁻⁴`.

### 2. Damped COR (e = 0.7, 0.9)

`gran/hertz/history` damps with a **constant `γ_n`** (older viscoelastic dashpot,
`√δ·v_n`), whereas DIRT / `granular` use the **Tsuji (1992) polynomial**. The two
damping *laws* differ, so their COR-vs-speed behaviour differs:

- constant-`γ_n` Hertz dashpot → COR **weakly velocity-dependent** (∝ v₀^{1/5});
- Tsuji → COR **velocity-independent**.

This is a documented, expected model distinction, not a code error (it is the
same reason DIRT itself moved off a linear damping ratio to Tsuji — see the bench
README). We calibrate the single free knob `γ_n` (via an independent RK4
integration of the contact ODE) so `gran/hertz/history` realizes **DIRT's
measured COR at the reference speed v₀ = 1.0 m/s**, confirm LAMMPS reproduces it,
and check the untuned `δ_max`/`t_c` still match.

At the reference speed **v₀ = 1.0 m/s** (calibration point):

| nominal e | γ_n | COR ghh | COR DIRT | δ_max ghh | δ_max DIRT | t_c ghh | t_c DIRT |
|---|---|---|---|---|---|---|---|
| 0.7 | 9.049e7 | 0.77480 | 0.77568 | 1.0617e-5 | 1.0631e-5 | 3.5826e-5 | 3.5826e-5 |
| 0.9 | 2.631e7 | 0.92760 | 0.92781 | 1.1386e-5 | 1.1389e-5 | 3.5064e-5 | 3.5064e-5 |

→ COR agrees to **≤ 9×10⁻⁴**, `δ_max`/`t_c` to **< 0.2 %**.

Across all four speeds the constant-`γ_n` COR spread is **0.104 (e=0.7)** and
**0.036 (e=0.9)**, versus DIRT/Tsuji ≈ 0 — the expected `v₀^{1/5}` model
difference. Full per-speed numbers are in the CSV; at off-reference speeds
`δ_max`/`t_c` diverge consistently with the COR divergence, so those are reported
but not asserted.

## Stated tolerances (all PASS)

| Check | Tolerance | Result |
|---|---|---|
| `Kn = (4/3)E*` mapping identity | exact (1e-12) | ✅ |
| Anchor COR ≈ 1 | \|COR−1\| ≤ 0.005 | ✅ (≤ 5e-4) |
| Anchor δ_max vs DIRT | ≤ 1 % | ✅ (< 0.02 %) |
| Anchor δ_max vs Hertz theory | ≤ 2 % | ✅ (< 0.5 %) |
| Anchor t_c vs DIRT | ≤ 3 % | ✅ (< 0.02 %) |
| Anchor t_c vs Hertz theory | ≤ 2 % | ✅ (< 1.5 %) |
| Damped COR vs DIRT @ v₀=1 | ≤ 0.01 | ✅ (≤ 9e-4) |
| Damped δ_max vs DIRT @ v₀=1 | ≤ 1 % | ✅ (< 0.2 %) |
| Damped t_c vs DIRT @ v₀=1 | ≤ 3 % | ✅ (< 0.1 %) |

No DIRT reference value or bench tolerance was altered; this is an additive,
independent overlay. `Kn` is derived from `E`, `ν` only (not fitted); the single
calibrated quantity is the damped `γ_n`, and the independent predictions that
validate the contact are the untuned `δ_max` and `t_c`.

## References

- LAMMPS `doc/src/pair_gran.rst` — `gran/hertz/history` force law + `Kn = 4G/(3(1−ν))` mapping.
- LAMMPS `doc/src/fix_wall_gran.rst` — `fix wall/gran hertz/history` wall contact.
- LAMMPS `src/GRANULAR/pair_gran_hertz_history.cpp`, `granular_model.cpp` (classic path).
- K.L. Johnson, *Contact Mechanics*, CUP 1985 (Hertz `δ_max`, `t_c` closed forms).
