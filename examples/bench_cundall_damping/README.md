# bench_cundall_damping — Cundall non-viscous global damping

Validates DIRT's `[[cundall]]` fix (`dirt_fixes`) — the non-viscous global
damping of **Cundall (1987)**, matching LAMMPS
[`fix damping/cundall`](https://docs.lammps.org/fix_damping_cundall.html) — against
its **exact analytical** effect, and cross-checks the linear part against the
LAMMPS binary running the native fix.

## The fix

Unlike velocity-proportional viscous damping (`[[viscous]]`, `F = −γv`), Cundall
damping scales a dimensionless *fraction* of the current force/torque,
component-by-component, keyed to the sign of the mechanical power:

```
F_k <- F_k · (1 − γ_l · sign(F_k · v_k))
T_k <- T_k · (1 − γ_a · sign(T_k · ω_k))
```

with `sign(x) = +1` for `x ≥ 0`, else `−1` (the exact LAMMPS convention). A
component's force is *reduced* when it acts along the motion and *amplified* when
it opposes it, so mechanical power is always non-positive and kinetic energy
monotonically dissipates. The coefficients are dimensionless and
velocity-independent — the standard tool for quasi-static DEM settling / energy
minimization. It is selectable in config alongside the existing viscous damping.

## The test

A single sphere in an empty box (no walls, no contacts) is:

* **launched straight up (+z) under gravity** — the linear driver, and
* **spun about +z with a constant applied torque** (opposing the spin; applied in
  `main.rs`, since there is no "angular gravity" fix) — the angular driver.

Because the force sign flips as the sphere passes its **apex** (rising → falling)
and the torque sign flips as the spin passes **zero**, each motion is
piecewise-constant-acceleration and every branch of the Cundall sign function is
exercised. The exact, mass-independent rates are

| phase | condition | rate |
|-------|-----------|------|
| rising    | `v_z > 0` (gravity opposes) | `a_up   = −g (1 + γ_l)` |
| falling   | `v_z < 0` (gravity along)   | `a_down = −g (1 − γ_l)` |
| spin-down | `ω_z > 0` (torque opposes)  | `α_down = T_z (1 + γ_a) / I` |
| spin-up   | `ω_z < 0` (torque along)    | `α_up   = T_z (1 − γ_a) / I` |

`main.rs` writes a self-describing CSV whose header carries the *actual* run
parameters (`g`, `I⁻¹`, applied torque, and `γ_l`/`γ_a` read back from the live
fix), so `sweep.py` validates against theory with zero duplicated constants.

## Running

```bash
# one demonstration case
cargo run --release --example bench_cundall_damping --no-default-features \
  --features precision-double -- examples/bench_cundall_damping/config.toml

# full PASS/FAIL sweep over gamma in {0.2, 0.5, 0.8}
python3 examples/bench_cundall_damping/sweep.py
```

`sweep.py` fits all four phase rates for every γ and requires each within **1 %**
of theory (the rates are exact, so the residual is just the single sign-flip
transition step). **PASS** ⇔ all four rates match for all γ.

![Cundall fitted rates](plots/cundall_rates.png)

*Fitted linear and angular rates vs the exact Cundall analytical rates across γ.
DIRT matches all four rates within the 1 % gate; the linear branch also overlays
the LAMMPS `fix damping/cundall` cross-check. Latest run: PASS.*

![Cundall damping traces](plots/cundall_traces.png)

*Velocity and angular-velocity traces for the swept γ values, showing the
piecewise-linear branches and the sign flip at the apex / zero spin. Latest run:
PASS.*

## Independent cross-check (LAMMPS)

If a LAMMPS binary is on `PATH`, the **linear** part of each case is also run in
LAMMPS with the native `fix damping/cundall` (single sphere, `fix gravity`
declared *before* the damping fix so gravity is damped) and its fitted
`a_up`/`a_down` must match theory too. This is genuine cross-code provenance, not
self-consistency. Latest run (`γ_l = γ_a = γ`):

```
 gamma   DIRT a_up    LMP a_up   DIRT a_dn    LMP a_dn
  0.20     -11.772     -11.772      -7.848      -7.848
  0.50     -14.715     -14.715      -4.905      -4.905
  0.80     -17.658     -17.658      -1.962      -1.962
RESULT: PASS
```

The angular branch has no LAMMPS analogue (no built-in constant body torque) and
is validated against theory (derived directly from the documented torque
formula) plus unit tests in `dirt_fixes`.

## Reference

P.A. Cundall, *Distinct element models of rock and soil structure* (1987); LAMMPS
`doc/src/fix_damping_cundall.rst` and `src/GRANULAR/fix_damping_cundall.cpp`; as
implemented in Yade-DEM and PFC.
