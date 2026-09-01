# CubeCL migration — CPU baseline

**Purpose.** This file records what the DIRT benchmark suite does *before* any
CubeCL work starts, so a future GPU port can be diffed against it rather than
against memory. It is a measurement artifact, not a design document; the plan of
record is [`docs/CUBECL-MIGRATION.md`](../docs/CUBECL-MIGRATION.md).

**What is and is not in here.** Every wall-clock and every physics number below
was observed on this machine on the date given. Nothing is estimated,
extrapolated, or copied from an earlier run. Where a bench was not executed for
this baseline, the row says so and cites `VALIDATION.md`/`VERIFICATION.md` as a
*documented* status rather than an observed one — those two are different kinds
of claim and are kept separate on purpose.

**No GPU on this box.** `elizabeth-hpc` has an integrated Radeon with no ROCm and
no CUDA. Everything here is the scalar CPU path. Nothing in this file is
evidence about GPU speed, and a future port must not be compared against it as
if it were.

---

## Provenance

| Item | Value |
|---|---|
| Date captured | 2026-09-01 (UTC) |
| Host | `elizabeth-hpc`, AMD Ryzen 9 9955HX (16C/32T), 29.9 GB RAM |
| OS / kernel | Ubuntu 26.04 LTS, Linux 7.0.0-27-generic |
| `rustc` | `1.96.1 (31fca3adb 2026-06-26)` |
| `cargo` | `1.96.1 (356927216 2026-06-26)` |
| dirt rev | `e7f6bc6d54066bdf58bbc70def1cf5cc290ebe0f` (branch `auto/cubecl-00-baseline`, identical tree to `main` tip `e7f6bc6`) |
| soil rev | `397c1a2510773cc268bfc5953e0a28cdf6e3cfcb` — the rev **pinned** in `Cargo.toml`. Note `~/projects/soil` `main` is at `aaa39caac2ff5066d998c30165ae0422c496899e`; the pin is what actually builds. |
| grass rev | `83198be92efa43320b32f986d7dd72550e3bc141` (resolved from `branch = "main"` in `Cargo.lock`; `Cargo.lock` is gitignored, so this is *not* pinned and can move) |
| Python driver | `$BENCH_PYTHON` = `~/.venvs/bench/bin/python3`, Python 3.14.4 |
| LAMMPS (cross-code legs) | `lmp_serial`, stable `22Jul2025_update4` |
| Build profile | `[profile.release]` `opt-level=3`, `lto="fat"`, `codegen-units=1`, `panic="abort"` |
| Parallelism | serial, single rank, single thread. There is no `rayon` in the stack; MPI is off in every run below. |

`source ~/projects/.build-env` is required before any of the commands in this
file (it puts `cargo`, `$BENCH_PYTHON` and `lmp_serial` on the path and fixes
bindgen's include path for the MPI bindings).

---

## Build baseline (observed)

Two builds, both clean, both through `qsub -q build`:

| Feature set | Command | Result |
|---|---|---|
| default (`mpi_backend`, `precision-double`) | `cargo build --release --examples` | **PASS**, `Finished release profile in 1m 15s` (cold `target/`) |
| sweep feature set (`--no-default-features --features precision-double`) | `cargo build --release --examples --no-default-features --features precision-double` | **PASS**, `Finished release profile in 56.31s` (warm registry, cold objects for this feature set) |

The second is the one that matters for the benchmark suite: 49 of the 51
`sweep.py` drivers under `examples/` build with
`--no-default-features --features precision-double`. Only
`bench_mpi_decomposition` and `bond_mpi_drift` use the default (MPI) feature
set.

Exact submission form used:

```
qsub -N cubecl00-build -q build -- bash -lc \
  "source ~/projects/.build-env && cd <repo> && cargo build --release --examples"
qsub -N cubecl00-build-nodefault -q build -- bash -lc \
  "source ~/projects/.build-env && cd <repo> && \
   cargo build --release --examples --no-default-features --features precision-double"
```

---

## The CI gate is currently red before it runs anything

`examples/ci_validation.py` is the repo's own regression gate. As of this rev it
**exits 2 without running a single sweep**, because its manifest check finds
eight `sweep.py` drivers that are neither in `FULL_SWEEPS` nor in
`DOCUMENTED_EXCLUSIONS`:

```
$ $BENCH_PYTHON examples/ci_validation.py --list ; echo $?
CI validation manifest is inconsistent:
  MISSING examples/bench_fiber_timestep/sweep.py
  MISSING examples/bench_haff_ensemble/sweep.py
  MISSING examples/bench_potyondy_cundall_bpm/sweep.py
  MISSING examples/bench_rod_shear_aspect_ratio/sweep.py
  MISSING examples/material_pair_table_validation/sweep.py
  MISSING examples/particleswith_migration_validation/sweep.py
  MISSING examples/simulation_fixture_validation/sweep.py
  MISSING examples/typed_schedule_validation/sweep.py
2
```

This matters for the migration: the intended "run the suite with the feature on
and off and show agreement" gate of stage 4 does not currently have a working
single entry point. Every measurement in this file was therefore taken by
invoking each `sweep.py` directly, not through `ci_validation.py`. Fixing the
manifest is proposed as a separate goal; it is deliberately **not** done here,
because this goal is forbidden from changing anything.

---

## How a bench is run

Every driver is a self-contained Python file that does `generate` (write
configs) → `start` (build the example, run each case) → `graph` (validate
against the reference and plot). Running it with no argument does all three and
exits non-zero if any gate fails. So the canonical queue command for any bench
`<b>` is:

```
qsub -N <b> -q sim -- bash -lc \
  "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/<b>/sweep.py"
```

with two exceptions that need the MPI (default) feature set and multiple ranks:
`bench_mpi_decomposition` and `bond_mpi_drift`. The drivers build the example
themselves, so a prior `qsub -q build` only saves time; it is not required.

`bench_guo2018_fiber_shear_cell` has **no `sweep.py`** — it is a set of audit
scripts (`campaign_preflight.py`, `evidence_contract.py`, …) rather than a
runnable benchmark, and has no single command.

---

## The five measured baselines

These five were chosen as the cheapest set that still spans the physics a GPU
port would touch first: linear normal contact, Hertzian normal contact against a
wall, a global damping fix, the tangential/Coulomb wall law, and the twisting
(rotational-history) law. All five were run through `qsub -q sim` on
2026-09-01 with a release build already warm, so the wall-clock below is
simulation + validation + plotting, with a no-op `cargo build` inside it.

Wall-clock was measured inside the job as `date +%s` around the driver, so it
excludes queue wait. All runs are single rank, single thread, `precision-double`.

| Bench | qsub job | cases | particles/case | steps/case | wall-clock | exit | checks |
|---|---|---|---|---|---|---|---|
| `bench_hooke_rebound` | 3599 | 20 | 2 spheres | 21 304 – 55 732 (`dt = 1.016633e-07 s`) | **3 s** | 0 | 71/71 |
| `bench_hertz_rebound` | 3600 | 20 | 1 sphere + 1 plane wall | 15 903 – 128 070 (default `dt`, not set in config) | **13 s** | 0 | 66/66 + 20 LAMMPS |
| `bench_cundall_damping` | 3601 | 3 | 1 sphere | 6 000 (`dt = 1e-04 s`) | **2 s** | 0 | 12 rates + 6 LAMMPS |
| `bench_twisting_friction` | 3602 | 6 | 2 spheres (one frozen) | 20 000 (`dt = 1e-05 s`) | **2 s** | 0 | 6×3 gates |
| `bench_sliding_friction` | 3603 | 6 | 1 sphere + 1 plane wall | 25 800 – 119 000 (`dt = 2e-06 s`) | **24 s** | 0 | 6×3 + 6 LAMMPS |

Case/step/particle counts were read out of the generated
`examples/<b>/sweep/<case>/config.toml` files, not from the drivers' source.

The exact commands, one per bench:

```
qsub -N bl-bench_hooke_rebound     -q sim -- bash -lc "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/bench_hooke_rebound/sweep.py"
qsub -N bl-bench_hertz_rebound     -q sim -- bash -lc "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/bench_hertz_rebound/sweep.py"
qsub -N bl-bench_cundall_damping   -q sim -- bash -lc "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/bench_cundall_damping/sweep.py"
qsub -N bl-bench_twisting_friction -q sim -- bash -lc "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/bench_twisting_friction/sweep.py"
qsub -N bl-bench_sliding_friction  -q sim -- bash -lc "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/bench_sliding_friction/sweep.py"
```

What follows is the **verbatim** driver output for each, with the pueue header
and cargo build noise removed and the absolute repo path replaced by `<repo>`.
It is reproduced in full rather than summarised because every asserted physics
number is a thing a GPU port can silently change, and a summary would decide in
advance which of them mattered.

### `bench_hooke_rebound` — linear-spring normal rebound (job 3599, 3 s)

```
Generated 20 configs under <repo>/examples/bench_hooke_rebound/sweep
  [ 1/20] v=0.5  COR=0.3   COR=0.2998
  [ 2/20] v=1.0  COR=0.3   COR=0.2998
  [ 3/20] v=2.0  COR=0.3   COR=0.2998
  [ 4/20] v=4.0  COR=0.3   COR=0.3001
  [ 5/20] v=0.5  COR=0.5   COR=0.4999
  [ 6/20] v=1.0  COR=0.5   COR=0.4999
  [ 7/20] v=2.0  COR=0.5   COR=0.4999
  [ 8/20] v=4.0  COR=0.5   COR=0.5002
  [ 9/20] v=0.5  COR=0.7   COR=0.7000
  [10/20] v=1.0  COR=0.7   COR=0.7000
  [11/20] v=2.0  COR=0.7   COR=0.7000
  [12/20] v=4.0  COR=0.7   COR=0.7000
  [13/20] v=0.5  COR=0.9   COR=0.9000
  [14/20] v=1.0  COR=0.9   COR=0.8999
  [15/20] v=2.0  COR=0.9   COR=0.9000
  [16/20] v=4.0  COR=0.9   COR=0.9000
  [17/20] v=0.5  COR=1.0   COR=1.0000
  [18/20] v=1.0  COR=1.0   COR=1.0000
  [19/20] v=2.0  COR=1.0   COR=1.0000
  [20/20] v=4.0  COR=1.0   COR=1.0000

DIRT: 20/20 cases -> <repo>/examples/bench_hooke_rebound/data/sweep_results.csv

========================================================================
Hooke (linear-spring) Normal-Contact Rebound Benchmark
========================================================================
  Two identical spheres  R=0.005 m  rho=2500.0 kg/m^3
  kn=1.000e+05 N/m  kt=2.857e+04 N/m   m_eff=6.5450e-04 kg  omega_0=1.2361e+04 rad/s
  Reference: exact linear damped-oscillator collision (theory only).
  Tolerances: COR 1%  contact-time 2%  overlap 2%

  e=0.30  v= 0.50 m/s:
    COR:    0.29985 vs 0.30000   (err  0.02%)  [PASS]
    t_c:    2.72153e-04 vs 2.72183e-04 s (err  0.01%)  [PASS]
    d_max:  2.54821e-05 vs 2.54916e-05 m (err  0.04%)  [PASS]
  e=0.30  v= 1.00 m/s:
    COR:    0.29985 vs 0.30000   (err  0.02%)  [PASS]
    t_c:    2.72153e-04 vs 2.72183e-04 s (err  0.01%)  [PASS]
    d_max:  5.09543e-05 vs 5.09831e-05 m (err  0.06%)  [PASS]
  e=0.30  v= 2.00 m/s:
    COR:    0.29985 vs 0.30000   (err  0.02%)  [PASS]
    t_c:    2.72153e-04 vs 2.72183e-04 s (err  0.01%)  [PASS]
    d_max:  1.01945e-04 vs 1.01966e-04 m (err  0.02%)  [PASS]
  e=0.30  v= 4.00 m/s:
    COR:    0.30012 vs 0.30000   (err  0.01%)  [PASS]
    t_c:    2.72051e-04 vs 2.72183e-04 s (err  0.05%)  [PASS]
    d_max:  2.03925e-04 vs 2.03932e-04 m (err  0.00%)  [PASS]
  e=0.50  v= 0.50 m/s:
    COR:    0.49990 vs 0.50000   (err  0.01%)  [PASS]
    t_c:    2.60258e-04 vs 2.60271e-04 s (err  0.00%)  [PASS]
    d_max:  3.00004e-05 vs 3.00066e-05 m (err  0.02%)  [PASS]
  e=0.50  v= 1.00 m/s:
    COR:    0.49990 vs 0.50000   (err  0.01%)  [PASS]
    t_c:    2.60258e-04 vs 2.60271e-04 s (err  0.00%)  [PASS]
    d_max:  5.99938e-05 vs 6.00133e-05 m (err  0.03%)  [PASS]
  e=0.50  v= 2.00 m/s:
    COR:    0.49990 vs 0.50000   (err  0.01%)  [PASS]
    t_c:    2.60258e-04 vs 2.60271e-04 s (err  0.00%)  [PASS]
    d_max:  1.20013e-04 vs 1.20027e-04 m (err  0.01%)  [PASS]
  e=0.50  v= 4.00 m/s:
    COR:    0.50017 vs 0.50000   (err  0.02%)  [PASS]
    t_c:    2.60156e-04 vs 2.60271e-04 s (err  0.04%)  [PASS]
    d_max:  2.40052e-04 vs 2.40053e-04 m (err  0.00%)  [PASS]
  e=0.70  v= 0.50 m/s:
    COR:    0.69996 vs 0.70000   (err  0.00%)  [PASS]
    t_c:    2.55785e-04 vs 2.55791e-04 s (err  0.00%)  [PASS]
    d_max:  3.42771e-05 vs 3.42805e-05 m (err  0.01%)  [PASS]
  e=0.70  v= 1.00 m/s:
    COR:    0.69996 vs 0.70000   (err  0.00%)  [PASS]
    t_c:    2.55785e-04 vs 2.55791e-04 s (err  0.00%)  [PASS]
    d_max:  6.85499e-05 vs 6.85610e-05 m (err  0.02%)  [PASS]
  e=0.70  v= 2.00 m/s:
    COR:    0.69996 vs 0.70000   (err  0.00%)  [PASS]
    t_c:    2.55785e-04 vs 2.55791e-04 s (err  0.00%)  [PASS]
    d_max:  1.37115e-04 vs 1.37122e-04 m (err  0.01%)  [PASS]
  e=0.70  v= 4.00 m/s:
    COR:    0.69996 vs 0.70000   (err  0.00%)  [PASS]
    t_c:    2.55785e-04 vs 2.55791e-04 s (err  0.00%)  [PASS]
    d_max:  2.74245e-04 vs 2.74244e-04 m (err  0.00%)  [PASS]
  e=0.90  v= 0.50 m/s:
    COR:    0.90002 vs 0.90000   (err  0.00%)  [PASS]
    t_c:    2.54260e-04 vs 2.54301e-04 s (err  0.02%)  [PASS]
    d_max:  3.84169e-05 vs 3.84179e-05 m (err  0.00%)  [PASS]
  e=0.90  v= 1.00 m/s:
    COR:    0.89995 vs 0.90000   (err  0.01%)  [PASS]
    t_c:    2.54362e-04 vs 2.54301e-04 s (err  0.02%)  [PASS]
    d_max:  7.68323e-05 vs 7.68358e-05 m (err  0.00%)  [PASS]
  e=0.90  v= 2.00 m/s:
    COR:    0.90002 vs 0.90000   (err  0.00%)  [PASS]
    t_c:    2.54260e-04 vs 2.54301e-04 s (err  0.02%)  [PASS]
    d_max:  1.53670e-04 vs 1.53672e-04 m (err  0.00%)  [PASS]
  e=0.90  v= 4.00 m/s:
    COR:    0.90002 vs 0.90000   (err  0.00%)  [PASS]
    t_c:    2.54260e-04 vs 2.54301e-04 s (err  0.02%)  [PASS]
    d_max:  3.07344e-04 vs 3.07343e-04 m (err  0.00%)  [PASS]
  e=1.00  v= 0.50 m/s:
    COR:    1.00000 vs 1.00000   (err  0.00%)  [PASS]
    t_c:    2.54158e-04 vs 2.54158e-04 s (err  0.00%)  [PASS]
    d_max:  4.04505e-05 vs 4.04505e-05 m (err  0.00%)  [PASS]
  e=1.00  v= 1.00 m/s:
    COR:    1.00000 vs 1.00000   (err  0.00%)  [PASS]
    t_c:    2.54158e-04 vs 2.54158e-04 s (err  0.00%)  [PASS]
    d_max:  8.09011e-05 vs 8.09011e-05 m (err  0.00%)  [PASS]
  e=1.00  v= 2.00 m/s:
    COR:    1.00000 vs 1.00000   (err  0.00%)  [PASS]
    t_c:    2.54158e-04 vs 2.54158e-04 s (err  0.00%)  [PASS]
    d_max:  1.61802e-04 vs 1.61802e-04 m (err  0.00%)  [PASS]
  e=1.00  v= 4.00 m/s:
    COR:    1.00000 vs 1.00000   (err  0.00%)  [PASS]
    t_c:    2.54158e-04 vs 2.54158e-04 s (err  0.00%)  [PASS]
    d_max:  3.23604e-04 vs 3.23604e-04 m (err  0.00%)  [PASS]

Velocity-independence (linear-contact signature):
  e=0.30:  COR spread 0.0003 (tol 0.005)  [PASS]   t_c spread 0.04% (tol 1%)  [PASS]
  e=0.50:  COR spread 0.0003 (tol 0.005)  [PASS]   t_c spread 0.04% (tol 1%)  [PASS]
  e=0.70:  COR spread 0.0000 (tol 0.005)  [PASS]   t_c spread 0.00% (tol 1%)  [PASS]
  e=0.90:  COR spread 0.0001 (tol 0.005)  [PASS]   t_c spread 0.04% (tol 1%)  [PASS]
  e=1.00:  COR spread 0.0000 (tol 0.005)  [PASS]   t_c spread 0.00% (tol 1%)  [PASS]

Completeness: 20/20 cases  [PASS]

Overall: 71/71 checks passed
ALL CHECKS PASSED
Saved: <repo>/examples/bench_hooke_rebound/plots/cor_validation.png
Saved: <repo>/examples/bench_hooke_rebound/plots/contact_duration.png
Saved: <repo>/examples/bench_hooke_rebound/plots/peak_overlap.png
### bench_hooke_rebound EXIT=0 WALLCLOCK_S=3
```

### `bench_hertz_rebound` — Hertz/Tsuji normal rebound vs ODE and LAMMPS (job 3600, 13 s)

```
Generated 20 DIRT configs under <repo>/examples/bench_hertz_rebound/sweep
LAMMPS: /home/suehr/.local/bin/lmp_serial
  [ 1/20] v0=0.1  COR=0.5    DIRT COR=0.6162
  [ 2/20] v0=0.5  COR=0.5    DIRT COR=0.6143
  [ 3/20] v0=1.0  COR=0.5    DIRT COR=0.6135
  [ 4/20] v0=2.0  COR=0.5    DIRT COR=0.6126
  [ 5/20] v0=0.1  COR=0.7    DIRT COR=0.7777
  [ 6/20] v0=0.5  COR=0.7    DIRT COR=0.7780
  [ 7/20] v0=1.0  COR=0.7    DIRT COR=0.7757
  [ 8/20] v0=2.0  COR=0.7    DIRT COR=0.7754
  [ 9/20] v0=0.1  COR=0.9    DIRT COR=0.9282
  [10/20] v0=0.5  COR=0.9    DIRT COR=0.9286
  [11/20] v0=1.0  COR=0.9    DIRT COR=0.9278
  [12/20] v0=2.0  COR=0.9    DIRT COR=0.9280
  [13/20] v0=0.1  COR=0.95   DIRT COR=0.9651
  [14/20] v0=0.5  COR=0.95   DIRT COR=0.9652
  [15/20] v0=1.0  COR=0.95   DIRT COR=0.9648
  [16/20] v0=2.0  COR=0.95   DIRT COR=0.9649
  [17/20] v0=0.1  COR=1.0    DIRT COR=1.0000
  [18/20] v0=0.5  COR=1.0    DIRT COR=1.0000
  [19/20] v0=1.0  COR=1.0    DIRT COR=1.0000
  [20/20] v0=2.0  COR=1.0    DIRT COR=1.0000

DIRT:   20/20 cases -> <repo>/examples/bench_hertz_rebound/data/sweep_results.csv
LAMMPS: 20/20 cases -> <repo>/examples/bench_hertz_rebound/data/lammps_results.csv

=================================================================
Hertz Contact Rebound Benchmark Validation
=================================================================
  E* = 3.678e+10 Pa
  m  = 1.308997e-03 kg
  R  = 5.0 mm

v0=0.1 m/s, COR_in=0.50:
  COR:     0.6162 vs Tsuji ODE 0.6161  (err=0.0%)  [PASS]
  t_c:     5.793e-05 vs Tsuji ODE 5.803e-05 s  (err=0.2%)  [PASS]
  d_max:   1.557e-06 vs Tsuji ODE 1.562e-06 m  (err=0.3%)  [PASS]
v0=0.5 m/s, COR_in=0.50:
  COR:     0.6143 vs Tsuji ODE 0.6161  (err=0.3%)  [PASS]
  t_c:     4.192e-05 vs Tsuji ODE 4.206e-05 s  (err=0.3%)  [PASS]
  d_max:   5.636e-06 vs Tsuji ODE 5.659e-06 m  (err=0.4%)  [PASS]
v0=1.0 m/s, COR_in=0.50:
  COR:     0.6135 vs Tsuji ODE 0.6161  (err=0.4%)  [PASS]
  t_c:     3.659e-05 vs Tsuji ODE 3.662e-05 s  (err=0.1%)  [PASS]
  d_max:   9.792e-06 vs Tsuji ODE 9.853e-06 m  (err=0.6%)  [PASS]
v0=2.0 m/s, COR_in=0.50:
  COR:     0.6126 vs Tsuji ODE 0.6161  (err=0.6%)  [PASS]
  t_c:     3.201e-05 vs Tsuji ODE 3.188e-05 s  (err=0.4%)  [PASS]
  d_max:   1.702e-05 vs Tsuji ODE 1.716e-05 m  (err=0.8%)  [PASS]
v0=0.1 m/s, COR_in=0.70:
  COR:     0.7777 vs Tsuji ODE 0.7777  (err=0.0%)  [PASS]
  t_c:     5.641e-05 vs Tsuji ODE 5.629e-05 s  (err=0.2%)  [PASS]
  d_max:   1.688e-06 vs Tsuji ODE 1.691e-06 m  (err=0.2%)  [PASS]
v0=0.5 m/s, COR_in=0.70:
  COR:     0.7780 vs Tsuji ODE 0.7777  (err=0.0%)  [PASS]
  t_c:     4.040e-05 vs Tsuji ODE 4.080e-05 s  (err=1.0%)  [PASS]
  d_max:   6.113e-06 vs Tsuji ODE 6.127e-06 m  (err=0.2%)  [PASS]
v0=1.0 m/s, COR_in=0.70:
  COR:     0.7757 vs Tsuji ODE 0.7777  (err=0.3%)  [PASS]
  t_c:     3.583e-05 vs Tsuji ODE 3.552e-05 s  (err=0.9%)  [PASS]
  d_max:   1.063e-05 vs Tsuji ODE 1.067e-05 m  (err=0.3%)  [PASS]
v0=2.0 m/s, COR_in=0.70:
  COR:     0.7754 vs Tsuji ODE 0.7777  (err=0.3%)  [PASS]
  t_c:     3.125e-05 vs Tsuji ODE 3.092e-05 s  (err=1.1%)  [PASS]
  d_max:   1.851e-05 vs Tsuji ODE 1.857e-05 m  (err=0.4%)  [PASS]
v0=0.1 m/s, COR_in=0.90:
  COR:     0.9282 vs Tsuji ODE 0.9284  (err=0.0%)  [PASS]
  t_c:     5.564e-05 vs Tsuji ODE 5.517e-05 s  (err=0.9%)  [PASS]
  d_max:   1.805e-06 vs Tsuji ODE 1.806e-06 m  (err=0.1%)  [PASS]
v0=0.5 m/s, COR_in=0.90:
  COR:     0.9286 vs Tsuji ODE 0.9284  (err=0.0%)  [PASS]
  t_c:     3.964e-05 vs Tsuji ODE 3.998e-05 s  (err=0.9%)  [PASS]
  d_max:   6.543e-06 vs Tsuji ODE 6.546e-06 m  (err=0.1%)  [PASS]
v0=1.0 m/s, COR_in=0.90:
  COR:     0.9278 vs Tsuji ODE 0.9284  (err=0.1%)  [PASS]
  t_c:     3.506e-05 vs Tsuji ODE 3.481e-05 s  (err=0.7%)  [PASS]
  d_max:   1.139e-05 vs Tsuji ODE 1.140e-05 m  (err=0.1%)  [PASS]
v0=2.0 m/s, COR_in=0.90:
  COR:     0.9280 vs Tsuji ODE 0.9284  (err=0.0%)  [PASS]
  t_c:     3.049e-05 vs Tsuji ODE 3.030e-05 s  (err=0.6%)  [PASS]
  d_max:   1.982e-05 vs Tsuji ODE 1.984e-05 m  (err=0.1%)  [PASS]
v0=0.1 m/s, COR_in=0.95:
  COR:     0.9651 vs Tsuji ODE 0.9652  (err=0.0%)  [PASS]
  t_c:     5.564e-05 vs Tsuji ODE 5.494e-05 s  (err=1.3%)  [PASS]
  d_max:   1.834e-06 vs Tsuji ODE 1.834e-06 m  (err=0.0%)  [PASS]
v0=0.5 m/s, COR_in=0.95:
  COR:     0.9652 vs Tsuji ODE 0.9652  (err=0.0%)  [PASS]
  t_c:     3.964e-05 vs Tsuji ODE 3.982e-05 s  (err=0.5%)  [PASS]
  d_max:   6.644e-06 vs Tsuji ODE 6.646e-06 m  (err=0.0%)  [PASS]
v0=1.0 m/s, COR_in=0.95:
  COR:     0.9648 vs Tsuji ODE 0.9652  (err=0.0%)  [PASS]
  t_c:     3.506e-05 vs Tsuji ODE 3.467e-05 s  (err=1.1%)  [PASS]
  d_max:   1.157e-05 vs Tsuji ODE 1.157e-05 m  (err=0.0%)  [PASS]
v0=2.0 m/s, COR_in=0.95:
  COR:     0.9649 vs Tsuji ODE 0.9652  (err=0.0%)  [PASS]
  t_c:     3.049e-05 vs Tsuji ODE 3.018e-05 s  (err=1.0%)  [PASS]
  d_max:   2.013e-05 vs Tsuji ODE 2.015e-05 m  (err=0.1%)  [PASS]
v0=0.1 m/s, COR_in=1.00:
  COR:     1.0000 vs Tsuji ODE 1.0000  (err=0.0%)  [PASS]
  t_c:     5.488e-05 vs Tsuji ODE 5.475e-05 s  (err=0.2%)  [PASS]
  d_max:   1.860e-06 vs Tsuji ODE 1.860e-06 m  (err=0.0%)  [PASS]
v0=0.5 m/s, COR_in=1.00:
  COR:     1.0000 vs Tsuji ODE 1.0000  (err=0.0%)  [PASS]
  t_c:     3.964e-05 vs Tsuji ODE 3.968e-05 s  (err=0.1%)  [PASS]
  d_max:   6.739e-06 vs Tsuji ODE 6.741e-06 m  (err=0.0%)  [PASS]
v0=1.0 m/s, COR_in=1.00:
  COR:     1.0000 vs Tsuji ODE 1.0000  (err=0.0%)  [PASS]
  t_c:     3.506e-05 vs Tsuji ODE 3.454e-05 s  (err=1.5%)  [PASS]
  d_max:   1.173e-05 vs Tsuji ODE 1.174e-05 m  (err=0.0%)  [PASS]
v0=2.0 m/s, COR_in=1.00:
  COR:     1.0000 vs Tsuji ODE 1.0000  (err=0.0%)  [PASS]
  t_c:     3.049e-05 vs Tsuji ODE 3.007e-05 s  (err=1.4%)  [PASS]
  d_max:   2.043e-05 vs Tsuji ODE 2.043e-05 m  (err=0.0%)  [PASS]

Velocity-independence (realized COR vs impact speed):
  nominal 0.50: realized 0.6142, spread 0.0037 (tol 0.01)  [PASS]
  nominal 0.70: realized 0.7767, spread 0.0026 (tol 0.01)  [PASS]
  nominal 0.90: realized 0.9282, spread 0.0008 (tol 0.01)  [PASS]
  nominal 0.95: realized 0.9650, spread 0.0003 (tol 0.01)  [PASS]
  nominal 1.00: realized 1.0000, spread 0.0001 (tol 0.01)  [PASS]

Completeness: 20/20 cases  [PASS]

Velocity-independence: 5/5 passed
COR vs ODE checks:     20/20 passed
Contact time checks:   20/20 passed
Overlap checks:        20/20 passed

Overall: 66/66 checks passed
ALL CHECKS PASSED

==================================================
Realized COR: DIRT vs LAMMPS (same nominal input, both Tsuji)
==================================================
     v0   COR |     DIRT   LAMMPS |     diff
    0.1   0.5 |   0.6162   0.6147 |  -0.0015
    0.5   0.5 |   0.6143   0.6142 |  -0.0001
    1.0   0.5 |   0.6135   0.6136 |  +0.0001
    2.0   0.5 |   0.6126   0.6127 |  +0.0001
    0.1   0.7 |   0.7777   0.7769 |  -0.0008
    0.5   0.7 |   0.7780   0.7780 |  +0.0000
    1.0   0.7 |   0.7757   0.7756 |  -0.0000
    2.0   0.7 |   0.7754   0.7754 |  +0.0001
    0.1   0.9 |   0.9282   0.9284 |  +0.0002
    0.5   0.9 |   0.9286   0.9286 |  +0.0000
    1.0   0.9 |   0.9278   0.9278 |  +0.0000
    2.0   0.9 |   0.9280   0.9280 |  +0.0000
    0.1  0.95 |   0.9651   0.9651 |  +0.0000
    0.5  0.95 |   0.9652   0.9652 |  +0.0000
    1.0  0.95 |   0.9648   0.9648 |  +0.0000
    2.0  0.95 |   0.9649   0.9649 |  +0.0000
    0.1   1.0 |   1.0000   0.9996 |  -0.0005
    0.5   1.0 |   1.0000   0.9996 |  -0.0004
    1.0   1.0 |   1.0000   0.9996 |  -0.0004
    2.0   1.0 |   1.0000   0.9996 |  -0.0004

Max |DIRT - LAMMPS| COR = 0.0015  (tol 0.005)  [PASS]

Saved: <repo>/examples/bench_hertz_rebound/plots/cor_validation.png
Saved: <repo>/examples/bench_hertz_rebound/plots/contact_duration.png
Saved: <repo>/examples/bench_hertz_rebound/plots/peak_overlap.png
### bench_hertz_rebound EXIT=0 WALLCLOCK_S=13
```

### `bench_cundall_damping` — Cundall non-viscous damping (job 3601, 2 s)

```
Generated 3 DIRT sweep configs under <repo>/examples/bench_cundall_damping/sweep
LAMMPS: /home/suehr/.local/bin/lmp_serial
  [1/3] gamma=0.2   DIRT a_up=-11.772/-11.772 a_dn=-7.848/-7.848 al_dn=-91.7/-91.7 al_up=-61.1/-61.1   LAMMPS a_up=-11.772 a_dn=-7.848
  [2/3] gamma=0.5   DIRT a_up=-14.715/-14.715 a_dn=-4.905/-4.905 al_dn=-114.6/-114.6 al_up=-38.2/-38.2   LAMMPS a_up=-14.715 a_dn=-4.905
  [3/3] gamma=0.8   DIRT a_up=-17.658/-17.658 a_dn=-1.962/-1.962 al_dn=-137.5/-137.5 al_up=-15.3/-15.3   LAMMPS a_up=-17.658 a_dn=-1.962

DIRT:   3/3 cases -> <repo>/examples/bench_cundall_damping/data/sweep.csv
LAMMPS: 3/3 cases -> <repo>/examples/bench_cundall_damping/data/sweep_lammps.csv


=== Cundall non-viscous damping validation ===
  g=9.81  V0=3.0  Omega0=5.0  Tz=-1e-06
  exact: a_up=-g(1+gl)  a_down=-g(1-gl)  alpha_down=Tz(1+ga)/I  alpha_up=Tz(1-ga)/I
   gamma              a_up            a_down          alpha_down            alpha_up  note
    0.20   -11.772/-11.772    -7.848/-7.848    -91.673/-91.673   -61.115/-61.115  
    0.50   -14.715/-14.715    -4.905/-4.905  -114.592/-114.592   -38.197/-38.197  
    0.80   -17.658/-17.658    -1.962/-1.962  -137.510/-137.510   -15.279/-15.279  

  tolerance: <= 1% relative error on all four rates, all gamma
RESULT: PASS

=== DIRT vs LAMMPS fix damping/cundall (linear a_up / a_down) ===
   gamma   DIRT a_up    LMP a_up   DIRT a_dn    LMP a_dn
    0.20     -11.772     -11.772      -7.848      -7.848
    0.50     -14.715     -14.715      -4.905      -4.905
    0.80     -17.658     -17.658      -1.962      -1.962

Figures -> <repo>/examples/bench_cundall_damping/plots/cundall_traces.png, cundall_rates.png
### bench_cundall_damping EXIT=0 WALLCLOCK_S=2
```

### `bench_twisting_friction` — constant and SDS twisting spin-down (job 3602, 2 s)

```
Generated 6 DIRT sweep configs under <repo>/examples/bench_twisting_friction/sweep
  [1/6] constant mu_tw=0.05   a_fit=122.625 a_pred=122.625 rel=0.00% perp=0.0e+00 drift=0.0e+00
  [2/6] constant mu_tw=0.1    a_fit=245.250 a_pred=245.250 rel=0.00% perp=0.0e+00 drift=0.0e+00
  [3/6] constant mu_tw=0.2    a_fit=490.501 a_pred=490.500 rel=0.00% perp=0.0e+00 drift=0.0e+00
  [4/6] sds      mu_tw=0.05   a_fit=122.625 a_pred=122.625 rel=0.00% perp=0.0e+00 drift=0.0e+00
  [5/6] sds      mu_tw=0.1    a_fit=245.250 a_pred=245.250 rel=0.00% perp=0.0e+00 drift=0.0e+00
  [6/6] sds      mu_tw=0.2    a_fit=490.501 a_pred=490.500 rel=0.00% perp=0.0e+00 drift=0.0e+00

DIRT: 6/6 cases -> <repo>/examples/bench_twisting_friction/data/sweep.csv


=== Twisting-friction spin-down validation ===
  R=0.005 m  equal spheres (r_eff = R/2)  g=9.81  omega0=8.0
  model: alpha = (5/4) mu_tw g / R   (exact)
      model  mu_tw     a_fit    a_pred  rel_err  max_perp     drift  note
   constant  0.050   122.625   122.625    0.00%  0.00e+00  0.00e+00  
   constant  0.100   245.250   245.250    0.00%  0.00e+00  0.00e+00  
   constant  0.200   490.501   490.500    0.00%  0.00e+00  0.00e+00  
        sds  0.050   122.625   122.625    0.00%  0.00e+00  0.00e+00  
        sds  0.100   245.250   245.250    0.00%  0.00e+00  0.00e+00  
        sds  0.200   490.501   490.500    0.00%  0.00e+00  0.00e+00  

  tolerances: slope <= 3% rel, omega_perp <= 0.1% of omega0, drift <= 10 um
RESULT: PASS

Figures -> <repo>/examples/bench_twisting_friction/plots/twist_spindown.png, spindown_vs_mu_tw.png
### bench_twisting_friction EXIT=0 WALLCLOCK_S=2
```

### `bench_sliding_friction` — slip-to-roll transition vs LAMMPS (job 3603, 24 s)

```
Generated 6 DIRT sweep configs under <repo>/examples/bench_sliding_friction/sweep
LAMMPS: /home/suehr/.local/bin/lmp_serial
  [ 1/6] mu_0.2_v0_1       a_fit=1.962  a_th=1.962  v_final=0.7143 (th 0.7143)   | LAMMPS a_fit=1.962  v_final=0.7141
  [ 2/6] mu_0.3_v0_1       a_fit=2.943  a_th=2.943  v_final=0.7143 (th 0.7143)   | LAMMPS a_fit=2.943  v_final=0.7139
  [ 3/6] mu_0.5_v0_1       a_fit=4.905  a_th=4.905  v_final=0.7143 (th 0.7143)   | LAMMPS a_fit=4.906  v_final=0.7137
  [ 4/6] mu_0.7_v0_1       a_fit=6.867  a_th=6.867  v_final=0.7143 (th 0.7143)   | LAMMPS a_fit=6.865  v_final=0.7135
  [ 5/6] mu_0.5_v0_0.5     a_fit=4.905  a_th=4.905  v_final=0.3571 (th 0.3571)   | LAMMPS a_fit=4.909  v_final=0.3566
  [ 6/6] mu_0.5_v0_1.5     a_fit=4.905  a_th=4.905  v_final=1.0714 (th 1.0714)   | LAMMPS a_fit=4.905  v_final=1.0720

DIRT:   6/6 cases -> <repo>/examples/bench_sliding_friction/data/sweep_summary.csv
LAMMPS: 6/6 cases -> <repo>/examples/bench_sliding_friction/data/lammps_results.csv


=== Sliding-friction validation ===
  g=9.81  R=0.005 m  (floor: dirt_wall z-plane at z=0)
     mu    v0    a_fit   a=mu g   err%    v_fin   5/7 v0   err%   t*meas     t*th   err%  note
   0.50  0.50    4.905    4.905    0.0   0.3571   0.3571    0.0   28.54m   29.12m    2.0  
   0.20  1.00    1.962    1.962    0.0   0.7143   0.7143    0.0  142.71m  145.62m    2.0  
   0.30  1.00    2.943    2.943    0.0   0.7143   0.7143    0.0   95.14m   97.08m    2.0  
   0.50  1.00    4.905    4.905    0.0   0.7143   0.7143    0.0   57.08m   58.25m    2.0  
   0.70  1.00    6.867    6.867    0.0   0.7143   0.7143    0.0   40.77m   41.61m    2.0  
   0.50  1.50    4.905    4.905    0.0   1.0714   1.0714    0.0   85.63m   87.37m    2.0  

  tolerances: a 8%, v_final 3%, t* 10%
RESULT: PASS

=== DIRT vs LAMMPS (fitted a, rolling plateau v_final) ===
     mu    v0   a_DIRT    a_LMP      d_a   vf_DIRT    vf_LMP     d_vf
   0.50  0.50    4.905    4.909   -0.004    0.3571    0.3566  +0.0006
   0.20  1.00    1.962    1.962   -0.000    0.7143    0.7141  +0.0002
   0.30  1.00    2.943    2.943   +0.000    0.7143    0.7139  +0.0003
   0.50  1.00    4.905    4.906   -0.001    0.7143    0.7137  +0.0006
   0.70  1.00    6.867    6.865   +0.002    0.7143    0.7135  +0.0008
   0.50  1.50    4.905    4.905   +0.000    1.0714    1.0720  -0.0006

=== DIRT vs LAMMPS cross-validation (Mindlin tangential) ===
     mu    v0   a_err%  vf_err%  note
   0.50  0.50     0.09     0.16  
   0.20  1.00     0.00     0.03  
   0.30  1.00     0.01     0.05  
   0.50  1.00     0.01     0.08  
   0.70  1.00     0.03     0.11  
   0.50  1.50     0.00     0.05  
  tolerance: 2% rel on a and v_final
XVAL: PASS

Figures -> <repo>/examples/bench_sliding_friction/plots/slip_to_roll.png, decel_vs_mu.png, vfinal_vs_v0.png
### bench_sliding_friction EXIT=0 WALLCLOCK_S=24
```

---

## Full bench inventory

Every directory matching `examples/bench_*`. "Observed" means it was executed on
this box on 2026-09-01 as part of capturing this baseline. The first batch of 12
cost probes ran under a 400 s `timeout` cap and the remaining 26 under 240 s;
the five detailed baselines ran under 900 s. Rows that hit their
cap say **TIMEOUT** and their true cost is only bounded below.

Cost probes captured wall-clock and exit status only — their full output was not
transcribed, so use the five detailed sections above for number-level diffing.

`cases` / `particles/case` / `steps/case` are read out of the `config.toml`
files each driver *generates*, not from the driver source. A `—` means that
driver does not configure its runs through a `config.toml` (it passes parameters
to the example binary directly), so no count could be read without executing
code — those cells are blank rather than guessed. `particles/case` is the sum of
`[[particles.insert]] count`, so it excludes particles read from a bond/data
file, which is why a few bonded-fiber benches show `—` for particles but a step
range.

Command for every row (except the MPI one, see below):

```
qsub -N <bench> -q sim -- bash -lc \
  "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/<bench>/sweep.py"
```

| Bench | wall-clock | result | run as | cases | particles/case | steps/case |
|---|---|---|---|---|---|---|
| `bench_angle_of_repose` | 400 s | **TIMEOUT** (hit cap) | cost probe | 27 | 1,200 | 100,000 |
| `bench_bond_breakage` | 22 s | **FAIL** (exit 1) | cost probe | 66 | — | 30,000–120,000 |
| `bench_chung_ooi_impact` | 2 s | **PASS** | cost probe | 10 | 1–2 | 3,286–21,908 |
| `bench_clump_haff_cooling` | 240 s | **TIMEOUT** (hit cap) | cost probe | — | — | — |
| `bench_clump_inertia_sampler` | 1 s | **PASS** | cost probe | — | — | — |
| `bench_clump_insertion_determinism` | 0 s | **PASS** | cost probe | — | — | — |
| `bench_column_collapse` | 400 s | **TIMEOUT** (hit cap) | cost probe | 33 | 110–1,100 | 80,000 |
| `bench_convergence` | 74 s | **PASS** | cost probe | 32 | 1–1,600 | 3,574–150,000 |
| `bench_cundall_damping` | 2 s | **PASS** | detailed | 3 | 1 | 6,000 |
| `bench_curtis_cantilever` | 11 s | **PASS** | cost probe | 5 | — | 350,000 |
| `bench_curtis_wet_fiber_breakage` | 1 s | **PASS** | cost probe | 4 | — | 4,500 |
| `bench_dmt_sjkr_cohesion` | 63 s | **PASS** | cost probe | 13 | 2 | 1,014,244 |
| `bench_fiber_crossover` | 2 s | **PASS** | cost probe | 7 | — | 30,000 |
| `bench_fiber_timestep` | 1 s | **PASS** | cost probe | — | — | — |
| `bench_granular_conductivity` | 121 s | **PASS** | cost probe | — | — | — |
| `bench_guo2018_fiber_shear_cell` | — | no `sweep.py` — not runnable as a bench | n/a | — | — | — |
| `bench_haff_ensemble` | 240 s | **TIMEOUT** (hit cap) | cost probe | — | — | — |
| `bench_hertz_rebound` | 13 s | **PASS** | detailed | 20 | 1 | 15,903–128,070 |
| `bench_hooke_rebound` | 3 s | **PASS** | detailed | 20 | 2 | 21,304–55,732 |
| `bench_hooke_wall_rebound` | 5 s | **PASS** | cost probe | 20 | 1 | 28,558–113,761 |
| `bench_hopper_beverloo` | 20 s | **PASS** | cost probe | 9 | 600–1,400 | 25,000–60,000 |
| `bench_jkr_adhesion` | 36 s | **PASS** | cost probe | 6 | 2 | 1,014,244 |
| `bench_kharaz_oblique` | 8 s | **PASS** | cost probe | 16 | 1 | 40,000 |
| `bench_lebc_shear` | 400 s | **TIMEOUT** (hit cap) | cost probe | 17 | 82–1,634 | 10,000–20,000 |
| `bench_liquid_bridge_cohesion` | 28 s | **PASS** | cost probe | 4 | 2–260 | 200–30,000 |
| `bench_marshall_twisting` | 2 s | **PASS** | cost probe | 3 | 2 | 40,000 |
| `bench_mdr_elastoplastic_normal` | 1 s | **PASS** | cost probe | — | — | — |
| `bench_mindlin_rescale_tangential` | 0 s | **PASS** | cost probe | — | — | — |
| `bench_mpi_decomposition` | 3 s | **PASS** | cost probe | — | — | — |
| `bench_nohistory_tangential` | 1 s | **PASS** | cost probe | — | — | — |
| `bench_oblique_impact` | 15 s | **PASS** | cost probe | 11 | 2 | 120,000 |
| `bench_plate_sinkage` | 240 s | **TIMEOUT** (hit cap) | cost probe | 4 | 2,400 | 45,000 |
| `bench_polydisperse_mixing` | 9 s | **PASS** | cost probe | 17 | 2 | 40,000–200,000 |
| `bench_potyondy_cundall_bpm` | 26 s | **PASS** | cost probe | — | — | — |
| `bench_restart_determinism` | 1 s | **PASS** | cost probe | — | — | — |
| `bench_rod_haff_cooling` | 240 s | **TIMEOUT** (hit cap) | cost probe | — | — | — |
| `bench_rod_shear_aspect_ratio` | 11 s | **PASS** | cost probe | 4 | 19 | 2,000 |
| `bench_rolling_decay` | 3 s | **PASS** | cost probe | 3 | 1 | 40,000 |
| `bench_sds_rolling` | 2 s | **PASS** | cost probe | 5 | 2 | 6,000–60,000 |
| `bench_sliding_friction` | 24 s | **PASS** | detailed | 6 | 1 | 25,800–119,000 |
| `bench_sphere_haff_cooling` | 24 s | **PASS** | cost probe | — | — | — |
| `bench_twisting_friction` | 2 s | **PASS** | detailed | 6 | 2 | 20,000 |
| `bench_wall_activate_by_name` | 0 s | **PASS** | cost probe | — | — | — |
| `bench_wall_twisting_parity` | 1 s | **PASS** | cost probe | — | — | — |

### The one that needs a different feature set

`bench_mpi_decomposition` builds with the **default** feature set
(`mpi_backend`, i.e. no `--no-default-features` in its driver) and compares
`1x1x1`, `2x1x1` and `2x2x1` decompositions. `ci_validation.py` excludes it from
the stock no-MPI green suite by documented decision. The command shape is the
same; the difference is inside `sweep.py`:

```
qsub -N bench_mpi_decomposition -q sim -- bash -lc \
  "source ~/projects/.build-env && cd <repo> && \$BENCH_PYTHON examples/bench_mpi_decomposition/sweep.py"
```

`bond_mpi_drift` is the other MPI-feature driver. It is not under
`examples/bench_*` so it is outside this inventory, but it is mentioned so a
reader does not conclude there is only one MPI path.

### `bench_guo2018_fiber_shear_cell` has no runnable command

It contains audit/contract scripts (`campaign_preflight.py`,
`evidence_contract.py`, `source_geometry_audit.py`, `status_contract.py`, …) and
a committed `data/` payload, but no `sweep.py`. There is no single command that
"runs" it and no pass/fail to record here.

---

## Observed known-fail and known-red items

Failures **as observed today**. This goal changed no crate source, so none of
these are new:

- **`examples/ci_validation.py` exits 2** before running a single sweep —
  manifest inconsistency, eight unaccounted drivers. Detailed above.
- **`bench_bond_breakage` exits 1** after 22 s. Observed tail:
  `realizations=60, bonds/run=10, max per-seed err=3.8%, KS D=0.075` /
  `empirical eps range: 0.00135 .. 0.00435; analytical median=0.00319` /
  `RESULT: CHECKS FAILED`. It sits in `ci_validation.py`'s `FULL_SWEEPS` and is
  **not** in `DOCUMENTED_EXCLUSIONS`, so the full suite cannot go green today
  even if the manifest were fixed.
- **`bench_column_collapse`** did not finish inside 400 s here, and is in any
  case a documented honest FAIL in `VALIDATION.md` (status "WITHHELD" — DIRT
  inserted only 60–72 % of requested particles while the LAMMPS comparison used
  full counts). It is explicitly excluded from the green CI suite for that
  reason.

`VALIDATION.md`'s summary table lists every other scientific bench as PASS, with
`bench_plate_sinkage` as DIAGNOSTIC. That is a **documented** status. Only the
rows in the inventory above with an observed wall-clock were re-confirmed today.

---

## Which benches are too expensive to be a routine regression gate

Measured, not guessed. Of the 43 runnable `bench_*` drivers, 36 completed and 7
hit their timeout cap. The distribution is very lopsided:

| band | count | benches |
|---|---|---|
| under 5 s | 19 | the single-contact, determinism and contract checks |
| 5–30 s | 13 | most cross-code and small-ensemble benches |
| 30–130 s | 4 | `bench_jkr_adhesion` 36 s, `bench_dmt_sjkr_cohesion` 63 s, `bench_convergence` 74 s, `bench_granular_conductivity` 121 s |
| exceeded its cap | 7 | see below |

Summed serial wall-clock of everything that completed: **538 s** (~9 min).

**Cheap enough to gate on every PR (32 benches, all ≤ 30 s).** Their summed
serial wall-clock is **244 s**, and on the `sim` queue's 8 slots they overlap to
well under a minute. This is the most useful fact here for the
migration: stage 4's "run the suite with the feature on and off and show
agreement" is affordable per-PR for the large majority of the suite, not just
nightly.

**Affordable but not free (4 benches, 30–130 s).** `bench_jkr_adhesion` and
`bench_dmt_sjkr_cohesion` each run 1,014,244 steps per case; `bench_convergence`
sweeps a resolution ladder to 1,600 particles and 150,000 steps;
`bench_granular_conductivity` is driver-configured. Keep these in a per-PR gate
only if the gate is allowed a couple of minutes.

**Too expensive to be a routine gate (7 benches).** Each exceeded its cap, so
these are lower bounds, not measurements of their true cost:

| bench | cap hit | why it is heavy |
|---|---|---|
| `bench_angle_of_repose` | > 400 s | 27 cases × 1,200 particles × 100,000 steps |
| `bench_column_collapse` | > 400 s | 33 cases × up to 1,100 particles × 80,000 steps (and a documented FAIL) |
| `bench_lebc_shear` | > 400 s | 17 cases × up to 1,634 particles, Lees-Edwards shear to steady state |
| `bench_plate_sinkage` | > 240 s | 4 cases × 2,400 particles × 45,000 steps |
| `bench_clump_haff_cooling` | > 240 s | multisphere ensemble, driver-configured |
| `bench_rod_haff_cooling` | > 240 s | multisphere ensemble, driver-configured |
| `bench_haff_ensemble` | > 240 s | many-seed ensemble, driver-configured |

These seven are exactly the many-particle collective benches — which is also to
say they are the ones a GPU port most needs to be checked against, and the ones
whose numbers are not in this file. **That is the main gap in this baseline.**
Capturing their asserted physics numbers needs a longer-running job than this
goal had, and is proposed as a followup.

Two further exclusions that are structural rather than about cost:

1. **`bench_guo2018_fiber_shear_cell`** — cannot gate anything; no runnable
   driver.
2. **`bench_mpi_decomposition`** (and `bond_mpi_drift`) — cheap (3 s observed)
   but needs the MPI feature set and multiple ranks, so it needs its own build.
   Already excluded from the stock no-MPI gate by documented decision; keep it a
   separate deliberately-invoked job.
3. **`bench_bond_breakage`** — cannot gate until it passes. A gate that is
   always red teaches reviewers to ignore it.

The thing that would make this suite unaffordable as a gate is *doubling* it:
stage 4 needs every bench run twice, feature on and feature off, and the
`cubecl-cpu` runtime may be nowhere near the scalar CPU path's speed. Whether
the doubled suite still fits in a PR gate cannot be determined until a kernel
exists. Do not infer it from these numbers.

---

## Caveats a future reader must not skip

- **The wall-clocks are not a performance model.** The box was simultaneously
  running five other agent sessions and other queued jobs, and up to 8 cost
  probes ran concurrently on the `sim` queue. The fleet throttle read `idle` at
  submission, but contention was not controlled and every number is a single
  sample, not a mean of repeats. Treat them as order-of-magnitude cost — good
  enough for "can this be a gate?", useless for "did the port get faster?".
- **A speedup claim against these numbers is invalid.** There is no GPU here.
  These are CPU numbers on a CPU-only code path. `docs/CUBECL-MIGRATION.md`
  forbids performance claims made on this box, and this file does not create an
  exception.
- **The physics numbers are the durable part.** COR, contact duration, peak
  overlap, fitted decelerations, spin-down rates, LAMMPS deltas — those are what
  a port must reproduce, and they are quoted verbatim above with the tolerance
  printed beside each. If a CubeCL path moves any of them beyond its stated
  tolerance, that is a physics change, and per the migration rules the response
  is to stop and report, not to widen the tolerance.
- **Only five benches have their numbers recorded.** The other 38 have a
  wall-clock and an exit status and nothing else. A port that changes a number
  in `bench_oblique_impact` will be caught by that bench's own gate, but this
  file will not tell you by how much.
- **`grass` is not pinned.** `Cargo.lock` is gitignored and `grass_*` resolves
  from `branch = "main"`, so re-running this baseline later may silently build
  against a different `grass`. The rev this baseline actually built against is
  in the provenance table. Pinning `grass` the way `soil_*` is pinned is a
  separate change, proposed as a followup, not done here.
- **Nothing in this goal touched crate source.** The diff is this file only. No
  `cubecl` dependency was added; no kernel was written; no physics was altered.
