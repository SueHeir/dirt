# GPU validation vs CPU

**Adapter:** Apple M5 Pro (Metal 4), via wgpu. f32-only (no f64), so the GPU is
effectively **single** precision (f32 storage + f32 accumulation).

**Verdict: PASS.** The GPU contact-force kernels reproduce dirt's real CPU force
code to ~1e-6 relative — far inside the 5e-3 tolerance, and even tighter than the
~1e-3 normally expected for f32-vs-f64. The previously-"not tested" code runs
correctly with no NaNs, errors, or kernel bugs.

## Tier 1 — per-evaluation force/torque vs REAL dirt CPU (`compare_cpu`)

GPU force-hook stack (soil `GpuState` + dirt `GranularForce`/`WallForce`) vs
`dirt_granular::hertz_mindlin_contact_force` + `dirt_wall::wall_contact_force`,
single force evaluation (no integration), f32 GPU vs f64 CPU. Pass = max
component-wise relative diff < 5e-3.

| config | max rel force | max rel torque |
|---|---|---|
| (a) head-on normal | 8.90e-7 | 0 |
| (b) overlap + sliding (tangential) | 8.38e-7 | 3.07e-7 |
| (c) overlap + spin (ω×r + torque) | 8.38e-7 | 2.87e-7 |
| (d) floor wall + slide + spin | 1.02e-6 | 4.28e-7 |

All ~1e-6 — about 10× f32 epsilon, i.e. the kernel arithmetic is correct; the
residual is pure f32 rounding of one evaluation. Normal, Mindlin tangential,
spin (ω×r lever + torque), and planar-wall paths all validated.

## `pile` — full resident GPU sim (qualitative)

2744 grains dropped under gravity into a box, stepped entirely on-device
(resident GpuState: integrate + rotation aux-DOF + cell list + Hertz/Mindlin +
walls). Physically correct:
- KE peaks on impact (5.47e-1) then decays to 5.96e-7 = **1.1e-6 of peak** —
  friction + damping dissipate energy as they should.
- min particle z = 0.0493 (r=0.05) — grains rest **on** the floor, **no
  tunneling**.
- pile height settles to a stable 1.275 m.

## Performance (Tier 3, full step CPU single-thread vs GPU resident)

| N | CPU ms/step | GPU ms/step | speedup |
|---|---|---|---|
| 8,000 | 0.967 | 0.144 | 6.7× |
| 64,000 | 8.238 | 0.655 | 12.6× |
| 216,000 | 28.086 | 1.778 | 15.8× |

Speedup grows with N (better GPU occupancy) — 15.8× at 216k grains vs
single-threaded CPU.

## Which CPU precision does the GPU match?

Per *evaluation*, the GPU matches even CPU-**double** to ~1e-6 (Tier 1 compares
to f64). Over a full *trajectory*, the GPU (f32) will track CPU-**single**, since
both share f32 storage; the f32-vs-f64 trajectory divergence is the
single-vs-double band already recorded in `cpu_precision_baseline.md`
(~1e-4–1e-3 for friction-heavy cases, ~1e-8 for normal/energy cases).

## Gap / next step (not done here)

Tier 1 validates kernel correctness on a *single force evaluation*. It does **not**
yet run a full GPU *trajectory* on the same configs as the CPU baseline examples
(e.g. drive a GPU `bench_hertz_rebound` drop and compare measured COR against
`validation/results/bench_hertz_rebound__precision-single.csv`). That end-to-end
trajectory check — driving each baseline example through the resident GpuState
and diffing the recorded metric within its Δ-band — is the remaining work toward
"fully validated GPU code." It needs a small GPU runner per example scenario; not
straightforward enough to include here.
