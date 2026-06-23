# CPU precision-validation baseline

Deterministic fingerprint of each example's output under each host-storage
precision. `signature` = sum of |numeric cells| in the output CSV; `Δ vs double`
is the relative difference of that signature from the double-precision run.
Mixed/single store positions as f32, so they bound what the f32 GPU should
reproduce. Raw outputs archived under `validation/results/`.

| example | double signature | mixed Δ vs double | single Δ vs double | rows |
|---|---|---|---|---|
| bench_hertz_rebound | 2502.808371 | 1.35e-08 | 1.35e-08 | 1 |
| bench_oblique_impact | 2791.96286 | 2.27e-04 | 2.27e-04 | 1 |
| bench_rolling_decay | 26206.57965 | 2.25e-05 | 2.24e-05 | 8536 |
| bench_sliding_friction | 4431298.354 | 5.93e-04 | 5.93e-04 | 45000 |
| bench_sphere_haff_cooling | 122150655.5 | 5.24e-10 | 6.92e-10 | 350 |
| bench_clump_haff_cooling | 122150305.9 | 3.27e-09 | 4.85e-09 | 350 |
| bench_rod_haff_cooling | 122150305.9 | 1.39e-09 | 5.24e-10 | 350 |
| bench_jkr_adhesion | 2500.013391 | 2.88e-13 | 3.16e-13 | 1 |

## Key physics final-states (the clean single-collision tests)

`mixed` and `single` are nearly identical: the f32 *storage* (positions/vels)
dominates the deviation from `double`, while the Accum type (f64 vs f32) barely
matters. So these f32 results are what the always-f32 GPU should reproduce.

| test | quantity | double | mixed | single |
|---|---|---|---|---|
| hertz_rebound | COR | 0.90166214 | 0.90164524 | 0.90164524 |
| oblique_impact | vt_rebound | 2.44111917 | 2.44308873 | 2.44308873 |
| jkr_adhesion | f_pulloff | 5.8904862e-3 | 5.8904862e-3 | 5.8904863e-3 |

## Using this for GPU validation

The GPU kernels are intrinsically f32. To validate them, run the same example
config on the GPU and compare its output against the **mixed/single** rows here
(via `validation/results/` or `cpu_precision_final_states.json`), using each
example's `Δ vs double` as the order-of-magnitude tolerance:

- hertz/jkr/haff-cooling: agree to ~1e-8 or tighter — GPU must match closely.
- friction/tangential (oblique, sliding, rolling): ~1e-4–1e-3 spread — looser band.

The `double` run is the theory anchor (gold reference); `mixed`/`single` are the
f32 reference the GPU is held to. Regenerate anytime: `python3
validation/precision_baseline.py`.

## Scope & deferred examples

This baseline covers the **contact-physics** benchmarks — single/few-body
collisions and granular-gas cooling — which directly exercise the normal and
tangential contact force the GPU kernels compute, and which run in seconds.

The **bulk/steady-state** benchmarks (`angle_of_repose`, `column_collapse`,
`hopper_beverloo`, `granular_conductivity`, `fiber_crossover`, `lebc_shear`,
`plate_sinkage`) are emergent-behaviour validations that take ~6–10 min per run,
so ~1–2 h across 3 precisions. They are deferred from the inline baseline; run
them on a longer budget with:

    python3 validation/precision_baseline.py bench_angle_of_repose ...

They append to `validation/results/` and merge into the summary via
`python3 validation/_summarize.py`.
