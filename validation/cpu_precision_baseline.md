# CPU precision-validation baseline

Deterministic fingerprint of each example output under each host-storage precision.
`signature` = sum of |numeric cells| in the output CSV; `Delta vs double` is
the relative difference of that signature from the double-precision run.
Mixed/single store positions as f32, so they bound what the f32 GPU should reproduce.
Raw completed outputs and non-OK status JSON files are archived under `validation/results/`.

Default contact run: `python3 validation/precision_baseline.py`.
Bulk/steady-state long run: `python3 validation/precision_baseline.py --set bulk --timeout 1200`.
Combined archive regeneration: `python3 validation/_summarize.py`.

![CPU precision deltas](plots/cpu_precision_deltas.png)

*Mixed/single precision signature deltas relative to double. Completed runs are plotted;
non-OK runs are listed below
instead of being silently dropped.*

## Contact physics

| example | double signature | mixed Delta vs double | single Delta vs double | rows |
|---|---|---|---|---|
| bench_hertz_rebound | 2502.860677 | 1.94e-08 | 1.93e-08 | 1 |
| bench_oblique_impact | 2794.960146 | 2.38e-04 | 2.38e-04 | 1 |
| bench_rolling_decay | 26206.57965 | 1.60e-05 | 1.59e-05 | 8536 |
| bench_sliding_friction | 4431192.354 | 6.23e-04 | 6.23e-04 | 45000 |
| bench_sphere_haff_cooling | 122150655.7 | 4.05e-10 | 3.61e-10 | 350 |
| bench_clump_haff_cooling | 809101972.7 | 3.66e-10 | 3.25e-10 | 900 |
| bench_rod_haff_cooling | 122150306.4 | 3.01e-10 | 3.16e-09 | 350 |
| bench_jkr_adhesion | 2500.013391 | 2.88e-13 | 3.16e-13 | 1 |

## Bulk and steady state

| example | double signature | mixed Delta vs double | single Delta vs double | rows |
|---|---|---|---|---|
| bench_angle_of_repose | 64.72983668 | 4.26e-03 | 7.88e-05 | 1200 |
| bench_column_collapse | 16.93652741 | 4.24e-02 | 4.07e-02 | 341 |
| bench_hopper_beverloo | 49596.75134 | 3.23e-04 | 9.92e-03 | 45 |
| bench_granular_conductivity | timeout | timeout | timeout | - |
| bench_fiber_crossover | 8985030.094 | 1.32e-11 | 1.37e-11 | 600 |
| bench_lebc_shear | 20642475.4 | 2.48e-04 | 7.37e-05 | 140 |
| bench_plate_sinkage | 592210.5737 | 7.58e-02 | 1.89e-03 | 12057 |

## Non-OK statuses

- `bench_granular_conductivity` `precision-double`: timeout after 1200 s.
- `bench_granular_conductivity` `precision-mixed`: timeout after 1200 s.
- `bench_granular_conductivity` `precision-single`: timeout after 1200 s.
