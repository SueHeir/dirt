# Hopper quiescence benchmark

Prototype of a DEM optimization aimed at quiescent/coherent regions, to
evaluate whether it is worth implementing as a LAMMPS fix. Everything lives
in this example — the DIRT crates are unmodified. The stock Hertz-Mindlin
contact plugin is replaced by a local copy (`contact.rs`) with optimization
hooks, plus a region state machine (`quiescence.rs`).

## The optimization: region velocity coherence

A coarse cell grid (`cell_size`, ~4 diameters) computes, every `check_every`
steps, each cell's mean velocity and the max deviation of any member from that
mean. Cells coherent for `quiet_checks` consecutive evaluations become:

- **plug** (|v̄| ≥ `v_sleep_tol`): members keep integrating, but their internal
  contacts *force-freeze* — once a pair's force is steady (relative change <
  `freeze_rel_tol`) the cached force and torques are applied directly each step
  and the Hertz-Mindlin evaluation, material lookups, and spring update are
  skipped. The contact *thaws* when the pair's relative position drifts more
  than `unfreeze_disp_frac × min(r1, r2)` from its freeze-time separation. The
  displacement tolerance sets the force error floor (Hertz: dF/F ≈ 1.5·dδ/δ).
- **sleeping** (|v̄| < `v_sleep_tol`): member particles stop integrating
  (velocity/force zeroed, like the stock `pin` fix) and pairs where *both*
  particles sleep are skipped entirely.

Wake-up: plug cells demote the moment any member's velocity deviates from the
stored cell mean (checked every step). Sleeping particles snapshot the partial
force they receive from gravity, walls, and awake neighbors; a relative
deviation > `wake_force_rel` wakes the cell — this is what propagates the
wake-front when the blocker is removed or when overburden grows during
filling.

## Scenario

Stage 1 ("filling", 120k steps): particles rain in through rate insertion
(1000 every 1500 steps, 20000 total) and settle on a blocker wall covering the
funnel outlet. Stage 2 ("flowing", 200k steps): the blocker is removed on the
first step (fixed step count → identical staging across variants); the bed
discharges into the deep catch basin below the outlet and re-settles there.

The column is tall (ceiling z = 0.45) so the settled bed (~z = 0.18) clears the
insertion zone, and the floor sits at z = −0.15 so the discharged pile (~0.08 m
deep) never backs up into the outlet at z = 0.06.

10× gravity (like the stock hopper example) accelerates settling; E = 10 MPa
(softened, standard DEM practice) allows dt = 2 µs.

## Running

```bash
cargo build --release --no-default-features --example hopper_quiescence
for v in baseline coherence; do
  ./target/release/examples/hopper_quiescence \
      examples/hopper_quiescence/config_$v.toml
done
```

Each run writes `config_<variant>_stats.csv`:
`step, elapsed_s, ke_total, n_atoms, n_discharged, n_asleep, n_plug, pairs,
frozen_pairs, skipped_pairs, top_z` (`frozen_pairs` counts plug-frozen
contacts; `top_z` is the tallest particle, a settled-bed height proxy).

`bench.sh` runs both variants and reports: **(A)** wall time split by phase
(filling vs emptying), **(B)** the settled fill height (max z at the
end-of-filling step) baseline vs coherence, and **(C)** the emptying-speed
discharge curve plus the step at which 90% has discharged — (B) and (C) are
the fidelity checks that the optimization does not change the physics.

## Important: tangential damping sign fix

The local contact copy deliberately diverges from stock `dem_granular` in one
line: the tangential dashpot is `+γₜ·vₜ` instead of `−γₜ·vₜ`. With this
codebase's convention (`vr = v_j − v_i`, `+ft` applied to particle `i`), the
stock sign *injects* energy during sliding (P = +γₜ|vₜ|²), which keeps dense
frictional beds fluidized forever — the stock hopper example never actually
reaches its KE settling threshold because of this. The single-particle
rebound benchmark passes because the *normal* dashpot sign is correct.
The same fix should be upstreamed to `dem_granular` (fused + standalone
tangential) before any quiescence work lands there. `dem_wall` is unaffected
(plane walls have no tangential friction).

To reproduce the stock behavior for comparison, set `USE_STOCK_CONTACT=1` in
the environment — the example then registers the stock
`HertzMindlinContactPlugin` instead of the local copy (the optimization and the
sign fix both inactive; CSV pair/freeze columns read zero).

## Measured results (M-series Mac, single rank, 5000 particles, 225k steps)

**A. Performance** — wall time by phase (`bench.sh` section A):

| variant            | fill [s] | empty [s] | total [s] | speedup |
|--------------------|---------:|----------:|----------:|--------:|
| baseline           |     34.2 |      84.5 |     118.7 |   1.00× |
| coherence          |     27.6 |      53.9 |      81.5 |   1.46× |

The win is larger in the empty phase (1.57×) than fill (1.24×): the settled bed
above the blocker and the re-piled discharge are mostly quiescent, while early
filling still has particles raining in. At the end of filling the coherent
variant has ~97% of particles asleep and skips ~94% of contact pairs; the
re-piled discharge sleeps fully (KE = 0, 97% pairs skipped) by the end.

**B + C. Validation** (the optimization must not change the physics):

- **Fill height** — tallest particle at the end of filling: baseline 0.10832 m
  vs coherence 0.10906 m, agreeing to 0.7 mm (below one particle radius).
- **Emptying speed** — discharge curves agree to within insertion-RNG noise
  (~1%); both variants reach 90% discharged at step 109000 and fully empty
  (5000/5000) by step 120000.

Key findings:

- **Region coherence delivers 1.46×** end-to-end by skipping sleeping pairs
  entirely (one branch on two mode bytes — no distance, no lookup, no model)
  and skipping integration for sleepers. The plug-freeze path caches forces in
  co-moving regions but, in this benchmark, the win is dominated by sleeping:
  active discharge is not coherent enough to form long-lived plugs.
- The speedup bound is the quiescent fraction of (particles × steps): a
  benchmark dominated by active flow gains nothing; a large settled bed with
  a slow local process would gain far more than 1.46×.
- Sleep/wake thresholds can stabilize marginal structures (legitimately frozen
  arches, or numerically over-stable ones; runs are not seed-deterministic), so
  threshold sensitivity deserves a sweep before trusting the method on
  arching-sensitive problems.

## Interpreting results for LAMMPS

- The wins concentrate where beds are static or co-moving: end of filling,
  post-discharge pile. During active discharge nothing freezes (correct
  behavior — no fidelity loss, but also no speedup).
- The optimization skips integration and *all* pair work in sleeping regions
  and would pair naturally with neighbor-list rebuild exemptions (not
  prototyped here — rebuild triggers still fire for coherent moving regions, a
  further available win in LAMMPS).
- The plug-freeze fast path maps onto LAMMPS `pair granular`'s history array
  (a cached force per contact), but on this problem it contributes little over
  sleeping; whether it is worth the per-pair bookkeeping depends on how much of
  the target workload is co-moving rather than static.
- The prototype is single-rank; an MPI version needs ghost exchange of the
  sleep/plug flags and rank-agreement on cell states.
</content>
</invoke>
