# GPU bonds — implementation status & validation (gpu_bonds_plan.md)

Branch `feat/gpu-bonds`. Phase A (persistent-bond primitive) and Phase B core
(elastic beam law) are implemented and validated on GPU; remaining beam features
(plasticity, Weibull thresholds, full ECS parity harness) and Phase C (MPI /
load-balance) are itemized below.

## Done

### Phase A — persistent-bond primitive (`dirt_gpu/src/bond_force.rs`, commit 4157839)
- `BondForce` `GpuForce` hook over a persistent CSR bond topology (`BondTopology`:
  per-point `offsets`/`partner`/`r0`), mirroring the CPU `BondStore`'s per-local-atom
  `Vec<Vec<BondEntry>>`.
- **In-place per-bond state** (no ping-pong — topology is fixed), **i-centric
  atomic-free** accumulation (each endpoint computes the mirror half), **flag-and-skip
  breakage**.
- Phase-A law: central bond `f = k_n·(L−r0)` (= bond-based peridynamics), critical
  |strain| breakage.
- `BondTopology::from_bond_store(&BondStore, &Atom)` — builds the GPU CSR from a CPU
  bond store, resolving `partner_tag → index` (drops MPI-cut partners; offsets span
  all GPU atoms, ghosts own no bonds).
- **Tests (GPU):** intact bond restores (symmetric, momentum-conserved); past-threshold
  bond breaks on both copies together (consistent f32 decision).

### Phase B core — elastic beam (`dirt_gpu/src/beam_bond_force.rs`, commit a8a19b2)
- `BeamBondForce`: normal `E·A/L`, history shear `G·A/L` (reproject ⊥ n̂ + integrate),
  twist `G·J/L`, bending `E·I/L` moments; per-channel critical damping; constant
  beam-stress breakage (`σ = Fₙ/A + 2|M_bend|r_b/J`, `τ = |F_t|/A + |M_tor|r_b/J`).
- Verified the **mirror-image algebra term by term** against the CPU single-owner
  `bond_force`: each endpoint's own half reproduces CPU's `+f/i −f/j`, `+M/i −M/j`,
  shear-torque-same-on-both — so it stays atomic-free with no canonical frame.
- **Tests (GPU):**
  - normal restoring (symmetric, linear momentum conserved);
  - bending transmits moment to the partner;
  - **force + torque parity vs the f64 CPU beam formula** under a combined
    normal+shear+bending+rotation load — `< 2e-3` relative (f32 vs f64).

All 7 `dirt_gpu` tests pass (`--features precision-double`).

## Benchmark finding — host-authoritative GPU is a dead end (2026-06-24)

LEBC contact+neighbor timing, `bench_lebc_shear` (1634 glass beads, 42k steps,
`DIRT_FORCE=cpu` vs `gpu`, identical config — harness committed):

| Contact force | Wall-clock |
|---|---|
| CPU Hertz-Mindlin | **8.6 s** |
| GPU host-authoritative (`GpuGranularForcePlugin`) | **156 s** |

**GPU is 18× *slower*.** `gpu_granular_force` re-uploads state, rebuilds the cell
list, and does **two blocking device-waits** to download force+torque *every step*
— ≈3.7 ms/step of pure sync/launch latency vs the CPU's ~0.2 ms/step total. At 1634
particles the M5 is also underfilled. **Decision: do not pursue the host-authoritative
path (it round-trips every step); full residency is the only route to a GPU speedup.**
Target scale is BPM aggregates (20–30 spheres/shape, ≤~10k particles total).

### Toward resident periodic + Lees–Edwards (started)
- **Periodic minimum-image (+ LE tilt) in the bond kernel** — `BeamBondConfig.{lx,ly,lz,tilt_xy}`;
  the bond vector uses the triclinic minimum image, so a BPM aggregate spanning a
  periodic boundary stays bonded. Bonds don't use the cell list, so this needs no
  cell-list change. Test: a bond across a periodic x-boundary matches the equivalent
  close pair (`< 1e-3`). LE velocity-offset on the damping term of a y-wrapped bond
  is not yet applied (orthogonal-periodic is fully correct).

## Remaining (toward full CPU parity)

1. **End-to-end ECS parity harness** — run the actual `dirt_bond::bond_force` system
   and `BeamBondForce` on the same multi-bond scene (chain / lattice) and compare
   forces & trajectory. The formula-parity test covers the per-bond math; this would
   close the loop against the real CPU code path.
2. **Plasticity** — axial + bending return-maps (`plasticity::update_axial/bending`):
   extra per-bond state (`eps_p_axial`, `theta_p_bend`, `*_max`) + the return-map in
   WGSL.
3. **Per-bond Weibull thresholds + breakage-criterion variants** — currently a single
   constant σ_max/τ_max in params; the CPU samples per-bond thresholds (seeded) and
   supports several criteria. Add a per-bond `thresholds` buffer (host-sampled).
4. **Contact exclusion on GPU** — the contact kernel must skip bonded pairs
   (parity with `bonds.are_excluded`); check point `i`'s bond slots for `partner == j`.
5. **Granular + beam torque composition** — `BeamBondForce` currently owns the aux-rate
   (torque) buffer; running it alongside `GranularForce` needs a torque-seed kernel in
   `soil_gpu` so both can accumulate (`+=`).
6. **Phase C — MPI / load balance** — per-step `partner_tag→index` re-resolution on
   migration (host-resolve + upload), and per-bond state migration on rebalance
   (download/migrate/re-upload, or host-authoritative through a rebalance).
7. **f32 break-timing** — decide tolerance vs compensated-f32 stress accumulation for
   the parity bar (plan open question 3).
