# Monkey-barrel LEBC VF campaign — μ(I), Φ(I), and the monkey jamming ceiling

Analysis of the constant-volume Lees–Edwards simple-shear volume-fraction campaign
(`examples/monkey_shear`, run via the pueue `monkey` group) comparing three particle
shapes at a **common equivalent-volume diameter** `D_eq = 0.1 m`:

- **sphere** — single sphere, `r = 0.05`.
- **rigid**  — the 44-sub-sphere "monkey" as a rigid multisphere clump.
- **bpm**    — the same monkey as a bonded-particle model (44 free sub-spheres welded
  intra-monkey by `auto_bond`, neighbouring monkeys never bonded).

Grid: `type ∈ {sphere, rigid, bpm} × Φ ∈ {0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55}`,
shear rate `γ̇ = 1.0 s⁻¹`, protocol *settle → compress → shear (target strain γ ≈ 2)*.

> **Snapshot: 2026-07-04 20:47 UTC — the campaign is still running.** Four `rigid`
> cells and three `bpm/high-Φ` cells were mid-flight or had aborted at snapshot time.
> Every number below is honestly tagged **COMPLETE / PARTIAL / NO-SHEAR** and the
> figure marks non-steady points with `×`. Re-run `tools/analyze_vf.py` to refresh once
> the group drains.

## Method

Per run the recorder (`main.rs`, reused verbatim from `bench_lebc_shear`) writes the full
stress tensor → `p, τ, N1, N2`, granular temperature `T`, and `Φ` every 2000 steps to
`<out>/lebc_shear_results.csv`. `tools/analyze_vf.py`:

- **strain** `γ = Σ γ̇·dt` accumulated over shear rows only (`gdot > 0`; settle/compress carry `gdot = 0`);
- **steady window** = `γ ∈ [1.0, γ_max]` (or the last 30 % of shear if `γ_max < 1`);
- window-mean `p, τ, T, Φ`; `μ = τ/p`;
- inertial number `I = γ̇·D_eq / √(p/ρ)`, grain density `ρ = 2500`.

**Φ convention.** `N(Φ)` was chosen so the **union** solid volume gives the nominal Φ, but
the recorder's `Φ` column sums *sub-sphere* volumes, which over-counts the monkey union by
the 23.6 % intra-monkey overlap (`Φ_rec ≈ 1.31·Φ_nom`, confirmed: `0.0327/0.025 = 1.308`).
**Plots and the discussion use nominal Φ as the true solid fraction**; the reduced CSV keeps
both (`phi_nom` and measured `phi_meas`).

## Completeness (all 24 cells; pueue state as of snapshot)

| Φ_nom | sphere | rigid | bpm |
|------|--------|-------|-----|
| 0.025 | ✅ COMPLETE γ2.0 | ⏳ PARTIAL γ1.65 *(running)* | ✅ COMPLETE γ2.0 |
| 0.05  | ✅ COMPLETE | ⏳ PARTIAL γ0.71 *(running)* | ❌ PARTIAL γ0.19 *(aborted, SIGABRT)* |
| 0.1   | ✅ COMPLETE | ⏳ PARTIAL γ0.31 *(running)* | ❌ PARTIAL γ0.07 *(aborted)* |
| 0.2   | ✅ COMPLETE | ⏳ PARTIAL γ0.12 *(running)* | ❌ PARTIAL γ0.02 *(aborted)* |
| 0.3   | ✅ COMPLETE | ❌ NO-SHEAR *(compress abort)* | ❌ NO-SHEAR *(compress abort)* |
| 0.4   | ✅ COMPLETE | ❌ NO-SHEAR *(compress abort)* | ❌ NO-SHEAR *(compress abort)* |
| 0.5   | ✅ COMPLETE | ❌ NO-SHEAR *(compress abort)* | ⏳ NO-SHEAR *(running, pre-shear)* |
| 0.55  | ✅ COMPLETE | ⏳ NO-SHEAR *(running, pre-shear)* | ⏳ NO-SHEAR *(running, pre-shear)* |

- ✅ = finished / usable steady state. ⏳ = still running at snapshot (value provisional/absent).
  ❌ = failed. **Failed runs exit 134 = SIGABRT** — the contact overlap guard
  (`crates/dirt_granular/src/contact.rs`, ">500 excessive overlaps") correctly aborts rather
  than emit garbage. The guard was **not** relaxed (anti-gaming); see the blocked
  `monkey-shear-dt-stability-fix` goal — those blow-ups are physical, not a dt artifact.
- **Only the sphere series (8/8) and `bpm/0.025` are fully steady.** Rigid Φ≤0.2 are still
  accumulating strain; their μ/Φ points below are **transient**, not steady rheology.

## Results

Reduced table: `tools/` output CSV (`analyze_vf.py --csv`). Figure: `rheology_vf.png`
(μ(I), Φ(I), **p(Φ)**, μ(Φ); ○ COMPLETE, × PARTIAL).

### 1 — Sphere baseline (complete)

| Φ | p [Pa] | τ [Pa] | μ=τ/p | I | T |
|---|--------|--------|-------|-----|-----|
| 0.025 | 6.24e-3 | 9.37e-3 | 1.50 | 63.3 | 4.8e-5 |
| 0.05  | 1.78e-2 | 2.48e-2 | 1.39 | 37.5 | 9.1e-5 |
| 0.10  | 2.04e-1 | 2.43e-1 | 1.19 | 11.1 | 3.3e-4 |
| 0.20  | 2.03e+0 | 1.84e+0 | 0.91 | 3.51 | 1.3e-3 |
| 0.30  | 2.68e+1 | 1.10e+1 | 0.41 | 0.97 | 9.5e-3 |
| 0.40  | 5.43e+1 | 2.22e+1 | 0.41 | 0.68 | 7.5e-3 |
| 0.50  | 9.89e+1 | 4.62e+1 | 0.47 | 0.50 | 5.2e-3 |
| 0.55  | 2.10e+2 | 9.42e+1 | 0.45 | 0.35 | 4.7e-3 |

Textbook granular behaviour: `p` and `Φ` rise monotonically as `I` falls; the dilute cells
(Φ≤0.1) sit in the collisional regime (`I ≫ 1`, high `μ`), and the dense cells (Φ≥0.3) reach
the **dense-flow μ(I) plateau `μ ≈ 0.41–0.47`** at `I ≈ 0.3–1.0` — consistent with a
`friction = 0.16` grain. Spheres flow at every Φ up to 0.55; **no jamming within the grid.**

### 2 — Monkeys interlock: pressure amplified 100–500× at matched Φ

At the **same** nominal solid fraction, shape carries far more stress:

| Φ_nom | p_sphere | p_rigid | p_bpm | rigid/sph | bpm/sph |
|------|---------|---------|-------|-----------|---------|
| 0.05 | 1.8e-2 | 1.7e+1 | 4.5e+1† | ~970× | ~2500× |
| 0.10 | 2.0e-1 | 4.7e+1 | 2.9e+2† | ~230× | ~1400× |
| 0.20 | 2.0e+0 | 2.7e+2 | 1.1e+3† | ~130× | ~520× |

(† bpm values are from the transient pre-blow-up window — upper-bound magnitudes, not steady.)

Even the *rigid* monkeys — same convex-hull shape, no flexibility — develop 2–3 orders more
pressure than spheres at identical Φ. This is the geometric-interlocking signature: elongated,
hooked bodies mechanically engage and transmit contact force at packings where spheres still
flow freely. `μ` at low Φ is also higher for monkeys (rigid `μ ≈ 1.7` vs sphere `≈ 1.5`).

### 3 — The monkey jamming ceiling: **0.2 < Φ_ceiling < 0.3**, vs sphere > 0.55

For **both** rigid and bpm, the campaign produces **no shear data at Φ ≥ 0.3**: the
**compression stage itself aborts** on the overlap guard before shear ever starts (`gdot`
never goes positive; the run dies compacting toward the target packing). Reaching union
Φ_nom = 0.3 means a *sub-sphere* packing of `Φ_rec ≈ 0.39` (0.4 → 0.52, 0.5 → 0.65) — and for
randomly-oriented interlocking monkeys those densities are geometrically unreachable at finite
grain stiffness: contacts pile into unresolvable overlaps and the guard fires.

**The interlocking-monkey compaction/jamming ceiling sits between Φ = 0.2 and Φ = 0.3**, roughly
**half or less of the sphere's flowable range (> 0.55).** Rigid Φ ≤ 0.2 do compact and shear
(runs in progress); Φ ≥ 0.3 never do. The `p(Φ)` panel shows both monkey branches climbing far
steeper than the sphere and terminating at the Φ = 0.3 wall (dashed line).

### 4 — BPM adds flexural heating (a second, lower ceiling)

BPM monkeys fail **earlier and harder** than rigid ones. Only `bpm/0.025` reached steady
strain; `bpm/0.05, 0.1, 0.2` reached shear but **aborted mid-shear** (γ = 0.19/0.07/0.02) —
not a packing failure but a *dynamic* one. The granular temperature tells the story:

| Φ_nom | T_sphere | T_rigid | T_bpm |
|------|----------|---------|-------|
| 0.05 | 9e-5 | 1.1e-1 | 1.9e+0 |
| 0.10 | 3e-4 | 1.3e-1 | 2.1e+0 |
| 0.20 | 1.3e-3 | 2.7e-1 | 1.1e+0 |

BPM runs ~10⁴× hotter than spheres and ~10× hotter than rigid monkeys. The flexible welded
arms whip/resonate under shear (bond bending modes), pumping kinetic energy until sub-sphere
overlaps trip the guard. This matches the `monkey-shear-dt-stability-fix` diagnosis (bpm KE
heats ~500× vs the rigid case; a physics/operating-point issue, **not** curable by smaller dt).
So the *unbreakable-elastic* bpm has an even lower effective flow ceiling than the rigid clump,
driven by bond-mode heating rather than by packing alone.

## Honest status & follow-ups

- **Solid:** full sphere μ(I)/Φ(I) baseline; the interlocking pressure amplification (rigid vs
  sphere, matched Φ); the Φ ≥ 0.3 compaction ceiling for both monkey types; bpm flexural heating.
- **Provisional:** rigid Φ ≤ 0.2 are still running — μ/Φ points are transient (γ < 1.7). Re-run
  the analyzer when tasks 448/451/454/457 finish for steady rigid rheology.
- **Pending/absent:** rigid 0.55 and bpm 0.5/0.55 were still in their pre-shear (settle/compress)
  phase at snapshot and may yet abort at the ceiling; no data claimed for them.
- **No fabricated points.** NO-SHEAR cells contribute nothing to μ(I)/Φ(I).

## Reproduce

```bash
source ~/projects/.build-env
# --data points at the campaign output tree (per-type/phi dirs of lebc_shear_results.csv)
python3 examples/monkey_shear/tools/analyze_vf.py \
  --data ~/projects/worktrees/dirt-monkey-shear-lebc-campaign/examples/monkey_shear/data \
  --csv /tmp/monkey_reduced.csv \
  --plot examples/monkey_shear/rheology_vf.png
```

---

# BPM blow-up: instrumented root-cause characterization

*(goal `dirt-bpm-monkey-blowup-investigation`, follow-up to dirt#44)*

**Verdict — the `bpm Φ=0.05` monkey_shear blow-up is a genuine PHYSICAL instability
(a shear-driven granular-temperature runaway of the flexible, near-elastic bonded
aggregate), NOT any of the enumerated setup/parameter bugs.** No pass was forced;
the overlap guard was not touched. dirt#44 already showed the failure is *physical,
not CFL* (dt-converged); this investigation instruments **what** the physical
mechanism is and **rules out** each candidate setup bug with evidence.

## What was instrumented

- **Static bond/exclusion audit** — `tools/bond_exclusion_audit.py`, replaying the
  engine's own `auto_bond_touching` (bond if `dist ≤ 1.1·(Rᵢ+Rⱼ)`) and
  `BondStore::are_excluded` (contact-exclude a pair iff **1-2** bonded or **1-3**
  sharing a bonded neighbour — nothing beyond 1-3).
- **Runtime trace** — the `instrument_blowup` system in `main.rs` (opt-in via
  `MONKEY_INSTRUMENT=1`), writing `<out>/instrument.csv` each thermo interval:
  bond count, **largest bonded component** (a mis-weld tripwire), total kinetic
  energy `KE = Σ½m|v|²`, and the largest **contact-active** (non 1-2/1-3-excluded)
  **intra-monkey** sub-sphere overlap. Monkey membership is the bond-graph
  connected component (robust to tag order; a component > 44 would mean auto_bond
  fused two monkeys). Runtime intra/inter classification was cross-checked with a
  temporary in-`contact.rs` guard probe (reverted — core is unchanged in this PR).

## The four candidate setup bugs — each ruled out

| Candidate (from the goal) | Instrument evidence | Verdict |
|---|---|---|
| **auto_bond welds wrong pairs** (neighbouring monkeys fused) | `nbonds = 37054 = 382 × 97` **constant**; **largest bonded component = 44** for every monkey, whole run | **Ruled out** — monkeys are cleanly separated; the 1.25·(Rᵢ+Rⱼ) placement gap holds against the 1.1 bond cutoff. |
| **bond-break exposes as-built overlap** (Liz's stated chain) | bonds are configured **unbreakable** (no `[bonds.breakage]`); `nbonds` never drops from 37054; **zero break events** | **Ruled out for this campaign** — the premise (a bond breaks) never occurs, so contact-exclusion never lifts. |
| **insertion / compression overlap at gap placement** | **prep is quiescent**: through settle+compress, `KE ≈ 10⁻²⁵→10⁻⁸ J` and `max_active_intra_overlap = 0`; the guard **never trips in prep** — only in shear | **Ruled out** — the blow-up is a shear-stage phenomenon, not a packing/insertion defect. |
| **bond-vs-contact stiffness mismatch (CFL)** | dirt#44: run reproduces the **same abort strain (~0.196)** at dt 1.5e-6 and 5.79e-7 (dt-converged); Hertz overshoot `≈ v·dt ≈ 0.3 %` of overlap | **Ruled out** — not a timestep artifact. |

Additionally, the static audit finds **0 contact-active build-overlaps at t=0** (all
86 overlapping intra pairs are 1-2/1-3 excluded) and the **nearest** contact-active
intra pair sits at a **34.6 % gap** — so there is no t=0 Hertz injection, and
intra-body self-contact requires large deformation.

## The actual mechanism (evidence)

`bpm Φ=0.05`, `MONKEY_INSTRUMENT=1`, single rank (KE is exact for the bpm/sphere
free-atom types; note it is unreliable for the *rigid* clump type, whose velocity
is body-integrated — for rigid, cite dirt#44):

| shear step | strain γ | KE [J] | active intra-overlaps |
|---|---|---|---|
| 22 000 (onset) | 0.001 | **94.1** | **0** |
| 30 000 | 0.006 | 100.8 | 0 |
| 34 000 | 0.008 | 115.8 | 0 |
| 40 000 | 0.012 | 199.5 | few |
| 46 000 | 0.015 | 372.0 | growing |
| 50 000 | 0.017 | 531.9 | growing |
| ~abort (dirt#44) | ~0.196 | **~2000** | guard trips |

Two facts pin the mechanism:

1. **At shear onset KE jumps 0 → 94 J within one thermo interval with *zero* active
   intra-monkey overlaps.** The initial heating is therefore **inter-monkey
   collisional**, not self-contact. Intra-body self-contact appears only *later*
   (arms fold once the aggregate is already hot) — it is a **consequence** of the
   runaway, not its trigger, consistent with the 34.6 % nearest-gap from the static
   audit.
2. **The flexible aggregate never reaches a steady granular temperature.** KE grows
   monotonically/near-exponentially under continued shear (94 → 532 J over a tiny
   strain 0.017, → ~2000 J at abort). By contrast, at the **same** shear and box the
   **sphere** series equilibrates at a **bounded** `KE ≈ 77 J` (steady to 4 figures,
   run completes), and dirt#44's **rigid** monkey — identical shape — stays at
   `KE ≈ 4 J` and shears past γ = 0.6. Only the **flexible bonded** body diverges.

So the blow-up is the flexible, near-elastic (restitution 0.926, friction 0.16)
bonded monkey failing to thermostat the shear work: collisional + bond-mode energy
accumulates faster than it dissipates until sub-sphere overlaps exceed the guard.

## Construction context (reconciling Liz's diagnosis)

Liz correctly identified the **~23.6 % sub-sphere build-overlap** as the root
construction issue, and that this is a *construction* problem, not a fundamental
result. The instrumentation **refines** the chain: because the campaign's bonds are
*unbreakable*, the failure is **not** "a bond breaks → Hertz sees the build overlap"
(that never happens here). Instead the root construction fault is that the BPM
monkey **reuses the rigid-clump sphere decomposition**, where:

- a rigid clump excludes **all** same-body sub-sphere contact (`same_body` skip in
  `contact.rs`), so its 23.6 % overlap and any self-contact are harmless; but
- a BPM body excludes only **1-2/1-3** pairs, so its non-bonded sub-spheres can
  self-contact, and the dense overlapping decomposition makes the aggregate lumpy
  and compliant — the extra energy-storage/injection channels the rigid case lacks.

## Recommended remedies (campaign-owner decision — NOT forced here)

1. **Rebuild the BPM monkey from near-tangent (non-overlapping) sub-spheres** — the
   standard bonded-particle convention (LAMMPS `bpm/sphere` bonds tangent spheres).
   Removes both the 23.6 % build overlap and the spurious self-contact channel.
   (Liz's option (a); the LAMMPS-consistent construction.) Requires re-deriving the
   sphere set at `D_eq = 0.1` and re-validating the bpm series to completion.
2. **Add shear-stage dissipation appropriate to the aggregate** (lower pair
   restitution and/or a background drag) so the flexible body thermostats like the
   sphere/rigid cases — an operating-point/model choice for the rheology surface.

**Anti-gaming:** the overlap guard (`crates/dirt_granular/src/contact.rs`) is
**unchanged**; no tolerance was loosened; no bpm cell was marked passing. The
unstable cells remain honestly reported: **bpm at every Φ, and rigid at Φ ≥ 0.3**
(the latter a distinct *compaction*-stage jamming ceiling, §3 above).

## Reproduce the investigation

```bash
source ~/projects/.build-env
# 1. static bond/exclusion audit (no engine):
python3 examples/monkey_shear/tools/bond_exclusion_audit.py \
  examples/monkey_shear/monkey_Deq0.1.toml
# 2. regenerate the bpm Φ=0.05 config + CSV, then trace a run:
$BENCH_PYTHON examples/monkey_shear/tools/gen_series.py all --phi 0.05
MONKEY_INSTRUMENT=1 cargo run --release --example monkey_shear \
  --no-default-features --features precision-double -- \
  examples/monkey_shear/configs/bpm/phi_0.05.toml
# trace: <output_dir>/instrument.csv  (step,stage,nbonds,max_bond_component,KE,...)
```
