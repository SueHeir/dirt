# Coherence validation (coherence_plan.md Phase 4)

> **Current main-branch status (2026-07-06):** this page records validation from
> the historical GPU/coherence branch. Current `main` has no `gpu_coherence`
> feature, resident GPU bridge, `GpuGranularResidentPlugin`, or
> `GpuGranularResidentMpiPlugin`. The results below are branch provenance, not a
> claim that main currently includes these GPU-resident coherence pieces.
>
> **Provenance:** historical GPU/coherence branch `origin/gpu`, with the coherence
> validation recorded in dirt commit `077f748` after the bridge work in `d59b9f8`.
> Main removed the GPU/coherence artifacts in `c631d97`; the validation command
> below is intentionally archival and is not runnable from current `main`.

Scheduler-mediated host↔device coherence: a CPU system added to a GPU-resident
config now syncs transparently (and is attributed/counted) instead of silently
dropping. On the historical branch this was implemented across three repos,
behind the `gpu_coherence` feature.

## What landed

| Layer | Change | Repo / commit |
|---|---|---|
| Access metadata | `SystemParam::access_kind()` → `System::accesses()` (Read/Write per resource) | grass `6bbb88a` |
| Coherence engine | `CoherenceRegistry` + `MirrorBridge` (3-state machine; run-loop hooks; attributed warning) | grass `448621b`, `ffc1ec5` |
| GPU bridge | `ResidentMirrorBridge` + coherence-aware resident step (lazy pull, reupload+reprime) | dirt `d59b9f8` |
| Design | `Atom`-only trigger (omega slaved); no registry `Res→ResMut` churn | dirt `5baf9a6` |

## How it works

The resident step marks its `Atom` mirror `DeviceDirty` after each window instead
of eagerly downloading. The scheduler, before running each system, reads that
system's declared access and — if it reads a `DeviceDirty` mirror — pulls the
device copy back to the host via the bridge (counted; warned once per system).
A host **write** marks the mirror `HostDirty`; the next resident window re-uploads
the host state and re-primes (policy A — contact history resets). The resident step
is *self-managed* (`Option<ResMut<CoherenceRegistry>>`) so its own `Atom` access
doesn't trip the auto hooks. With the feature off, no registry is registered and
the resident step uses the original eager per-window download — unchanged.

## Validation: `coherence_validate` example

Historical command on `origin/gpu` only:

```bash
cargo run -p dirt_granular --example coherence_validate --release \
   --no-default-features --features precision-double,gpu_coherence
```

Four resident runs of the same 216-grain drop (window 50 × 20 ticks, Apple M5 Pro):

| Run | coherence | CPU writer | CPU reader |
|---|---|---|---|
| eager_base | off | – | – |
| coh_base | on | – | each tick |
| coh_writer | on | vel ×0.9/tick | each tick |
| eager_writer | off | vel ×0.9/tick | – |

**Results (all PASS):**

```
A  max|coh_base - eager_base|     : 0.000e0    (lazy pull == eager download — bit-faithful)
B  forced device->host syncs      : 20         (one per reader tick, attributed)
C  max|coh_writer - coh_base|      : 1.681e-1   (CPU writer RESPECTED under coherence)
D  max|eager_writer - eager_base|  : 0.000e0    (CPU writer SILENTLY DROPPED without it)
```

**C vs D is the headline:** the identical CPU velocity-damp writer changes the
trajectory under coherence (1.7e-1) but is dropped bit-for-bit by the eager
resident path (0.0). Attribution fires by name:

```
[coherence] system `coherence_validate::reader_system` forced a device→host sync — residency lost this tick
[coherence] resident GPU re-primed from host-modified state (contact history reset)
```

- **A (bit-faithful):** the lazy device→host pull reproduces the eager download exactly.
- **B (attributed + counted):** every host read is a counted, named sync.
- **Cost:** syncs scale with host *consumers*, not windows — a config with no
  per-tick reader pays zero per-window downloads (vs eager's download every window).

Unit coverage: grass `coherence` module (4 tests: device-dirty pull, coherent
no-op, write→HostDirty + take clears, self-managed skip) + access metadata (2
tests); dirt_granular 33 tests green with `gpu_coherence` on and off.

## Limitations (documented, MVP-acceptable)

- **Omega-only consumers:** a host system that touches *only* registry-stored
  mirrored data (omega) without `Atom` won't trigger a sync — it must also take
  `Res`/`ResMut<Atom>`. pos/vel is the universal currency, so this is rare.
- **Direct post-loop reads:** under coherence the host `Atom` is only fresh after a
  *system* reads it. Code that reads `Atom` directly after the run loop (not via a
  system) sees the last-synced state — use a readout system at `PostFinalIntegration`
  (as the validation does), which the standard dump path already is.
- **Position-teleporting writes:** `reupload_locals` assumes positions stay within
  the existing grid bounds (true for velocity/force edits); a teleport would also
  need a `set_grid`.
- **MPI resident variant:** on the historical branch, coherence was wired for the
  single-GPU `GpuGranularResidentPlugin`; `GpuGranularResidentMpiPlugin` (which
  re-primes every tick anyway) was unchanged.
