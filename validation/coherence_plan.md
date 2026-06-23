# Design plan: scheduler-mediated host↔device coherence

Status: **approved for implementation** (2026-06-23). This is the spec a `/goal`
run should follow. It continues the GPU/MPI residency roadmap (steps 1–5 complete).

## Problem

In a GPU-resident config the device owns the canonical trajectory between window
syncs; the host `Atom`/`DemAtom` arrays are a mirror refreshed each window
(`crates/dirt_granular/src/gpu_resident.rs:145-152`). The sync is one-directional:
device→host every window, host→device never (for owned atoms). Consequences today:

- A CPU system that **reads** `Atom` *before* the resident step in a tick sees
  stale (one-window-old) data.
- A CPU system that **writes** `Atom` is **silently dropped** — the next tick
  continues from device buffers and the following download clobbers it. No crash,
  no error, just wrong physics.
- The plugin `provides()/requires()` guard only catches conflicting *plugins*, not
  a hand-written system that pokes `Atom` directly.

## Goal

The scheduler owns all host↔device transfers. Any CPU system can be added to a
resident config without silent loss: transfers are inserted automatically based on
each system's declared access, **counted**, and a **warning naming the offending
system** is emitted when a transfer forces residency loss
(*"system `X` forced a device→host sync — residency lost this tick"*). Nothing
crashes regardless of what you add.

## Non-goals

- GPU kernels are unchanged.
- Contact history (Mindlin tangential springs), the seed force, and the on-device
  cell list are **not** mirrored — they are derived device-only state with no
  representation in `Atom`. See Decision 1.
- No multi-GPU peer-to-peer transfer.

## Locked decisions

1. **Writer policy = (A)** — a host write to a `DeviceDirty` mirror triggers
   upload + **re-prime**, with a loud warning that the contact history was reset
   (a small physics discontinuity). We may switch to **(B)** (forbid writers in
   resident mode / debug-panic) later; the policy lives behind one enum so the
   swap is cheap. We are explicitly *not* doing (C) (mirroring history).
2. **Granularity = per-resource** for the MVP (whole `Atom`, whole `DemAtom`). A
   vel-only reader still forces a full download. Per-field (`pos`/`vel`/`force`)
   tracking is a later optimization, out of scope here.
3. **Rollout = feature flag `gpu_coherence`, default-off** until Phase 4 is green,
   then flip the default. Phase 1 (access metadata) is harmless and lands on
   `main` immediately, independent of the flag.

## Design refinement (2026-06-23, during implementation)

`DemAtom` (which carries the mirrored `omega`) is **not** a top-level scheduler
resource — it lives inside `AtomDataRegistry`, accessed via
`Res<AtomDataRegistry>` + `expect_mut` (interior mutability through a *shared*
ref). The access tracker therefore cannot tell a registry **read** from a
registry **write** — exactly the Phase 0 hazard. Making it honest would mean
converting ~10 `Res<AtomDataRegistry>` writers to `ResMut` across stable physics
crates (contact, wall, rotation, bond, clump, insert) — none of which even
*coexist* with the GPU mirror (the resident plugin replaces the host force /
integrator / rotation). That churn is high-risk and mostly moot.

Resolution: **track `Atom` as the single coherence trigger.** Verified `Atom` is
a plain SoA struct with no interior mutability → it is honestly classified
(`Res<Atom>` = read, `ResMut<Atom>` = write). The bridge syncs **pos + vel
(`Atom`) and omega (`DemAtom`) atomically** on any `Atom` access, so omega rides
along. `AtomDataRegistry` is *not* a separate trigger, so no `Res→ResMut` churn.

This works because pos/vel is the universal currency: essentially every host
system that touches particle state takes `Res<Atom>`/`ResMut<Atom>`. **Limitation
(documented, MVP-acceptable):** a host system that reads/writes *only*
registry-stored mirrored data (e.g. omega) without touching `Atom` won't trigger
a sync; such a system must also take `Res`/`ResMut<Atom>` or call the explicit
`mark_host_dirty` escape hatch. Phase 0 thus reduces to *verifying* Atom is
honest (it is) rather than a churn.

## Coherence model

Mirrored state: `Atom` (pos/vel, the trigger) + `DemAtom` omega (slaved to Atom
syncs via the bridge). State machine per mirror:

- `Coherent` — host and device agree.
- `DeviceDirty` — device authoritative (resident stepper just advanced).
- `HostDirty` — a CPU system wrote the host copy.

Scheduler-driven transitions, applied around each system's `run`:

| Event | Action |
|---|---|
| Resident stepper finishes | mark `DeviceDirty` (no eager download) |
| System **reads** a `DeviceDirty` mirror | download first → `Coherent`; `forced_syncs += 1`; warn with system name |
| System **writes** a mirror | download-if-needed before run; after run → `HostDirty` |
| Resident stepper enters with `HostDirty` | upload + re-prime → `forced_reprimes += 1`; warn "history reset" (policy A) |

CPU-only runs register no mirror, so the per-system check is a single `is_empty()`
branch — zero overhead.

## Integration seams (verified in `grass_scheduler/src/lib.rs`)

- `SystemParam` trait (line 275); `resource_type_id()` (292) and `is_optional()`
  (298) are the existing per-param metadata seam to extend.
- `Res::retrieve` does `borrow()` (305); `ResMut::retrieve` does `borrow_mut()`
  (326); `Option<Res>` (452), `Option<ResMut>` (476), `Local` (415). The
  read/write distinction exists here but is **not** surfaced to the scheduler yet.
- `System::prepare` (161) aggregates per-param `resource_type_id` into
  `self.indices`; `System::run` (146) / `name()` (179).
- Run loop calls `entry.system.run(resources)` at lines 709, 1000, 1897, 1930 —
  the hook points for the pre-run coherence check.

## Implementation phases

### Phase 0 — honest access declarations (prerequisite; `dirt`, maybe `soil`)

The tracker infers read/write from `Res` vs `ResMut`. But `contact_force_core`
takes `Res<AtomDataRegistry>` and mutates through it via `expect_mut` (interior
mutability through a shared ref — `crates/dirt_granular/src/contact.rs`). That
would be mis-classified as a *reader* and skip required uploads, reintroducing the
exact silent bug. **Audit every system that mutates a tracked resource through a
shared ref** and either switch to `ResMut` or add an explicit write declaration.
Must land before Phase 2 changes behavior.

Deliverable: a list of offending systems + the fix; tests still green.

### Phase 1 — access metadata in the scheduler (pure addition, no behavior change)

- `enum AccessKind { None, Read, Write }`.
- `SystemParam::access_kind() -> AccessKind` (default `None`); impl: `Res` /
  `Option<Res>` → `Read`; `ResMut` / `Option<ResMut>` → `Write`; `Local` → `None`.
- Extend `System::prepare` to build `accesses: Vec<(usize, AccessKind)>` next to
  `indices`; add `System::accesses() -> &[(usize, AccessKind)]` (default `&[]`).

Deliverable + test: a sample multi-param system reports the correct
`(resource_index, kind)` set. Ships harmlessly on `main`.

### Phase 2 — coherence registry + bridge trait (`grass_scheduler`, wgpu-free)

- `trait MirrorBridge { fn download(&self, res: &[RefCell<Box<dyn Any>>]);
  fn upload(&self, res: &[RefCell<Box<dyn Any>>]); }` — the bridge holds the
  resource indices it needs and borrows those cells itself, so wgpu never enters
  the scheduler crate.
- `CoherenceRegistry` resource: per-mirror `{ resource_index, state, syncs,
  reprimes, bridge: Box<dyn MirrorBridge> }`, plus a `register_mirror` API.
- Run-loop hook: before `entry.system.run`, call
  `registry.ensure_coherent(system.accesses(), system.name(), resources)`; after
  run, mark written mirrors `HostDirty`. The check runs **before** `retrieve`, so
  the resource cells are free — no aliasing with the system's own borrows.
- Warning is rate-limited (once per offending system per phase) and honors
  `SIM_SUPPRESS_WARNINGS`.

Deliverable + test: registry state transitions (Coherent→DeviceDirty→sync→Coherent;
write→HostDirty→reprime) with a fake in-memory bridge; counters correct.

### Phase 3 — GPU bridge + resident integration (`dirt_gpu` / `dirt_granular`)

- `impl MirrorBridge` for the resident GPU: `download` = `GpuState::download_*`
  → `Atom`/`DemAtom`; `upload` = upload host slices + re-prime.
- Register the mirror in `GpuGranularResidentPlugin::build` and the MPI variant.
- `gpu_granular_resident_step`: replace the unconditional download with `mark
  DeviceDirty`; add an upload-on-`HostDirty` path at entry (re-prime, policy A).
  Old behavior preserved when `gpu_coherence` is off.

### Phase 4 — validation

- **Reader:** resident gas + a CPU diagnostic reading `Atom` *before* Force →
  assert it sees current (not stale) data; `forced_syncs` +1; warning names it.
- **Writer:** resident gas + a CPU velocity-rescale → assert it now actually moves
  the trajectory (was silently dropped); `forced_reprimes` +1; cross-check against
  a pure-CPU thermostat run for physical correctness.
- **Cost:** CPU-only run → zero forced syncs, no warnings, no measurable overhead.
- **Bit-exactness:** resident + coherence with *no* CPU systems == resident
  without coherence (download just moves from eager to lazy-at-I/O).

## Rollout order

1. Phase 1 metadata → `main` immediately (harmless).
2. Phase 0 audit/fix → `main`.
3. Phases 2–3 behind `gpu_coherence` (default off).
4. Phase 4 green → flip `gpu_coherence` default on.

## Open risk to watch

The bridge's `upload` + re-prime (policy A) must reproduce a valid resident state
from host arrays alone. The MPI resident step
(`crates/dirt_granular/src/gpu_resident_mpi.rs`) already does ghost-slice upload +
`run_steps` re-prime every tick and is bit-exact at `window=1`, so it is the
template for the upload path. Confirm history behavior under a *local* (not just
ghost) re-upload during Phase 3.
