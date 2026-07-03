# bench_restart_determinism — restart continuity + run-to-run determinism

**Gated PASS/FAIL.** Proves two properties scientists must trust before they
believe a long or interrupted production run:

- **(a) Restart continuity.** A run checkpointed at step *N* via `[restart]` and
  resumed in a **fresh process** reproduces the uninterrupted run to `<= 1e-9`
  relative on **positions and velocities** at the final step. (In the gated
  configuration the reproduction is in fact **bit-exact** — see below.)
- **(b) Run-to-run determinism.** Two independent, identically-configured
  single-rank runs are **bit-identical** (byte-for-byte equal dump files).

## The system

~300 frictional glass spheres in a fully periodic 20 mm box, given a **seeded**
Gaussian velocity field, colliding and cooling (`restitution < 1`). It is
deliberately dense and frictional: at any instant many contacts are live and
carry **tangential spring history** in the per-contact `ContactHistoryStore`, so
a restart that dropped any per-atom or per-contact state would visibly diverge.
The chaos of a granular gas makes it an unforgiving reproduction target.

Everything is deterministic single-rank: seeded insertion + seeded velocities,
deterministic Hertz–Mindlin contact and velocity-Verlet integration, and
**deterministic neighboring** (`[neighbor] every = 1, sort_every = 0`). The
neighbor list, its rebuild schedule, and the cache-locality atom re-sort are not
part of a restart file; if their *timing* differed between the resumed and
uninterrupted runs, the (physics-identical) pair-summation **order** would
differ and floating-point associativity would seed a ~1e-16 perturbation that
the chaotic gas amplifies. Rebuilding every step and disabling the re-sort
removes that confound **without changing the physics**, isolating exactly what a
restart must preserve: the per-atom + per-contact state. (With atom sorting on,
continuity still holds — but only to the ~1e-9 FP-associativity floor, not
bit-exactly.)

## Protocol (four processes)

| run | steps | role |
|-----|-------|------|
| A  | `0 .. TOTAL`            | uninterrupted reference |
| B1 | `0 .. CKPT`, `save_at_end` | write a checkpoint restart at step CKPT |
| B2 | `[restart] read`, `CKPT .. TOTAL` | resume in a fresh process |
| C  | `0 .. TOTAL`            | independent twin of A |

`sweep.py` composes the four configs from the shared physics prefix in
`config.toml` (everything after the `SWEEP CONTROL` marker is generated
per-run), runs each in its own process, and compares final binary dumps.

## Checks (all must pass)

1. `checkpoint written` — a restart file was produced at the checkpoint.
2. `restart read on resume` — B2 genuinely loaded it (`loaded N atoms from step N`).
3. `dynamics non-trivial` — the gas actually moved and stayed finite (the test isn't vacuous).
4. `continuity positions` — A vs B, per-atom by tag, `<= 1e-9` relative.
5. `continuity velocities` — A vs B, per-atom by tag, `<= 1e-9` relative.
6. `determinism (final bytes)` — `sha256(A final dump) == sha256(C final dump)`.
7. `determinism (initial bytes)` — seeded ICs are identical byte-for-byte.

## Run

```bash
python3 examples/bench_restart_determinism/sweep.py     # exits 0 iff ALL CHECKS PASSED
# or, via the harness:
~/projects/automation/bin/run-bench.sh examples/bench_restart_determinism
```

## What this bench found

On first authoring, checks (a) **failed** (continuity ~1e-4 m / ~0.1 m/s) while
determinism (b) passed bit-exactly. Bisecting by friction isolated the cause to
the **tangential contact history**: `soil_core::AtomDataRegistry::unpack_all_from_restart`
calls the **push-based** `AtomData::unpack` on stores that setup had already
sized to `nlocal` (via `truncate_all`), so restored records landed at indices
`[nlocal .. 2·nlocal)` and every atom read a stale empty slot — silently
dropping restored rotational (`DemAtom`) and per-contact tangential state. With
`friction = 0` (no tangential history and no spin) restart was already
bit-exact; only history-carrying state was lost.

**Fix (soil):** clear each store (`truncate(0)`) before the push-based unpack so
the restored records stay index-aligned with their atoms. With that fix this
bench is bit-exact on continuity. See the companion soil PR on branch
`auto/dirt-restart-determinism-bench`.

> **Cross-repo dependency:** the continuity checks pass only against a soil that
> contains the `unpack_all_from_restart` fix. Against unfixed soil, checks 4–5
> fail by design (that is the regression this bench guards). Complements the
> soil-tier `soil-restart-roundtrip-test`.
