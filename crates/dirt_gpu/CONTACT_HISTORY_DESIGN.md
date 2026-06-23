# GPU Contact History + Tangential Friction + Rotation — Design (Milestone #1)

Status: **design for review** (implementation follows once approved).

## Goal

Add Mindlin tangential friction with **persistent per-contact spring history** and
**rotational dynamics** to the resident GPU DEM stepper (`dirt_gpu`), matching
`dirt_granular::hertz_mindlin_contact_force` within f32 tolerance. Today the GPU
force is normal-only Hertz and translation-only — frictionless, so no realistic
granular behavior. This milestone makes it real DEM.

**In scope (v1):** tangential Mindlin friction + spring history, Coulomb cap,
torque, rotational velocity-Verlet (omega).
**Deferred:** rolling/twisting springs (history slots reserved but not computed),
cohesion/adhesion, walls (#2), large-coordinate frame (#3).

## The CPU model (what we must reproduce)

Per `tangential.rs` + `contact.rs`:
- Per-atom variable-length list of contacts: `(partner_tag, spring[7], active)`.
- Spring stored **canonically** (lower-tag frame); a `sign = ±1` converts to the
  local (i,j) frame each step. `spring[0..3]`=tangential, `[3..6]`=rolling, `[6]`=twist.
- Per step, per pair (half-list, processed once):
  1. `s = sign * stored`; remove normal component: `s -= (s·n) n`.
  2. integrate: `s += v_t * dt`, where `v_t = v_rel - (v_rel·n) n` and
     `v_rel = vc_j - vc_i`, `vc = vel + omega × (±r n)`.
  3. spring Coulomb cap: if `k_t|s| > μ|F_n|`, scale `s`.
  4. `F_t = k_t s + γ_t v_t`, `γ_t = 2 β √(5/6) √(k_t m_r)`; cap `|F_t| ≤ μ|F_n|`.
  5. torque `τ_i = (r1 n) × F_t`, `τ_j = (r2 n) × F_t`; force `+F_t` on i, `−F_t` on j.
  6. write back updated canonical spring; mark active. Prune inactive after the loop.
- `k_t = 8 G* √(R* δ)`.

## Key GPU design decisions

### 1. i-centric full kernel ⇒ no atomics, each atom owns its half of each spring
Our force kernel is **full** (thread `i` scans all neighbors `j`, accumulates force
on `i` only; the pair is evaluated again by thread `j`). So each atom independently
stores **its own** spring for each partner, in its **local (i→j) frame** (no
canonical/tag/sign needed). The two halves stay exact mirror images:
`n_ji=−n_ij`, `v_rel` flips, so `v_t_ji=−v_t_ij` ⇒ `spring_ji=−spring_ij` and
`F_t_on_j=−F_t_on_i` — provably identical to the CPU's single canonical store, but
**with zero atomics and zero cross-thread writes** (thread `i` only writes atom `i`'s
slots). This is the central simplification.

### 2. Per-atom fixed-capacity slots, ping-pong buffers (automatic prune)
Each atom has `MAX_CONTACTS` slots. Each step the kernel **reads old slots, writes
new slots** into a second buffer, then we swap (ping-pong). Contacts not seen this
step simply aren't written ⇒ pruned for free; new contacts start with spring 0.
- `MAX_CONTACTS` = 32 (dense sphere packing contacts ≈ 12; 32 is safe headroom;
  tunable). Overflow ⇒ drop + a debug counter (should never trigger in practice).

### 3. Key by persistent atom index (valid on single device)
The GPU cell list builds a *sort index* (`sorted_atoms`) but does **not** permute the
atom arrays, and there's no MPI migration on a single device — so **atom index `i`
is stable across steps**. We key contact slots by partner **index** (not tag), and
search is a short linear scan of `MAX_CONTACTS`.
> Assumption to revisit if we ever add GPU spatial sorting or multi-GPU migration:
> then switch to tag keying + an index remap after each reorder.

### 4. Rotation
Add resident `omega` (angular velocity), `torque`, and `inv_inertia` (from
`DemAtom`). The force kernel computes contact velocity with `omega × r n` (needed for
`v_t`), accumulates `torque[i]` (i-centric, no atomics), and we add rotational
velocity-Verlet half-kicks (`omega += 0.5 dt · inv_inertia · torque`). Orientation
quaternion is **not** needed for force on spheres → skip in v1 (add later for output).

## Data structures (new GPU buffers)

| Buffer | Type | Size | Notes |
|--------|------|------|-------|
| `omega` | f32 | 3N | resident, read_write |
| `torque` | f32 | 3N | written by force kernel (i-centric overwrite) |
| `inv_inertia` | f32 | N | uploaded once (DemAtom: 1/(0.4 m r²)) |
| `contact_partner[2]` | u32 | N·MAX_CONTACTS ×2 | ping-pong; `SENTINEL=0xFFFFFFFF`=empty |
| `contact_spring[2]` | f32 | 3·N·MAX_CONTACTS ×2 | tangential spring (reserve 7/slot if adding rolling/twist) |

Memory at N=216k, MC=32: partner ~28MB×2, spring(3) ~83MB×2 ≈ 220MB. Fine on M5
unified memory; MC=16 halves it. Spring stored as 3/slot in v1 (not 7) to save memory;
widen to 7 when rolling/twist land.

## Kernel changes

1. **`hooke_force` → full contact kernel** (rename to `contact_force`): per neighbor
   `j` in contact, compute Hertz normal (as now) **plus** Mindlin tangential using the
   old spring looked up from this atom's old slots (search by partner index; 0 if new).
   Steps 1–5 of the CPU model, in the local frame (no sign). Accumulate `force[i]` and
   `torque[i]`. Append `(j, new_spring)` to the new-slot list. After the neighbor scan,
   flush new slots to `contact_partner_new[i]`/`contact_spring_new[i]`, padding with
   `SENTINEL`.
2. **`integrate_initial`/`integrate_final`**: extend to also half-kick `omega` with
   `inv_inertia·torque`. (Translation unchanged.)
3. Host swaps `contact_*` old/new each step (just swaps bind groups or buffer refs).

### Bind-group note
We're already near Metal's storage-buffer limit; adding omega/torque/inv_inertia +
2×(partner,spring) is ~7 more storage buffers. We request `adapter.limits()` already,
which on M5 is far above the default 8, so this is fine — but verify at bring-up.

## Step sequence (resident, per step)
```
integrate_initial (translate + rotate)
clear_cells → assign → prefix_sum → scatter        (cell list)
contact_force (normal + tangential + torque + spring write-back to NEW slots)
integrate_final (translate + rotate)
swap contact buffers (old ↔ new)
```
All in one submit; zero per-step host transfers (unchanged from current residency).

## Params additions
Add `g_eff: f32` (effective shear modulus) and `mu: f32` (friction coefficient).
Single material in v1 (scalars), like `e_eff`/`beta` today.

## Validation plan
1. **Single-step vs CPU**: one GPU step vs a CPU velocity-Verlet step using a CPU
   reference of the *same* tangential model (or dirt's `hertz_mindlin` directly) on a
   packing with shearing contacts (nonzero relative tangential velocity + spin) →
   match force, torque, omega within f32 tol.
2. **History persistence**: run K steps; a sustained tangential contact must build up
   spring force and hit the Coulomb cap, matching CPU within tol over the first several
   steps (before chaotic divergence).
3. **Stability**: 500+ resident steps stay finite/bounded.
4. **Physics check**: a friction-sensitive scenario (two spheres in sustained sliding
   contact, or — ideally — angle-of-repose-style) behaves qualitatively right
   (frictional, not frictionless).
5. **Regression**: with `mu=0` and no spin, results reduce to the current
   normal-only Hertz (sanity).

## Risks / open questions
- **`MAX_CONTACTS` overflow** in very dense/polydisperse packings → add a device
  overflow counter; tune MC.
- **Exact spring-rotate match**: the "remove normal component then integrate" order
  and the two Coulomb caps (spring cap *and* force cap) must match CPU exactly — main
  source of subtle drift if mis-ordered.
- **Memory** at large N with 7-slot springs; start at 3-slot (tangential only).
- **Index-keying assumption** (decision #3) — fine now, revisit if GPU sort/migration
  is added.
- **Interaction with #3 (coordinate frame)**: independent — history stores springs,
  not positions; #3 only changes how the kernel reads positions. Do #3 after.

## Out of scope (later milestones)
Rolling/twisting springs (slots reserved), cohesion/adhesion, polydisperse cell
efficiency, walls + gravity (**#2**), large-coordinate frame (**#3**), quaternion
orientation, multi-material tables.
