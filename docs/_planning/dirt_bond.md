# dirt_bond planning

## Purpose

`dirt_bond` provides the Bonded Particle Model (BPM) for DIRT: it treats each
bonded pair as a solid cylindrical elastic beam resisting four independent
deformation channels (axial, shear, twist, bending), with optional per-channel
viscous damping, configurable breakage criteria (beam-stress, strain, or
linear-interaction-envelope families), and piecewise-linear plasticity on the
bending and axial channels. Bonds are created at setup (auto-bonding touching
pairs or loading from a LAMMPS data file) and suppressed from the normal
granular contact law via `soil_core`'s `BondStore`; on breakage the pair is
removed from the store and Hertz–Mindlin contact resumes automatically.

---

## Public surface to document

| Item | Kind | File:line |
|------|------|-----------|
| `DemBondPlugin` | Plugin struct | `src/lib.rs:521` |
| `BondConfig` | Deserialized `[bonds]` section | `src/lib.rs:194` |
| `BondHistoryStore` / `BondHistoryEntry` | Per-atom bond history (implements `AtomData`) | `src/lib.rs:337–436` |
| `BondBreakage` | Active criterion resource | `src/lib.rs:448` |
| `BondPlasticity` | Active plasticity resource | `src/lib.rs:483` |
| `BondMetrics` | Per-step diagnostics, thermo keys | `src/lib.rs:500` |
| `breakage::BreakageCriterion` | Breakage trait | `src/breakage.rs:256` |
| `breakage::BreakageConfig` | TOML-deserializable criterion enum | `src/breakage.rs:560` |
| `breakage::ThresholdDistribution` | `Constant` / `Weibull` / `CrackBand` | `src/breakage.rs:177` |
| `breakage::BondGeom`, `BondLoads`, `BondKinematics`, `BondThresholds` | Snapshots passed to criterion | `src/breakage.rs:102–147` |
| `breakage::BreakMode` | Which failure branch tripped | `src/breakage.rs:150` |
| `breakage::per_bond_uniform_samples` | MPI-stable RNG draw | `src/breakage.rs:69` |
| `plasticity::PlasticityConfig` | Top-level `[bonds.plasticity]` | `src/plasticity.rs:79` |
| `plasticity::BendingPlasticityConfig` | Bending variants | `src/plasticity.rs:89` |
| `plasticity::AxialPlasticityConfig` | Axial variant | `src/plasticity.rs:136` |
| `plasticity::BondPlasticityModel` | Runtime model | `src/plasticity.rs:289` |
| `plasticity::update_bending` / `update_axial` | Return-map functions | `src/plasticity.rs:376 / 432` |
| Setup systems: `auto_bond_touching`, `load_bonds_from_file`, `extend_ghost_cutoff_for_bonds`, `init_breakage`, `init_plasticity`, `init_bond_history` | Setup systems (PostSetup) | `src/lib.rs:662–933` |
| Force system: `bond_force` | Per-step update (Force schedule slot) | `src/lib.rs:937` |
| Thermo systems: `zero_bond_metrics`, `output_bond_metrics` | PreForce / PostForce | `src/lib.rs:1307–1349` |

---

## Config/TOML schema

All keys live under `[bonds]`.

### Bond geometry

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `auto_bond` | bool | `false` | Bond every pair within `bond_tolerance × (R_i + R_j)` at setup |
| `bond_tolerance` | f64 | `1.001` | Multiplier on sum-of-radii for auto-bond eligibility |
| `bond_radius_ratio` | f64 | `1.0` | Bond cylinder radius as a multiple of `min(R_i, R_j)` |
| `ghost_cutoff_multiplier` | f64 | `2.5` | Scales max bond r₀ to extend MPI ghost cutoff; set to `0.0` to disable on single-rank runs |

### Stiffness — material mode (preferred)

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `youngs_modulus` | Option\<f64\> | `None` | Young's modulus *E* (Pa); when set, derives `K_n = E·A/L` and `K_bend = E·I/L` |
| `shear_modulus` | Option\<f64\> | `None` | Shear modulus *G* (Pa); when set, derives `K_t = G·A/L` and `K_tor = G·J/L` |

### Stiffness — direct override (used when matching modulus is absent)

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `normal_stiffness` | f64 | `0.0` | `K_n` (N/m) |
| `shear_stiffness` | f64 | `0.0` | `K_t` (N/m) |
| `twist_stiffness` | f64 | `0.0` | `K_tor` (N·m/rad) |
| `bending_stiffness` | f64 | `0.0` | `K_bend` (N·m/rad) |

If both material mode and direct mode are set for the same channel, material
mode wins (`src/lib.rs:1034–1037`).

### Damping — critical-damping ratios

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `beta_normal` | f64 | `0.0` | Critical-damping ratio β for normal channel; `γ_n = 2β√(m*·K_n)` |
| `beta_shear` | f64 | `0.0` | Same for shear |
| `beta_twist` | f64 | `0.0` | Same for twist; uses reduced MOI: `γ_tor = 2β√(I*·K_tor)` |
| `beta_bending` | f64 | `0.0` | Same for bending |

### Damping — raw overrides (bypass β calculation)

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `normal_damping` | Option\<f64\> | `None` | `γ_n` (N·s/m) directly |
| `shear_damping` | Option\<f64\> | `None` | `γ_t` (N·s/m) |
| `twist_damping` | Option\<f64\> | `None` | `γ_tor` (N·m·s/rad) |
| `bending_damping` | Option\<f64\> | `None` | `γ_bend` (N·m·s/rad) |

### Breakage

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `seed` | u64 | `0` | Global seed for per-bond threshold RNG; same seed + same topology → identical breakage pattern regardless of MPI decomposition |
| `breakage` | Option\<[bonds.breakage]\> | `None` | Breakage criterion sub-table; absent → bonds never break |

`[bonds.breakage]` fields:

- `kind` (required): one of `"unbreakable"`, `"axial_force"`, `"axial_stress"`,
  `"axial_strain"`, `"combined_stress"`, `"combined_strain"`,
  `"interaction_linear_force"`, `"interaction_linear_stress"`,
  `"interaction_linear_strain"`
- Two-branch families (`axial_*`, `combined_*`): `tensile` (required) and
  `shear` (optional, absent → `∞` threshold). Each is a `ThresholdDistribution`:
  - `{ kind = "constant", value = <f64> }` — same value for every bond
  - `{ kind = "weibull", mean = <f64>, m = <f64>, l_calib = <f64>, l_min = <f64 default 0.0> }` — length-scaled 2-parameter Weibull
  - `{ kind = "crack_band", value_ref = <f64>, l_ref = <f64>, eps_yield = <f64 default 0.0>, l_min = <f64 default 0.0> }` — Bažant crack-band deterministic rescaling
- `InteractionLinear*` families: four optional channels `axial`, `shear`,
  `bending`, `twist`, each a `ThresholdDistribution`; absent channels drop out
  of the sum

### Plasticity

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `plasticity` | Option\<[bonds.plasticity]\> | `None` | Plasticity sub-table; absent → all channels elastic |

`[bonds.plasticity.bending]` fields:

- `kind`: `"guo_bending"`, `"guo_trilinear"`, or `"piecewise"`
- `yield_stress` (f64, required for `guo_bending` / `guo_trilinear`): material yield stress σ₀ (Pa)
- `"guo_trilinear"` requires `[bonds].youngs_modulus` to be set (derives ε_e = σ₀/E)
- `"piecewise"`: `breakpoint_strains` (Vec\<f64\>, extreme-fibre strain), `slope_multipliers` (Vec\<f64\>), `length_calibration` (Option\<f64\> m, crack-band regularization)

`[bonds.plasticity.axial]` fields:

- `kind`: `"piecewise"` only
- `breakpoint_strains` (Vec\<f64\>, axial strain ε = (L−L₀)/L₀)
- `slope_multipliers` (Vec\<f64\>)
- `length_calibration` (Option\<f64\> m)

### Explicit bond file

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `file` | Option\<String\> | `None` | Path to LAMMPS data file with a `Bonds` section |
| `format` | Option\<String\> | `None` | File format; only `"lammps_data"` is supported |

---

## Key behaviors, invariants, and gotchas

### Bond/contact exclusion

`DemBondPlugin::build` (re-)calls `app.add_plugins(soil_core::BondPlugin)`,
which registers `BondStore` as an `AtomData` and installs the exclusion hook
into the granular contact loop. The contract is: any atom pair present in
`BondStore::bonds` on the local rank will be **skipped by the Hertz–Mindlin
force** (src/lib.rs:598). On breakage, the pair is removed from `BondStore` on
**both** ends (`bond_force`, `src/lib.rs:1280–1303`) and normal contact
immediately resumes on the next step. Double-checking: `bond_force` only
removes bonds from the store *after* the force pass for that step, so a bond
that breaks this step still contributes forces for one final step before it
disappears.

### When auto-bonding runs

`auto_bond_touching` is registered in `ScheduleSetupSet::PostSetup` and gated by
`first_stage_only()` (`src/lib.rs:605–612`). It runs once, before neighbor
setup bakes the ghost cutoff. The O(N²) all-pairs loop over local atoms is
acceptable at setup time. Bonds wrap periodically via minimum-image on axes
flagged periodic in `Domain` (`src/lib.rs:689–698`). **MPI caveat**: at setup
time all atoms are on rank 0 (pre-migration), so auto-bonding on a multi-rank
run creates all bonds on rank 0, which then get distributed during the first
migration. This is fine as long as `ghost_cutoff` is extended first.

### Ghost cutoff extension (MPI)

`extend_ghost_cutoff_for_bonds` runs in PostSetup, ordered after both bond
creation systems and before `neighbor_setup` (`src/lib.rs:619–629`). It
finds the global max `r₀` via `all_reduce_min_f64(-local_max)` and bumps
`domain.ghost_cutoff` to `max_r0 × ghost_cutoff_multiplier` if that exceeds
the current value. The default multiplier of 2.5 covers: 1× for the bond
itself + 2× for shared-neighbour 1-3 exclusion + 0.5× safety margin for bond
stretch (`src/lib.rs:204–211`). If `bond_missing` is nonzero at run time,
raise this multiplier. Setting it to `0.0` disables the extension entirely
(safe only on single-process runs).

### Breakage is irreversible

Once a bond breaks, the pair is removed from `BondStore` and `BondHistoryStore`
on the owning rank (`src/lib.rs:1285–1303`). No mechanism re-creates a broken
bond. Breakage is detected in the `bond_force` system and the removal is
deferred to after the full force pass for that step (bonds_to_break collected,
then applied after the loop). Breakage counts accumulate in `BondMetrics` and
persist until the simulation ends (`total_bonds_broken` only grows).

### MPI bond consistency

Per-bond failure thresholds are drawn deterministically from `(min(tag_i,tag_j),
max(tag_i,tag_j), seed)` via `per_bond_uniform_samples` (`src/breakage.rs:69`).
The tag pair is canonicalized (lo, hi) and mixed via two SplitMix64 rounds before
seeding a `SmallRng`. This guarantees that two MPI ranks visiting the same bond
from different sides compute identical thresholds. Breaking: the lower-tag atom
owns the bond computation (`src/lib.rs:1002`), so on a multi-rank run only the
owning rank detects the break and removes it from its local `BondStore`. The
partner rank carries its mirror entry until the next migration/ghost-update
cycle reconverges the bond lists. This is by design: the `BondStore` is an
`AtomData` and its ghost-communication wire format handles this.

### `bond_type` is stored but not consumed

The `bond_type` field parsed from LAMMPS data files is stored in `BondEntry`
and carried through the atom data system, but the `bond_force` system ignores
it — all bonds use the same global `BondConfig` parameters. It is a reserved
slot for future per-type parameter dispatch (`src/lib.rs:216`, README.md note).

### Axial plasticity sign conventions

`update_axial` handles both tension (`eps_axial > 0`) and compression (`<0`)
symmetrically using the signed `eps_e` and a `signum` projection for the
return-map direction (`src/plasticity.rs:432–459`). The `eps_max_axial`
tracker records `|ε_axial|`, not the signed value, so kinematic hardening
freezes the cap at the largest absolute strain seen so far (`src/plasticity.rs:440`).

### Bending plasticity uses extreme-fibre strain

Breakpoint strains in `[bonds.plasticity.bending]` are always in extreme-fibre
strain `ε = r_b · θ_bend / l_b`, not in angle units. The conversion factor
`scale = l_b / r_b` is applied inside `evaluate_piecewise` (`src/plasticity.rs:183–186`).
This makes the material parameters geometry-independent.

### Crack-band regularization

Both plasticity and breakage thresholds support crack-band length rescaling.
For plasticity: the post-yield strain extents `(ε[i] − ε[0])` scale by
`l_calib / l_bond` and post-yield slopes divide by the same factor, preserving
`∫ F du` = plastic energy per bond (`src/plasticity.rs:249–285`). For breakage:
`CrackBand` rescales the threshold as `eps_yield + (value_ref − eps_yield) ×
l_ref / max(l_bond, l_min)` (`src/breakage.rs:243–248`). Use `crack_band`
thresholds with the strain-family criteria (`axial_strain`, `combined_strain`,
`interaction_linear_strain`) for mesh-invariant fracture energy.

---

## Tutorial outline: cantilever validation walkthrough

The `bond_cantilever` example (`examples/bond_cantilever/`) is the cleanest
validation target. A writeup would walk through:

1. **Setup the scene**: 10 glass-like spheres (radius 1 mm) spaced 2 mm along x;
   CSV insertion via `[[particles.insert]]`; material `"bpm"` with `E = 1 GPa`,
   `ν = 0.25`.

2. **Bond configuration**: `[bonds]` section — `auto_bond = true`,
   `bond_tolerance = 1.001`, `bond_radius_ratio = 1.0`, material mode with
   `youngs_modulus = 1.0e9` and `shear_modulus = 4.0e8`, `beta_* = 1.0` (critical
   damping on all channels). No breakage, no plasticity.

3. **Anchor via `[[freeze]]`**: The leftmost sphere (tag 0, x = 0) belongs to
   group `"anchor"`. `[[freeze]]` zeros velocity, force, omega, torque, and
   angular momentum every step at PostForce. This is essential: without zeroing
   rotation, bond shear produces a lever-arm torque that spins the anchor,
   creating positive feedback (the chain swings unstably about the anchor).

4. **Plugin chain in `main.rs`**: `CorePlugins + GranularDefaultPlugins +
   FixesPlugin + GravityPlugin + DemBondPlugin`. Note that `DemBondPlugin`
   must come after `GranularDefaultPlugins` so the contact exclusion hook has
   a contact system to hook into.

5. **Theoretical prediction**: Euler–Bernoulli beam theory, uniform load
   (self-weight per unit length `q = ρ A g`), cantilever deflection
   `δ = q L⁴ / (8 E I)`. For `E = 1 GPa`, `r_b = 1 mm`, `L = 18 mm` (9 bonds),
   `ρ = 2500 kg/m³`: predicted tip deflection ≈ -9.56 × 10⁻⁷ m.

6. **Running and checking**: `cargo run --release --example bond_cantilever
   --no-default-features -- examples/bond_cantilever/config.toml`. The recorder
   writes `data/cantilever.csv` with columns `step, t, tip_z, max_strain,
   bond_count, bond_missing`. PASS criteria: `bond_missing = 0`, tip stays
   bounded (no blowup), steady-state `tip_z ≈ −1.0 × 10⁻⁶ m` (within ~5% of
   theory given discrete-element vs. continuum approximation).

7. **Common failure modes**: non-zero `bond_missing` means `ghost_cutoff` is too
   small (raise `ghost_cutoff_multiplier`); tip blowup with wrong sign conventions
   (the README documents the historical sign bug); omitting `[[freeze]]` rotation
   zeroing (chain swings freely about the anchor).

---

## Doc gaps

1. **Breakage criterion reference table in docs**: `docs/src/physics/bonds.md`
   mentions only `combined_stress` and `Weibull`. The full criterion menu
   (`axial_force/stress/strain`, `combined_strain`, `interaction_linear_*`) is
   documented only in `src/breakage.rs` and the crate README. Should appear in
   the physics page with the same table as `src/breakage.rs:16–26`.

2. **Crack-band threshold explanation in docs**: `CrackBand` is not mentioned in
   `docs/src/physics/bonds.md` at all. This is the correct regularization for
   mesh-refinement invariance with strain criteria; it warrants its own paragraph.

3. **Plasticity variants not shown in docs**: `docs/src/physics/bonds.md` mentions
   `guo_bending` and piecewise axial, but not `guo_trilinear` or the
   `length_calibration` crack-band option on the plasticity envelopes.

4. **MPI ghost cutoff extension not documented in physics page**: The
   `ghost_cutoff_multiplier` key and why 2.5 is needed (1× bond + 2× shared-
   neighbour + 0.5× stretch) appears only in `BondConfig` rustdoc and the crate
   README, not in `docs/src/physics/bonds.md`.

5. **`bond_type` status note**: The docs page already has the "stored but unused"
   note (last sentence of Configuration section), but there is no roadmap for
   when it will be consumed or how per-type parameters would be specified.

6. **Sign convention rationale**: The "lower tag owns" convention is documented
   in the physics page but the lever-arm torque calculation (`τ_shear = (L/2) n̂
   × F_t`, same sign on both particles) and why twist/bending moments get +M on i
   and -M on j is not explained in the docs page (it is in `src/lib.rs:38–56`).

7. **No validation benchmark for breakage or plasticity in examples/**: The
   `fiber_bond_breakage` and `fiber_bond` examples have data and figures but lack
   a `sweep.py` with generate/start/graph matching the canonical example layout
   (see `dirt-example` skill). They are also not referenced from
   `docs/src/physics/bonds.md` other than a brief mention.

8. **`are_excluded` / contact exclusion internals**: The actual mechanism by
   which `BondStore` suppresses Hertz contact is in `soil_core`, not in
   `dirt_bond`. The docs page says contact is "suppressed via the substrate's
   `BondStore`" but does not explain how to verify this at run time (check that
   `bond_fiber_tensile_overlap` produces no extra force on bonded overlapping
   spheres).

---

## Suggested placement

The existing `docs/src/physics/bonds.md` is the right home. It needs to be
expanded, not restructured. Suggested additions in order:

1. A **Breakage criterion reference** section after "Breakage and plasticity",
   with the full 8-variant table from `src/breakage.rs:16–26` and a paragraph
   on threshold distributions (`constant`, `weibull`, `crack_band`) with
   crack-band regularization explained.

2. A **Plasticity variants** subsection under the plasticity paragraph, adding
   `guo_trilinear`, `piecewise` (both channels), and the `length_calibration`
   crack-band option on the envelopes.

3. A **MPI notes** section covering `ghost_cutoff_multiplier` and the
   `bond_missing` diagnostic.

4. A **Full TOML reference** table (as above) replacing or supplementing the
   current partial `[bonds]` example block.

5. Expand **Worked examples** to add `fiber_bond` (cantilever bending +
   bending plasticity) and `fiber_bond_breakage` (axial-strain / combined-
   stress / Weibull scatter), with PASS/FAIL criteria for each.
