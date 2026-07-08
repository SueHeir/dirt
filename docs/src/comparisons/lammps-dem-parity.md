# LAMMPS vs DIRT — DEM Feature Parity

This page tabulates the granular / DEM contact features implemented in
[LAMMPS](https://docs.lammps.org/) against those implemented in DIRT, marking
what each side **has** and **lacks**, with source citations for **both** codes.
It is a living reference used to seed follow-up implementation and validation
goals.

## Scope & method

- **LAMMPS reference** is the read-only clone at `~/projects/reference/lammps`.
  The authoritative prose is `doc/src/pair_granular.rst` and
  `doc/src/fix_wall_gran.rst`; the implementation lives in `src/GRANULAR/`
  (sub-model classes registered via the `GranSubModStyle(...)` macros in the
  `gran_sub_mod_*.h` headers), with bonded and rigid particles in `src/BPM/`
  and `src/RIGID/`. LAMMPS's modern granular pipeline is `pair_style granular`
  driving a `GranularModel` (`src/GRANULAR/granular_model.cpp`) that composes
  pluggable **normal / damping / tangential / rolling / twisting / heat**
  sub-models.
- **DIRT** paths are relative to the `dirt` crate root. DIRT selects contact
  sub-models by string keys read from `dirt_atom::MaterialTable`
  (`contact_model`, `adhesion_model`, `rolling_model`, `twisting_model`) rather
  than a Rust enum; the granular crate branches on those strings
  (`crates/dirt_granular/src/contact.rs:77, 254, 541, 606`).

Every row cites an exact file (and, where useful, a line or symbol). **No
feature is listed as present without a source on that side.** Absences were
verified by exhaustive grep and are annotated with what was searched.

### Legend

| Mark | Meaning |
|------|---------|
| ✅ | Implemented |
| ⚠️ | Implemented but **restricted** (subset of geometries/paths, or a simplified form) |
| ❌ | Not implemented (searched, not found) |

---

## Normal contact models

| Feature | LAMMPS | DIRT |
|---|---|---|
| **Hooke** (linear spring `F=k_n δ`) | ✅ `hooke` — `pair_granular.rst:102`; `GranSubModNormalHooke` `src/GRANULAR/gran_sub_mod_normal.cpp:160` | ✅ `contact_model="hooke"` — `crates/dirt_granular/src/contact.rs:722` (`hooke_contact_force`, normal at 813–847) |
| **Hertz** (`F ∝ R_eff^½ δ^{3/2}`) | ✅ `hertz` — `pair_granular.rst:117`; `gran_sub_mod_normal.cpp:188` | ✅ default — `crates/dirt_granular/src/contact.rs:190` (`contact_force_core`, assembly 391–430) |
| **Hertz/material** (from `E, ν`) | ✅ `hertz/material` — `pair_granular.rst:128`; `gran_sub_mod_normal.cpp:197` | ✅ DIRT's Hertz path already takes `E*, ν` from the `MaterialTable` (effective modulus at `crates/dirt_granular/src/contact.rs` normal assembly) — no separate keyword |
| **DMT** cohesion (Hertz − `4πγR`) | ✅ `dmt` — `pair_granular.rst:141`; `GranSubModNormalDMT` `gran_sub_mod_normal.cpp:288` | ⚠️ `adhesion_model="dmt"`, **Hertz path only** — `crates/dirt_granular/src/contact.rs:254, 392` (limitation `crates/dirt_granular/src/lib.rs:35`) |
| **JKR** adhesion (contact-radius solve, tensile to pull-off) | ✅ `jkr` — `pair_granular.rst:149`; `GranSubModNormalJKR` `gran_sub_mod_normal.cpp:389` | ⚠️ `adhesion_model="jkr"`, **simplified explicit** `F_adh≈³⁄₂πγR*`, **Hertz path only** — `crates/dirt_granular/src/contact.rs:397` (pull-off range 263–272) |
| **SJKR** (area-proportional cohesion) | ❌ — not a `pair granular` model | ✅ both Hertz + Hooke paths — `crates/dirt_granular/src/contact.rs:408` (Hertz), `845` (Hooke) |
| **MDR** elastic-plastic adhesive (Zunker & Kamrin) | ✅ `mdr` — `pair_granular.rst:186`; `GranSubModNormalMDR` `gran_sub_mod_normal.cpp:540` (+ companion `fix_granular_mdr.cpp`) | ❌ — searched `crates/` for `mdr`; not found |
| `none` (null normal) | ✅ `gran_sub_mod_normal.cpp:133` | — (n/a) |

**Notes.** DIRT's JKR is a *simplified* explicit pull-off force with an extended
interaction range, not LAMMPS's full contact-radius root-solve; treat as
approximate. DMT/JKR `surface_energy` is **silently ignored on the Hooke path**
in DIRT (documented, `crates/dirt_granular/src/lib.rs:35`). DIRT has **SJKR**
that LAMMPS `pair granular` lacks; LAMMPS has **MDR** that DIRT lacks.

## Normal damping

LAMMPS exposes damping as an independent sub-model (`F_{n,damp}=−η_n v_{n,rel}`,
`pair_granular.rst:315`); default `viscoelastic` (`:368`). DIRT bakes normal
damping into each contact model (viscoelastic for Hertz, linear for Hooke),
calibrated from restitution.

| Feature | LAMMPS | DIRT |
|---|---|---|
| `velocity` (`η=η₀`) | ✅ `pair_granular.rst:346`; `gran_sub_mod_damping.h:56` | ✅ Hooke linear damping `γ_n v_n` — `crates/dirt_granular/src/contact.rs:845` |
| `mass_velocity` (`η=η₀ m_eff`) | ✅ `pair_granular.rst:356`; `.h:64` | ✅ effectively (mass-scaled) via Hooke path — same cite |
| `viscoelastic` (`η=η₀ a m_eff`) | ✅ `pair_granular.rst:368`; `.h:72` | ✅ Hertz viscoelastic `β√(S_n m_r)` — `crates/dirt_granular/src/contact.rs:345` |
| `tsuji` (restitution-calibrated) | ✅ `pair_granular.rst:381`; `.h:80` | ✅ DIRT's Hertz damping is restitution-calibrated (β from COR) — `crates/dirt_granular/src/contact.rs:345` |
| `coeff_restitution` | ✅ `pair_granular.rst:406`; `.h:89` | ⚠️ same intent (COR→β) but not a separate selectable keyword |
| `mdr` damping | ✅ `pair_granular.rst:431`; `.h:97` | ❌ (no MDR) |
| **Cundall non-viscous** global damping (`γ_l, γ_a`) | ✅ `fix_damping_cundall.cpp` — `doc/src/fix_damping_cundall.rst` | ❌ — searched `crates/dirt_fixes`; not found (DIRT has *viscous* velocity damping instead, below) |
| **Viscous velocity damping** fix (`F=−γv`) | ✅ via `fix viscous` (MISC) | ✅ `crates/dirt_fixes/src/lib.rs:186` (`ViscousDef`, `apply_viscous` 533) |

## Tangential models

| Feature | LAMMPS | DIRT |
|---|---|---|
| `linear_nohistory` (velocity Coulomb, no history) | ✅ `pair_granular.rst:488`; `.h:59` | ✅ selectable `tangential_model = "linear_nohistory"` — history-free `F_t = -min(μ\|F_n\|, γ_t\|v_t\|) t̂`, zero accumulated ξ — `crates/dirt_granular/src/contact.rs` (tangential block); validated by `examples/bench_nohistory_tangential` |
| `linear_history` (spring on accumulated ξ) | ✅ `pair_granular.rst:550`; `.h:71` | ✅ Hooke-path linear `k_t` + Coulomb cap — `crates/dirt_granular/src/contact.rs:722` (tangential block) |
| **Mindlin** (adds contact-radius `a`, `k_t=8G*`) | ✅ `pair_granular.rst:613`; `.h:98` | ✅ incremental spring-history + Coulomb cap `μ|F_n|`, damping `γ_t=2β√(5/6)√(k_t m_r)` — `crates/dirt_granular/src/contact.rs:432` (`k_t` ~348, damping 484) |
| `mindlin/force` (history stored as elastic force) | ✅ `pair_granular.rst:645`; `.h:112` | ❌ — DIRT stores tangential *displacement* history (`crates/dirt_granular/src/tangential.rs:49` `ContactHistoryStore`), not the force form |
| `mindlin_rescale` / `mindlin_rescale/force` (rescale history on unloading) | ✅ `pair_granular.rst:675, 700`; `.h:119, 126` | ✅ selectable `tangential_model = "mindlin_rescale"` or `"mindlin_rescale/force"` (`"mindlin_rescale_force"` alias); scales tangential history by `a/a_prev` on unloading — `crates/dirt_granular/src/contact.rs`; validated by `examples/bench_mindlin_rescale_tangential` |
| `linear_history_classic`, `mindlin_classic` (legacy, source-only) | ✅ `gran_sub_mod_tangential.h:83, 91` (undocumented) | ❌ (n/a) |
| Per-contact tangential history store | ✅ implicit in sub-models | ✅ canonical-frame store, 8 f64/contact (tangential/rolling/twist plus previous contact radius for rescale) — `crates/dirt_granular/src/tangential.rs:49` |

## Rolling resistance

| Feature | LAMMPS | DIRT |
|---|---|---|
| `none` (default) | ✅ `pair_granular.rst:722`; `.h:37` | ✅ (default is **constant-torque**, see below) |
| **SDS** rolling (spring-dashpot-slider, capped by `μ_roll F_n`) | ✅ `sds` — `pair_granular.rst:716`; `GranSubModRollingSDS` `.h:45` | ✅ `rolling_model="sds"` — `crates/dirt_granular/src/contact.rs:541` (Hooke path 958) |
| **Constant-torque** rolling (`τ=μ_r|F_n|R*`) | ❌ — not a `pair granular` rolling model | ✅ default — `crates/dirt_granular/src/contact.rs:585` |
| Rolling history slot | ✅ implicit | ✅ `crates/dirt_granular/src/tangential.rs:44` (`[3..6]`) |

## Twisting friction

| Feature | LAMMPS | DIRT |
|---|---|---|
| `none` (default) | ✅ `pair_granular.rst:784`; `.h:38` | ✅ (default is **constant-torque**) |
| **SDS** twisting (`τ=−k ξ−γΩ`, capped `μ_tw F_n`) | ✅ `sds` — `pair_granular.rst:786`; `.h:58` | ✅ `twisting_model="sds"` — `crates/dirt_granular/src/contact.rs:606` |
| **Marshall** twisting (coeffs derived from tangential model) | ✅ `marshall` — `pair_granular.rst:818`; `GranSubModTwistingMarshall` `.h:46` | ✅ `twisting_model="marshall"` — `crates/dirt_granular/src/contact.rs:614` (`k_twist=½k_t a²`, `γ_twist=½γ_t a²`, `μ_twist=⅔ a μ_t`; Hooke path 1071). Validated by `examples/bench_marshall_twisting` |
| **Constant-torque** twisting (`τ=μ_tw|F_n|R*`) | ❌ | ✅ default — `crates/dirt_granular/src/contact.rs:638` |
| Twisting history slot | ✅ implicit | ✅ `crates/dirt_granular/src/tangential.rs:46` (`[6]`) |

## Heat / thermal

| Feature | LAMMPS | DIRT |
|---|---|---|
| Contact **heat conduction**, radius model `Q=2 k_s a ΔT` | ✅ `heat radius` — `pair_granular.rst:879`; `GranSubModHeatRadius` `gran_sub_mod_heat.h:46` | ❌ — searched `crates/` for `heat_conduction\|thermal_conduct\|conduction`; **not found** |
| Contact heat conduction, area model `Q=h_s π a² ΔT` | ✅ `heat area` — `pair_granular.rst:892`; `.h:58` | ❌ — not found |
| Per-atom **heat-flow integration** (temperature update) | ✅ `fix_heat_flow.cpp` — `doc/src/fix_heat_flow.rst` | ❌ — not found |
| External **heat source** | ✅ `fix_add_heat.cpp` | ❌ — not found |
| **Wall temperature** boundary | ✅ `fix wall/gran ... temperature` — `fix_wall_gran.rst:196` | ⚠️ **stub only** — `temperature: Option<f64>` field with no transfer physics, `crates/dirt_wall/src/lib.rs:303` (future-hook comment `:44`) |
| Granular *temperature* (velocity-fluctuation **diagnostic**, not a thermostat) | ✅ via `compute` machinery | ✅ opt-in output — `crates/dirt_granular/src/granular_temp.rs:33` (**not** in `GranularDefaultPlugins`, `lib.rs:112`) |

**Thermal is DIRT's largest gap:** no contact conduction, no heat-flow
integrator, no wall-temperature transfer. The only thermal artifact is a
velocity-fluctuation *diagnostic* (unrelated to real heat transfer).

## Walls (`fix wall/gran` vs `dirt_wall`)

| Feature | LAMMPS | DIRT |
|---|---|---|
| **Plane** wall (axis-aligned) | ✅ `xplane/yplane/zplane` — `fix_wall_gran.rst:156` | ✅ `WallPlane` — `crates/dirt_wall/src/lib.rs:351` |
| **Cylinder** wall | ✅ `zcylinder` (Z-axis only) — `fix_wall_gran.rst:158` | ✅ `WallCylinder` X/Y/Z + finite axial bounds, inside/outside — `crates/dirt_wall/src/lib.rs:421` |
| **Sphere** wall | ❌ — not a `fix wall/gran` primitive | ✅ `WallSphere` — `crates/dirt_wall/src/lib.rs:460` |
| **Region** wall (arbitrary geometry) | ✅ `fix wall/gran/region` — `fix_wall_gran_region.cpp`; `fix_wall_gran.rst:77` | ✅ `WallRegion` (block/sphere/cylinder/cone/union) — `crates/dirt_wall/src/lib.rs:484` |
| Wall **normal**: Hooke *and* Hertz | ✅ `hooke`/`hooke/history`/`hertz/history`/`granular` — `fix_wall_gran.rst:91` | ✅ Hertz and Hooke — `crates/dirt_wall/src/lib.rs:1175` (`wall_normal_force`); Hooke wall rebound validated by `examples/bench_hooke_wall_rebound` |
| Wall full sub-model set (rolling/twisting/heat via `granular`) | ✅ `fstyle granular` — `fix_wall_gran.rst:113` | ⚠️ partial (see below) |
| Wall **tangential** (Mindlin sliding), all shapes | ✅ | ✅ `wall_tangential_force` — `crates/dirt_wall/src/lib.rs:926` |
| Wall **rolling** (constant + SDS), all shapes | ✅ (via `granular`) | ✅ `wall_rolling_torque` — `crates/dirt_wall/src/lib.rs:1004` |
| Wall **twisting** | ✅ (via `granular`) | ⚠️ **plane walls only** — `crates/dirt_wall/src/lib.rs:1204` (note `:23`) |
| Wall **JKR/DMT** adhesion | ✅ (via `granular`) | ⚠️ **plane walls only** (SJKR on all) — `crates/dirt_wall/src/lib.rs:1173` (note `:33`) |
| Wall **heat/temperature** transfer | ✅ `temperature` keyword — `fix_wall_gran.rst:196` | ❌ stub field only (see Thermal) |
| Wall **motion** | ✅ `wiggle` (oscillate), `shear` (constant tangential) — `fix_wall_gran.rst:165, 186` | ✅ Static / ConstantVelocity / Oscillate / **Servo** (plane only) — `crates/dirt_wall/src/lib.rs:311` (`ServoDef` 206) |
| Per-atom wall **contact output** | ✅ `contacts` keyword — `fix_wall_gran.rst:224` | ⚠️ via general contact-analysis (geometry) — `crates/dirt_contact_analysis` |

**Notes.** DIRT walls add **sphere** primitives and a **servo** (force-controlled)
motion mode LAMMPS lacks; LAMMPS walls support a **Hooke** normal and
**temperature** boundary DIRT lacks. DIRT's curved/region walls are missing
twisting and JKR/DMT (plane-only).

## Bonded particles (BPM ↔ dirt_bond)

LAMMPS handles bonded/cohesive solids in `src/BPM/` (not `pair granular`); DIRT
has a first-class beam-theory bond model.

| Feature | LAMMPS (BPM) | DIRT (`dirt_bond`) |
|---|---|---|
| Bond model | ✅ `bond_bpm_spring` (normal spring), `bond_bpm_spring_plastic`, `bond_bpm_rotational` — `doc/src/bond_bpm_rotational.rst` | ✅ beam-theory (cylindrical) 4-channel: normal/shear/twist/bending — `crates/dirt_bond/src/lib.rs:521` (`DemBondPlugin`, `bond_force` 938) |
| Stiffness from material (`E, G`) vs direct | ⚠️ rotational bond uses explicit spring constants | ✅ both Material (`E,G`) and Direct modes — `crates/dirt_bond/src/lib.rs:214` |
| Per-channel viscous damping | ✅ (bond damping coeff) | ✅ critical-damping ratio per channel — `crates/dirt_bond/src/lib.rs:238` |
| Breakage criteria | ✅ stress/strain break, `fix_update_special_bonds` | ✅ rich set: axial/combined/interaction force·stress·strain — `crates/dirt_bond/src/breakage.rs:256` (enum 562) |
| Stochastic break thresholds (Weibull) | ❌ — not built-in | ✅ length-scaled Weibull inverse-CDF — `crates/dirt_bond/src/breakage.rs:179` |
| Bond **plasticity** | ✅ `bond_bpm_spring_plastic` (axial) | ✅ bending (Guo EPP/trilinear/piecewise) + axial hardening — `crates/dirt_bond/src/plasticity.rs:79` |
| Auto-bond touching / load from file | ✅ bonds created by `create_bonds`/read_data | ✅ `auto_bond_touching` / `load_bonds_from_file` — `crates/dirt_bond/src/lib.rs:662, 720` |
| Bond metrics/output | ✅ `compute_nbond_atom` | ✅ `BondMetrics` — `crates/dirt_bond/src/lib.rs:501` |

DIRT's bond model is notably richer on **breakage distributions (Weibull)** and
**bending plasticity**; LAMMPS BPM offers a `rotational` moment-carrying bond
family with a mature special-bonds exclusion pipeline.

## Rigid bodies / clumps (RIGID ↔ dirt_clump)

| Feature | LAMMPS (RIGID) | DIRT (`dirt_clump`) |
|---|---|---|
| Multisphere rigid clumps | ✅ `fix rigid` (+ `fix_rigid_nh`) — `doc/src/fix_rigid.rst` | ✅ `ClumpPlugin` + `MultisphereBody` — `crates/dirt_clump/src/lib.rs:308` |
| Clump definition / insertion | ✅ molecule template + `fix pour` | ✅ TOML `ClumpDef`/`ClumpInsertConfig` — `crates/dirt_clump/src/lib.rs:155, 164` |
| Rigid-body integration | ✅ (NVE/NH rigid) | ✅ angular-momentum + Richardson-iteration Euler — `crates/dirt_clump/src/body.rs:153`, `lib.rs:378` |
| Intra-body contact exclusion | ✅ (special bonds / rigid exclusion) | ✅ `same_body` guard — `crates/dirt_clump/src/lib.rs:20` |

Rough parity on rigid clumps.

## Fixes / integrators / diagnostics

| Feature | LAMMPS | DIRT |
|---|---|---|
| Velocity-Verlet translational integrator | ✅ `fix nve` | ✅ `VelocityVerletPlugin` (soil) — `crates/dirt_granular/src/lib.rs:186` |
| Rotational (quaternion) integrator | ✅ `fix nve/sphere` | ✅ `RotationalDynamicsPlugin` — `crates/dirt_granular/src/rotational.rs:49` |
| Gravity / body force | ✅ `fix gravity` | ✅ `GravityPlugin` — `crates/dirt_fixes/src/lib.rs:633` |
| Add/Set force | ✅ `fix addforce`/`setforce` | ✅ `AddForceDef`/`SetForceDef` — `crates/dirt_fixes/src/lib.rs:66, 96` |
| Prescribed motion | ✅ `fix move` | ✅ `MoveLinearDef` — `crates/dirt_fixes/src/lib.rs:129` |
| Freeze particles | ✅ `fix freeze` — `fix_freeze.cpp` | ✅ `FreezeDef` — `crates/dirt_fixes/src/lib.rs:165` |
| Displacement-limited NVE | ✅ `fix nve/limit` | ✅ `NveLimitDef` — `crates/dirt_fixes/src/lib.rs:219` |
| Particle insertion / pour | ✅ `fix pour` — `fix_pour.cpp` | ✅ `DemAtomInsertPlugin` (+ clump insert) |
| **Synchronized Verlet** (polydisperse frictional) | ✅ `synchronized_verlet` — `pair_granular.rst:854` | ❌ — not found |
| Fabric / anisotropy tensor | ✅ `compute_fabric.cpp` — `doc/src/compute_fabric.rst` | ✅ fabric tensor — `crates/dirt_contact_analysis/src/lib.rs` |
| Per-atom contact count | ✅ `compute_contact_atom.cpp` | ✅ coordination number / rattlers — `crates/dirt_contact_analysis/src/lib.rs` |
| Virial / stress | ✅ built-in | ✅ `VirialStressPlugin` — `crates/dirt_granular/src/contact.rs:67` |
| Flow-rate measurement plane | ⚠️ via region counts | ✅ `MeasurePlanes` — `crates/dirt_measure_plane/src/lib.rs:199` |

---

## Gap summary (seeds for follow-up goals)

**DIRT lacks (→ candidate `implement` goals):**

1. **Contact heat conduction** — LAMMPS `heat radius`/`heat area`
   (`pair_granular.rst:879, 892`) + per-atom heat-flow integration
   (`fix_heat_flow.cpp`). DIRT has only a wall-temperature *stub* and a
   velocity-fluctuation diagnostic. **Largest gap.**
2. **MDR elastic-plastic normal model** (`gran_sub_mod_normal.cpp:540` +
   `fix_granular_mdr.cpp`) — DIRT has no plastic normal contact.
3. **Mindlin `force`-form history without unloading rescale** (`mindlin/force`,
   `pair_granular.rst:645`) — DIRT now has the unloading-rescale displacement and
   force-history variants, but not the standalone `mindlin/force` keyword.
4. ~~**Marshall twisting** (coeffs derived from the tangential model,
   `pair_granular.rst:818`)~~ — ✅ **CLOSED**: `twisting_model="marshall"`
   (`crates/dirt_granular/src/contact.rs:614`) derives `k_twist=½k_t a²`,
   `γ_twist=½γ_t a²`, `μ_twist=⅔ a μ_t` from the active tangential model;
   validated by `examples/bench_marshall_twisting` (spin-down vs the analytical
   Marshall torque, 0.02% error). DIRT now has all three twist forms
   (constant / SDS / Marshall).
5. ~~**`linear_nohistory` tangential** (`pair_granular.rst:488`) — DIRT is always
   history-based.~~ **Closed:** `tangential_model = "linear_nohistory"` adds the
   history-free velocity-Coulomb law (`examples/bench_nohistory_tangential`).
6. **Cundall non-viscous global damping** (`fix_damping_cundall.cpp`).
7. **Curved/region-wall twisting + JKR/DMT**
   (currently plane-only, `crates/dirt_wall/src/lib.rs:1204, 1173`).
8. **Wall temperature boundary** transfer (stub today,
   `crates/dirt_wall/src/lib.rs:303`).
9. **Full JKR contact-radius solve** and **JKR/DMT on the Hooke path**
   (DIRT's JKR is simplified & Hertz-only, `crates/dirt_granular/src/contact.rs:397`;
   Hooke-path adhesion ignored, `crates/dirt_granular/src/lib.rs:35`).

**DIRT has that LAMMPS `pair granular` lacks (validate / keep):**

- **SJKR** area-proportional cohesion (`crates/dirt_granular/src/contact.rs:408`).
- **Constant-torque** rolling & twisting defaults
  (`crates/dirt_granular/src/contact.rs:585, 638`).
- **Sphere** wall primitive and **servo** (force-controlled) wall motion
  (`crates/dirt_wall/src/lib.rs:460, 206`).
- **Weibull** stochastic bond-break thresholds and **bending plasticity**
  (`crates/dirt_bond/src/breakage.rs:179`, `plasticity.rs:79`).

**Recommended validation targets (→ candidate `validate` goals):** run the
equivalent LAMMPS case (binary is on PATH) for DIRT's Hertz normal, Mindlin
tangential, SDS rolling, and DMT cohesion, and compare against the documented
formulas in `pair_granular.rst` for external provenance.

---

*Sources: LAMMPS `~/projects/reference/lammps` (`doc/src/pair_granular.rst`,
`doc/src/fix_wall_gran.rst`, `src/GRANULAR/`, `src/BPM/`, `src/RIGID/`); DIRT
crates `dirt_granular`, `dirt_wall`, `dirt_bond`, `dirt_clump`, `dirt_fixes`,
`dirt_contact_analysis`, `dirt_measure_plane`. Verified by grep; absences
annotated with the search performed.*
