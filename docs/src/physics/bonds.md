# Parallel Bonds

The **Bonded Particle Model (BPM)** glues particle pairs together with an
elastic beam that resists four independent deformations — for fibers,
agglomerates, and cohesive solids that can stretch, bend, twist, yield, and
break. Bonds live in the `dirt_bond` crate; add `DemBondPlugin` to enable them.

## Bonds vs. contact

A bonded pair is governed by the **bond** law, not the granular contact law:
while a bond is intact, granular contact between that pair is **suppressed**
(via the substrate's `BondStore`). When a bond breaks, the pair reverts to
ordinary Hertz–Mindlin contact. `DemBondPlugin` therefore needs the granular
contact plugin present so that non-bonded pairs still interact:

```rust
app.add_plugins(CorePlugins)
   .add_plugins(GranularDefaultPlugins) // Hertz–Mindlin contact + Verlet
   .add_plugins(DemBondPlugin);         // bonds; contact suppressed on bonded pairs
```

## The four deformation channels

Each bond is a solid cylindrical beam of radius `r_b`, resisting:

| Channel | Stiffness (beam form) | Physical meaning |
|---|---|---|
| Normal (extension/compression) | `K_n = E · A / L` | stretching along the bond axis **n̂** |
| Shear | `K_t = G · A / L` | sliding perpendicular to **n̂** |
| Twist (torsion) | `K_tor = G · J / L` | rotating about **n̂** |
| Bending | `K_bend = E · I / L` | relative rotation perpendicular to **n̂** |

with `A = π r_b²`, `J = ½ π r_b⁴`, `I = ¼ π r_b⁴ = ½ J`, and `L = r₀` the
equilibrium bond length. The bond radius per pair is `r_b = bond_radius_ratio ·
min(R_i, R_j)`; `bond_radius_ratio = 1.0` makes bonds as wide as the smaller
particle.

The normal force is `F_n = (K_n δ + γ_n v_n) n̂` with `δ = |r_ij| − r₀`. The
shear force is history-dependent — the shear displacement `Δs` is re-projected
perpendicular to the current `n̂` each step — and is evaluated at the bond
mid-point, producing a lever-arm torque `τ_shear = (L/2) n̂ × F_t` on both
particles.

### Sign convention: lower tag owns

Shear force and both moments are applied as **`+F` on atom `i` (the lower tag)
and `−F` on atom `j`**. This is the LIGGGHTS/Fortran-BPM convention: it damps
the relative motion symmetrically and keeps the bond's accounting consistent no
matter which particle the neighbor list visits first. The per-bond shear
history is integrated across steps and re-projected each step (its lifecycle is
created at bonding and destroyed at breakage).

## Damping

Per-channel viscous damping comes from a **critical-damping ratio** `β ∈ [0, 1]`
(critical = 1.0):

```text
γ   = 2 β √( m* · K_eff )      for F_n, F_t
γ_M = 2 β √( I* · K_eff )      for M_tor, M_bend
```

using the reduced mass `m* = m_i m_j / (m_i + m_j)` and reduced moment of
inertia `I*` of the bonded pair. Each channel accepts an optional raw-`γ`
override that bypasses the β-based calculation.

## Breakage and plasticity

A bond breaks when either combined stress at the extreme fibre exceeds its
limit:

```text
σ = F_n / A  +  2 |M_bend| r_b / J     →  break if σ > σ_max   (tensile)
τ = |F_t| / A  +  |M_tor| r_b / J      →  break if τ > τ_max   (shear)
```

The combined-stress criterion above is one of **nine** breakage `kind`s, set in
the `[bonds.breakage]` sub-table (absent → bonds never break):

| `kind` | Family | Trips on |
|---|---|---|
| `unbreakable` | — | never (explicit) |
| `axial_force` | two-branch | axial force / shear force |
| `axial_stress` | two-branch | axial stress / shear stress |
| `axial_strain` | two-branch | axial strain / shear strain |
| `combined_stress` | two-branch | beam-stress at the extreme fibre (tensile / shear) |
| `combined_strain` | two-branch | beam-strain at the extreme fibre |
| `interaction_linear_force` | envelope | weighted sum of axial/shear/bending/twist forces ≥ 1 |
| `interaction_linear_stress` | envelope | same, in stress |
| `interaction_linear_strain` | envelope | same, in strain |

The two-branch families take a `tensile` threshold (required) and a `shear`
threshold (optional; absent → ∞). The `interaction_linear_*` families take up to
four optional channels (`axial`, `shear`, `bending`, `twist`); absent channels
drop out of the sum. Each threshold is a **distribution**:

- `{ kind = "constant", value = … }` — the same strength for every bond.
- `{ kind = "weibull", mean, m, l_calib, l_min }` — a length-scaled
  two-parameter Weibull, so a population breaks with realistic scatter (`m` is
  the Weibull modulus, `l_calib` the calibration length).
- `{ kind = "crack_band", value_ref, l_ref, eps_yield, l_min }` — a deterministic
  Bažant crack-band rescaling of the threshold by bond length.

**Crack-band regularization** makes fracture energy mesh-invariant: pair a
`crack_band` threshold with a *strain* criterion (`axial_strain`,
`combined_strain`, `interaction_linear_strain`) and the per-bond strength scales
with bond length so that the dissipated energy per crack stays constant as you
refine the particle size.

### Plasticity

Bending and axial channels can yield independently, set in
`[bonds.plasticity.bending]` and `[bonds.plasticity.axial]`:

- `guo_bending` — Guo (2018) elastic-perfectly-plastic bending, yield moment
  `M^p = (4/3) σ_0 r_b³` from the material yield stress `yield_stress`.
- `guo_trilinear` — a trilinear extension; requires `[bonds].youngs_modulus` to
  be set (it derives the elastic strain `ε_e = σ_0/E`).
- `piecewise` — an arbitrary piecewise-linear curve from `breakpoint_strains` and
  `slope_multipliers`, with an optional `length_calibration` for crack-band
  regularization. Available on **both** bending and axial channels.

Bending breakpoint strains are always **extreme-fibre** strains
(`ε = r_b · θ_bend / l_b`), not angles, which makes the material parameters
geometry-independent. Axial plasticity handles tension and compression
symmetrically and tracks `|ε_axial|` for kinematic hardening.

## Configuration

You parametrise stiffness one of two ways:

- **Material mode** (paper-standard): give `youngs_modulus` *E* and
  `shear_modulus` *G*; per-bond stiffnesses derive from beam theory.
- **Direct mode**: give `normal_stiffness`, `shear_stiffness`,
  `twist_stiffness`, `bending_stiffness` directly. Used when E/G are not set.

If both are set, material mode wins for whichever channels E/G apply to.

Bonds are created either by **auto-bonding** at setup (`auto_bond = true` bonds
every pair within `bond_tolerance × (R_i + R_j)`) or loaded from an explicit
bond file (`auto_bond = false` with a `file` / `format`).

```toml
[bonds]
auto_bond = true
bond_tolerance = 1.001
bond_radius_ratio = 1.0
ghost_cutoff_multiplier = 2.5   # extend MPI ghost reach for bonded pairs + 1-3 exclusion

# Material mode
youngs_modulus = 1.0e9      # E (Pa)
shear_modulus  = 4.0e8      # G (Pa)

# Damping ratios (critical = 1.0)
beta_normal  = 0.05
beta_shear   = 0.05
beta_twist   = 0.05
beta_bending = 0.05

# Breakage (combined-stress with Weibull tensile threshold)
seed = 0
[bonds.breakage]
kind    = "combined_stress"
tensile = { kind = "weibull", mean = 5.0e7, m = 8.0, l_calib = 0.020 }
shear   = { kind = "constant", value = 3.0e7 }

# Plasticity — bending and axial channels independently configurable
[bonds.plasticity.bending]
kind         = "guo_bending"
yield_stress = 1.23e8

[bonds.plasticity.axial]
kind                = "piecewise"
breakpoint_strains  = [0.01, 0.02, 0.03]
slope_multipliers   = [0.5, 0.1, 0.0]
```

> The config also carries a `bond_type` field that is *stored but currently
> unused* — a placeholder for future per-bond type dispatch.

Auto-bonding runs once at `PostSetup`, before the neighbor setup bakes the ghost
cutoff. At that point — on a multi-rank run — all atoms are still on rank 0
(pre-migration), so every bond is created on rank 0 and then distributed during
the first exchange. Bonds wrap correctly across periodic boundaries (minimum
image on the periodic axes).

## Bond/contact exclusion

While a bond is intact, the granular contact law **skips** that pair: the
substrate's `BondStore` is registered as an `AtomData` and the Hertz–Mindlin
loop consults it. `DemBondPlugin` re-adds `soil_core::BondPlugin` to install this
hook, which is why it must come *after* `GranularDefaultPlugins`. On breakage the
pair is removed from `BondStore` on both ends and ordinary contact resumes on the
next step. (A bond that breaks this step still contributes its force for that one
final step before it disappears — removal is deferred to after the force pass.)

## MPI: ghost cutoff and reproducibility

Two distributed-run concerns are worth understanding:

- **Ghost cutoff.** A bond couples two atoms that may sit on different ranks, so
  the owning rank must see its partner as a ghost. `ghost_cutoff_multiplier`
  (default `2.5`) bumps `domain.ghost_cutoff` to `max_r0 × multiplier`. The 2.5
  budget covers 1× for the bond itself, 2× for shared-neighbour (1-3) exclusion,
  and a 0.5× margin for stretch. If the `bond_missing` thermo key is nonzero,
  the cutoff is too small — raise the multiplier. Set it to `0.0` only on a
  single-process run.
- **Reproducible breakage.** Per-bond failure thresholds are drawn from the
  canonicalized tag pair `(min(tag_i, tag_j), max(…), seed)`, so two ranks
  visiting the same bond from opposite sides compute *identical* thresholds. The
  same `seed` plus the same topology gives the same breakage pattern regardless
  of how the domain is decomposed — bit-for-bit across rank counts.

## Thermo output

`DemBondPlugin` publishes three keys each thermo step:

- `bond_strain` — mean axial strain `δ/r₀` over all live bonds (0 if none).
- `bonds_broken` — cumulative bonds broken since the run started.
- `bond_missing` — bonds skipped this step because the partner atom was not
  visible (ghost cutoff too small, or an MPI cut runs through a bond).

The last key is the diagnostic for distributed runs: if `bond_missing` is
nonzero, raise `ghost_cutoff_multiplier`.

## Worked examples

- `bond_cantilever` — a 10-sphere bonded chain anchored by `[[freeze]]`, bending
  under gravity and settling via critical bond damping.
- `bond_fiber_tensile` / `fiber_bond` — fiber tensile and crossover tests.
- `bench_fiber_crossover` — explicit bonds loaded from a file.

See the [Validation chapter](../reference/validation.md) and the related
[fix freeze-vs-pin](fixes.md#freeze-vs-pin-the-dirtsoil-split) discussion.
