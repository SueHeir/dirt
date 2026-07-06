# SPH glass-sphere calibration

The complete set of DEM simulations that calibrate the **canonical glass-sphere
material** into the closures the **dev_soil_sph continuum model** consumes. Glass beads
are the *validation* material — chosen because we have experimental data at every
rung (cross-code shear stress vs Φ from Fortran/LIGGGHTS/LAMMPS, column-collapse
runout scaling from Lube/Lajeunesse, angle of repose, hopper). Each deliverable
below is owned by its own agent and lives in its own subfolder.

## Canonical material (shared by every sim)

```toml
[[dem.materials]]
name           = "glass"
youngs_mod     = 7.0e7    # softened from ~65 GPa real glass (rigid-grain limit)
poisson_ratio  = 0.245    # glass
restitution    = 0.926    # measured glass–glass COR (exact-COR damping → realized = nominal)
friction       = 0.16     # measured glass–glass sliding μ_p
# rolling_friction = <pinned by the repose deliverable>

# NOTE: density is a per-insert property, not a material property — it goes on
# the particle insert block, NOT in [[dem.materials]]:
[[particles.insert]]
density        = 2500.0   # kg/m³ (glass)
```

## Deliverables (each = one DEM sim, one agent, one SPH closure)

| # | DEM experiment | SPH closure it informs | basis bench |
|---|---|---|---|
| 1 | **Shear rheology** (Lees–Edwards rheometer) | `μ(I)`, `Φ(I)`, `ρ_c = Φ_max ρ_s` | `bench_lebc_shear` |
| 2 | **Isotropic compression** | bulk modulus `K` / EOS slope | (new) |
| 3 | **Angle of repose** | rolling friction `μ_r` (pins macroscopic friction) | `bench_angle_of_repose` |
| 4 | **Enduring-contact stress** (dense rheometer residual) | `σ_contact(Φ)` = `p_DEM − p_KT` (the rate-independent branch) | `bench_lebc_shear` |
| 5 | **Granular-temperature dissipation** (Haff cooling) | dissipation coefficient / cooling rate | `bench_sphere_haff_cooling` |
| 6 | **Granular-temperature conductivity** (vibro-fluidized bed) | `κ(Φ)` (fluctuation-energy conduction) | `bench_granular_conductivity` |
| 7 | **Column collapse** (macro validation) | end-to-end runout scaling vs experiment | `bench_column_collapse` |
| 8 | **Cooperativity length** (velocity-fluctuation correlation under shear) | nonlocal amplitude `A` in `ξ(μ)=A·d/√(μ−μ_s)` + the `g∝√T` bridge (dev_soil_sph's nonlocal/NGF branch) | `bench_lebc_shear` |

(1)–(3) feed the **μ(I)/critical-state** rheology; (4)–(6) feed the
**granular-temperature de-fluidization** two-branch model; (7) is the macro
validation gate the calibrated SPH must also reproduce; (8) feeds the **nonlocal
cooperativity** branch (creep below yield) — `coop_amplitude` in MaterialParams,
but remains exploratory until it has an independent bounded gate.

## Automation harness: bounded smoke gates vs. full calibration runs

Several of these deliverables are **long-form calibration sweeps** — multi-case
(Φ, γ̇) grids or deep-bed runs whose scientific point is a fitted closure, not a
single fast check. Run in full they take far longer than the 1800 s per-bench cap
the hourly automation harness enforces (it runs each `sweep.py`), so they
**timed out every run and validated nothing**. Following the same pattern already
adopted for `bench_lebc_shear` / `bench_granular_conductivity` / `bench_hopper_beverloo`,
`sweep.py` with **no argument** now runs a **bounded smoke gate** with a genuine
PASS/FAIL, and the full calibration is preserved behind an explicit `full` command:

| # | deliverable | `sweep.py` (no-arg) gate | full run | asserted criterion |
|---|---|---|---|---|
| 1 | shear rheology | bounded gate | `sweep.py full` | steady, physical μ = \|σ_xy\|/P ∈ [0.15, 0.90] per case |
| 4 | enduring contact | bounded gate | `sweep.py full` | σ_contact = p_DEM − p_KT(Φ,T) self-consistent + branch opens with Φ |
| 6 | conductivity | bounded gate | `sweep.py full` | Fourier-law κ*(Φ) > 0, finite, order-unity vs KT |
| 8 | cooperativity length | **excluded** | `sweep.py --run` | no fixed target for `A` or `g∝√T` is justified yet |

The full runs (config.toml sweeps, fits, and plots) are **unchanged** and remain
the scientific calibration; the gates are additive fast breakage checks, never a
substitute for them. Deliverables 2 and 5 already complete inside the cap and gate
directly; 3 and 7 are tracked separately (angle-of-repose and column-collapse
honesty goals). Deliverable 8 is kept as a documented measurement rig, but its
default `sweep.py` path reports `SKIPPED` when no data exists rather than turning
missing exploratory data into a zero-duration FAIL. Each gated subfolder README
documents its own gate under **“Bounded smoke gate (CI)”**.

## Layout (per subfolder, following the dirt-example conventions)

```
SPH_glass_sphere_calibration/
  README.md
  01_shear_rheology/      main.rs · config.toml · sweep.py · plots/   (data/ sweep/ gitignored)
  02_compressibility/
  03_angle_of_repose/
  04_enduring_contact/
  05_cooling_dissipation/
  06_conductivity/
  07_column_collapse/
  08_cooperativity_length/
  calibration.yaml        # the assembled closure consumed by sph_constitutive::MaterialParams
```

Each subfolder uses the canonical material above, runs at the glass-sphere
parameters, and emits its closure into the shared `calibration.yaml`.
