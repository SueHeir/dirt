# Angle-of-Repose — Rolling-Friction (μ_r) Calibration Gate

**DEM formation study, DIRT step 03 (not an SPH solver calibration).** Forms a static granular heap and
measures its **angle of repose** θ_r while sweeping the **rolling friction μ_r**
at a **fixed sliding friction μ_p = 0.16** (the measured glass value). The goal is
to characterize whether this protocol reaches the historical glass repose interval
**[22, 26]°**. It does not establish a transferable canonical-material μ_r, and
the `graph` command intentionally reports `CALIBRATION: WITHHELD` even when its
local formation-study checks pass.  The local criteria are reported separately;
they are not silently converted into material validation.

A repository-authored JSON transcription, a Crossref bibliographic lookup, and
an in-band solver result are not independent experimental validation. Closure
requires an independently reviewed primary measurement matching the material,
surface condition, vessel, deposition, release, and measurement definition,
plus a separate cross-substrate validation of a defined SPH constitutive
quantity. No such evidence is committed here.

## Formation restitution and the non-transfer boundary

The calibration configuration forms its deposit at **e = 0.4**, while canonical
glass is listed at e = 0.926. That is a different formation protocol. Although
friction is important to static stability, restitution can change the packing and
contact network reached during deposition; therefore this study **does not claim
that a μ_r fitted at e = 0.4 transfers to e = 0.926**. A transfer claim requires
a separately completed, protocol-matched e = 0.926 validation.

The lift-the-cylinder protocol at e = 0.4 is the proven `bench_angle_of_repose`
protocol; that bench confirms e = 0.4 + the lift sequence gives clean, settled
deposits (it does not bounce apart).

The independent LAMMPS SDS formulation and the static sliding-feasibility
boundary are recorded in [LAMMPS_SDS_ORACLE.md](LAMMPS_SDS_ORACLE.md).  In
particular, `atan(0.16)=9.09°` is far below the retained 22--26° glass band;
rolling resistance is a capped couple, not an extra sliding-force allowance.
This is why the gate remains fail-closed pending a qualified, cross-code checked
physical result rather than changing a target or estimator.

> **STATUS (this calibration): FAIL pending independent review.** DIRT now uses
> the LAMMPS SDS pseudo-force dimensions, including an explicit conversion of
> the prior torque-law coefficients. The complete ensemble retains every
> original physical gate; no calibrated μ_r is reported unless it passes them.

### Wall-SDS sign correction (2026-07-18)

The floor is a rolling contact, so its SDS kinematics must be mapped separately
from particle-particle contacts. LAMMPS defines the wall pseudo-velocity as
`vrl = R (omega x n)` using the wall-to-particle normal. DIRT uses that same
normal for walls (but the opposite normal for particle pairs); the previous
wall implementation incorrectly retained the pair-contact minus sign. That
made the floor SDS spring/dashpot inject rolling power in a simple planar
contact. The wall path now uses the LAMMPS wall sign, with a focused regression
that requires negative rolling power for a spinning sphere on a plane.

An isolated, regenerated DIRT `mu_r=0.30`, replicate-0 case then met the
unchanged confined and released rest requirements at steps 19,000 and 35,600,
respectively, retaining all 1,200 particles. The shared fitter measured
**3.8866°**, not the 22--26° target. This validates the sign and removes the
previous non-settling failure; it is adverse calibration evidence, not a full
six-by-two campaign, an experimental comparison, or a pinned `mu_r`.

### Cross-substrate hand-off boundary

This DIRT experiment does **not** currently produce an SPH constitutive input.
`μ_r` is a DEM rolling-contact parameter, while the checked `dev_soil_sph`
glass-bead material interface exposes continuum parameters such as `μ_s`,
`μ_2`, `I_0`, density, and bulk modulus—not a rolling-friction parameter.
Consequently `../calibration.yaml` records `mu_r_pinned: null` and a withheld
status. This prevents an unqualified DEM number from being represented as a
cross-substrate calibration. A future hand-off needs (a) a defined SPH model
term that actually consumes the quantity, (b) a protocol-matched experimental
record, and (c) a cross-substrate validation of the resulting SPH prediction.
The current LAMMPS disagreement is adverse evidence, not a substitute for any
of those requirements.

### Handoff metadata correction (2026-07-18)

The shared calibration header previously described `μ_r=0.10` as a canonical
glass value while its machine-readable handoff correctly contained
`mu_r_pinned: null`. Those statements are incompatible: a reader can copy a
headline value even when no consuming SPH field exists. The header now omits
rolling friction and states the same withheld boundary as the structured field.
The regression rejects either a numerical `μ_r=0.10` claim or loss of the
explicit omission notice. This is provenance protection, not a numerical
calibration, a replacement for an experiment, or a change to the 22--26° gate.
It was implemented and checked with AI assistance; the Crossref lookup checks
bibliographic identity only, not experimental extraction or model transfer.

### Qualification and external-reference boundary

The executable records a deposit only if (1) the confined column first meets the
2 mm/s global-speed criterion and (2) the released heap later meets the 1 cm/s
criterion. It exits non-zero instead of auto-lifting an unsettled column or
treating the 150,000-step cap as rest. It writes a solver-authored JSON witness
with both events, their steps/speeds, and N=1200; the sweep rejects a CSV without
that witness rather than parsing log text. After a complete qualified campaign,
the driver writes a local SHA-256 ledger binding each fitted CSV to its exact
config and witness; `graph` re-hashes those artifacts and rejects mixed, stale,
or edited local evidence. This is provenance protection only, not an external
validation receipt.

The 22–26° band is an acceptance target, not a protocol-matched external
reference receipt in this repository. Zhou et al. (1999), *Physica A* 269,
536–553, doi:10.1016/S0378-4371(99)00183-1 establishes that rolling resistance
can affect simulated pile formation, but it does not validate this monodisperse,
e=0.4 lift protocol or its numerical target. `sweep.py external` replays a
separately generated, solver-authored DIRT release in LAMMPS 22 Jul 2025 Update 4
and records hashes for the source CSV, generated LAMMPS data file, input, and
final dump. The snapshot includes position, radius, linear velocity, and angular
velocity. It also uses an explicit **zero-contact-history release**: DIRT clears
its pair Mindlin/SDS and wall histories at lift because LAMMPS `read_data` cannot
receive those records. The normal calibration never clears its formation history.
This isolates a defined post-lift diagnostic rather than misrepresenting an
unmatched formed-contact state. It is still an implementation check,
not a glass-band pass condition or experimental reference: a future closure
needs a citable glass-bead dataset matching diameter, surface condition,
deposition method, and formation protocol.

For a record to be admitted as that future match, its reported range and every
raw angle must each be anchored to a fragment of the same DOI URL checked by
the audit (for example `https://doi.org/<doi>#table-2-row-1`). Solver CSVs,
generated plots, and free-text table labels are rejected as source locators.
This is an adversarial provenance check, not experimental validation.

Before launching that replay, the adapter independently re-parses the
solver-authored CSV and the generated LAMMPS data file and requires every
position, radius, linear velocity, and angular velocity to survive decimal
serialization within `1e-9` absolute units. This guards the initial-value
transfer; it is not an agreement tolerance and cannot make a DIRT--LAMMPS
comparison pass.

`graph` preserves the original monotonicity, spread, completeness, and 22--26°
gates. A protocol-matched primary record is reported as a transferability
caveat, rather than being made into an additional locally authored numerical
pass criterion. The original calibration contract is therefore neither relaxed
nor silently made impossible by a new evidence-file requirement.
That record must be human-extracted from a primary source and match the frozen
material, vessel, deposition, release, measurement definition, and retained
band. A would-be matched record must include a page/figure/table locator, its
extraction method, and locatable finite numerical observations.  The record's
band must bound those cited observations, so a DOI with legitimate metadata
cannot be paired with a separately invented 22--26° target.  It must also name
the source-reported interval and a separate locator for that interval; the
executable rejects a locally widened or narrowed `band_deg`.  At closure time its DOI
metadata is also checked live against Crossref. A network failure or a metadata
mismatch fails closed. Crossref establishes bibliographic identity only—it does
not certify source extraction, protocol equivalence, material transfer, or the
calibration. The committed Zhou metadata record is deliberately incompatible and has
no numerical observation, so it cannot be repurposed as a target. This is a
scientific provenance requirement, not a solver-derived pass criterion.

Run `python3 sweep.py reference-audit` to independently recheck every committed
reference-record bibliography. As of this revision it verifies the Zhou negative
control only; it does not create a matched experimental reference.

### Cross-code receipt contract (2026-07-16)

`sweep.py external` is a diagnostic-receipt command, not a passing validation
command. It may write a hash-bound replay receipt after both solvers reach their
separately qualified rest states, but it exits non-zero with **CROSS-CODE
VALIDATION: NOT ASSESSED**. No predeclared, independently justified DIRT--LAMMPS
angle-agreement criterion exists for this formation protocol, so process
completion cannot be presented as cross-code agreement. Run `sweep.py
replay-generate --mu-r 0.30 --rep 0`, execute that generated DIRT case, then
run `sweep.py external --mu-r 0.30 --rep 0`. The receipt records the two fitted
observables and the declared zero-history provenance; it does not grade their
difference.

For the independently rerun `mu_r=0.30`, replicate 0 witness, DIRT fitted
0.0000° and LAMMPS fitted 20.5123°. This qualitative discrepancy is retained as
adverse evidence, not converted into a tolerance or a calibration result.

### Replay qualification correction (2026-07-15)

Earlier LAMMPS replay inputs advanced a fixed 150,000 post-lift steps. A fixed
duration is not a comparable endpoint: it did not establish LAMMPS's global
heap-rest speed before fitting, whereas DIRT records only a qualified rest state.
The replay now polls LAMMPS every 2,000 steps and writes/fits a final dump only
when its independently computed global maximum particle speed is below the same
1 cm/s threshold used by DIRT; otherwise it fails. This changes neither the
DIRT protocol, target band, fitting window, material parameters, nor acceptance
tolerances. It supersedes the prior fixed-time cross-code number as a qualified
implementation comparison. AI assistance was used for this correction; it is
software validation, not experimental validation or a calibration closure.

### Normal-force boundary correction (2026-07-15)

The DIRT protocol uses `limit_damping = true`: non-cohesive separating contacts
are repulsive-only, rather than allowing viscous normal damping to become
tensile. The LAMMPS replay had omitted its corresponding `limit_damping`
keyword on both grain--grain and grain--wall models. This revision makes that
choice explicit in the DIRT configuration and in every LAMMPS granular model.
It changes neither the glass band, the rolling coefficients, the estimator, nor
any tolerance. Earlier replay numbers are therefore not a matched
normal-contact comparison and are retained only as superseded diagnostics until
the corrected replay is rerun.

AI assistance was used to implement and document this harness. It does not
substitute for experimental provenance, a cross-code reproduction, or physics
sign-off.

### SDS frame-update threshold correction (2026-07-15)

The external LAMMPS source audit also found a discrete update condition that was
not represented in DIRT: `GranSubModRollingSDS` calls its
projection-and-rescale routine only when `abs(history·n) * k_roll >
1e-10 * mu_r * abs(F_n)`. DIRT had performed that transport on every contact
step. DIRT now applies the same threshold for grain--grain and grain--wall SDS
contacts; a focused wall regression covers the sub-threshold no-update boundary.
This is a source-derived constitutive-integration correction, not a fitted
change to the target, estimator, seeds, material values, case set, or tolerance.
It has not yet been used to claim a calibrated value or cross-code agreement.

### Active Hertz-path correction and independent rerun (2026-07-15)

The preceding threshold rule had been added to DIRT's Hooke-contact SDS path,
but this campaign declares `contact_model = "hertz"`.  The active Hertz path
still transported its rolling history at every step.  It now applies the same
LAMMPS source-derived `EPSILON * mu_roll * |F_n|` guard before a
projection/rescale update.  This is a correction to the active discrete
constitutive update, not a parameter fit and does not change the sweep,
glass band, estimator, seeds, or tolerances.

An isolated precision-double rerun of the qualified `mu_r=0.30`, replicate 0
case reached fill rest at step 11,600 (1.6491e-3 m/s) and heap rest at step
28,600 (8.2082e-3 m/s).  Its common fitter measured **5.3923°**; the resulting
complete-state LAMMPS 22 Jul 2025 Update 4 replay measured **19.9839°**.  This
remains an adverse implementation result, not a calibration, external
reference, or physics validation.  The complete 12-case campaign was not
rerun after this correction, so the calibration gate remains FAIL and no
calibrated `mu_r` is reported. AI assistance was used for implementation and
testing; the stated limits remain material.

### Tangential-contact declaration correction (2026-07-15)

The DIRT configuration previously relied on its plain `history` tangential
default. The campaign now explicitly selects DIRT's
`tangential_model = "mindlin_rescale"`, its corresponding unloading-rescaled
history law.  This is a model-declaration correction, not a fit: the geometry,
seeds, material values, rolling sweep, angle estimator, glass band, replicate
spread limit, and monotonicity slack are unchanged.  A new campaign and a fresh
complete-state LAMMPS replay are required before claiming calibration or
cross-code agreement.  AI assistance was used for this implementation; it is
not experimental validation.

### Cross-code tangential-law mapping correction (2026-07-16)

The DIRT campaign declares `tangential_model = "mindlin_rescale"`, but the
generated LAMMPS replay had used plain `tangential mindlin` for its grain pair
and all three wall contacts. LAMMPS documents `mindlin_rescale` as a distinct
unloading-history law, so the earlier replay was not a constitutively matched
comparison. The replay template now declares `mindlin_rescale` at every one of
those four contacts; `sweep.py lammps-mapping-check` fails if that mapping
regresses. This alters no material parameter, target band, seed, fitted window,
case set, or tolerance. The old DIRT/LAMMPS discrepancy is superseded, not
resolved: a new qualified DIRT case and independently settled LAMMPS replay are
still required, as is a protocol-matched experimental glass reference. AI
assistance was used to locate and implement the correction; it is not a
calibration result.

### Force-history tangential-law correction (2026-07-17)

The campaign previously selected `mindlin_rescale`, the displacement-history
variant. The installed LAMMPS granular documentation describes that variant as
a historical compatibility model that effectively rescales tangential loading
twice during unloading, and recommends `mindlin_rescale/force` instead. The
generated DIRT configurations and all four LAMMPS grain/wall declarations now
use `mindlin_rescale/force`; this is a constitutive-model correction, not a
change to the angle band, replicate count, fitting window, or tolerance.

An isolated zero-history `mu_r=0.30`, replicate-0 replay from the same
solver-written pre-lift state reached the declared rest criteria in both
solvers. The shared estimator returned DIRT **5.0233°** and LAMMPS
**21.0542°**. This is still adverse evidence and does not establish either
cross-code agreement or a glass calibration. The numerical acceptance gate and
the requirement for a protocol-matched external glass record remain fail
closed; `mu_r_pinned` remains null.

### Wall force-history correction (2026-07-17)

The preceding declaration was still incomplete in DIRT: particle--particle
contacts honored `mindlin_rescale/force`, but every wall contact used the
displacement-history Mindlin branch. The floor is an active, load-bearing wall
in this experiment, so this made the stated DIRT/LAMMPS mapping false precisely
where base sliding and rolling resistance shape the heap. `dirt_wall` now
stores the selected wall contact's elastic force plus its previous Hertz contact
radius, applies the same unloading rescale as the granular force-history path,
and reconstructs the capped elastic force after a Coulomb-limited step.

This was checked against the official LAMMPS `pair_granular` documentation on
2026-07-17: it distinguishes accumulated-force `mindlin_rescale/force` from
the displacement-history variant and recommends the force-history form. The
change does not alter the 22--26 degree band, sweep points, seeds, fitting
window, rest criterion, spread limit, or monotonicity slack. It has **not**
been used to claim a pass: no full post-correction ensemble, protocol-matched
experimental glass record, or DIRT/LAMMPS agreement exists yet. AI assistance
was used for the implementation and regression; scientific interpretation and
any eventual calibration closure require independent review.

That rerun has now been completed for the qualified `mu_r=0.30`, replicate 0
case: DIRT reached fill rest at step 12,000 (1.898e-3 m/s) and heap rest at step
27,200 (9.598e-4 m/s). The corrected LAMMPS 22 Jul 2025 Update 4 replay reached
its own rest criterion and the common fitter measured **DIRT 0.0000° versus
LAMMPS 20.5123°**. This remains an adverse cross-code result, not a glass-band
pass, an experimental validation, or a calibrated parameter. A complete current
campaign and protocol-matched primary glass measurement are still absent.

### SDS pseudo-force sign correction (2026-07-17)

An independent source audit of LAMMPS `GranSubModRollingSDS::calculate_forces`
found that its spring--dashpot pseudo-force is
`F_roll = -k_roll xi_roll - gamma_roll v_roll`.  DIRT had accumulated the same
`xi_roll += v_roll dt` history but applied its spring with the opposite sign in
both grain--grain and grain--wall paths.  It now uses the LAMMPS force sign and
the corresponding negative slip-history reconstruction.  This is a
constitutive correction, not a change to the sweep, material values, target
band, estimator, replicate gate, or tolerance. Focused DIRT wall/granular
tests and the analytical estimator pass, but no post-correction campaign or
cross-code agreement is claimed here. AI assistance was used for the code and
source audit; a full qualified ensemble and protocol-matched primary glass
record remain required for calibration closure.

### Independent rescue reproduction (2026-07-16)

An isolated precision-double reproduction regenerated the deterministic
`mu_r=0.30`, replicate-0 case from base seed `20260706`.  DIRT again reached
the qualified fill event at step 12,000 (1.898e-3 m/s) and the released-heap
event at step 27,200 (9.598e-4 m/s).  A separately generated zero-history
release from that same pre-lift state reached DIRT heap rest at step 26,400
(8.602e-3 m/s); the LAMMPS 22 Jul 2025 Update 4 replay then produced the same
fitted comparison, **DIRT 0.0000°; LAMMPS 20.5123°**.  This is a deliberately
adversarial implementation check, not a fitted agreement threshold.  It rules
out closing this calibration from the current run, but does not identify which
remaining constitutive or integration difference causes the discrepancy.

The existing analytic-cone estimator regression (24.000000° recovered) and
LAMMPS tangential-law mapping audit (4/4 `mindlin_rescale` declarations) pass.
They validate the analysis and declared contact-law mapping only; they do not
turn the failed cross-code result into experimental evidence. AI assistance was
used to run and document this reproduction. Human review is still required for
scientific interpretation and any future calibration claim.

### Executable protocol boundary

`config.toml` and the configurations generated by `sweep.py generate` are the
only supported inputs for this gate.  Historical hand-edited `pin_*`, `smoke_*`,
and `test_*` inputs have been removed: they combined a different restitution,
particle population/size distribution, viscous drain, rolling-law coefficients,
and machine-local output path.  Leaving those runnable beside the gate made it
too easy to mistake a result from another physical protocol for calibration
evidence.  They are neither negative controls nor evidence for this campaign.

The remaining executable protocol is deliberately narrow: fixed `mu_p=0.16`,
the declared e=0.4 formation protocol, 1200 monodisperse particles, and the
SDS coefficients stated below.  Changing any of those requires a new
protocol-matched external reference and a new, separately reviewed campaign;
it cannot be used to satisfy this gate retrospectively.

## Physics

A loose column of monodisperse spheres is confined inside a thin cylinder on a
flat floor and allowed to settle. The cylinder is then removed ("lifted") and the
column slumps into a conical heap. The heap stops growing when the surface slope
reaches the angle at which gravity along the slope is balanced by inter-particle
friction — the angle of repose:

```
θ_r = atan(slope of the heap surface)
```

measured by fitting the settled surface height `h(r)` against radial distance `r`
on the straight sloping flank, `θ_r = atan(−dh/dr)`.

The solver writes a snapshot only when the maximum particle speed is below
2×10⁻³ m/s both immediately before release and after collapse.  Using the same
bound on both sides prevents a still-relaxing post-release slope from being
reported as a static angle; it does not alter the monotonicity, spread, or
22–26° acceptance criteria.

There is **no exact θ_r**. It depends on μ_p, rolling friction μ_r, restitution,
polydispersity, and the protocol. The calibration sweeps **μ_r at fixed μ_p**, so
the expected behaviour is:

- θ_r **increases monotonically** with μ_r,
- whether any μ_r lands θ_r in the measured glass band **[22, 26]°** is
  reported as a formation-protocol diagnostic, not a canonical-material closure,
- the heap is **reproducible**: independent random packs, with seeds recorded in
  a manifest, give θ_r with a small but non-zero run-to-run spread.

## Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Young's modulus E | 1.0 × 10⁷ | Pa |
| Poisson's ratio ν | 0.25 | — |
| Restitution e | 0.4 (formation aid) | — |
| Sliding friction μ_p (FIXED) | 0.16 (measured glass) | — |
| Rolling Coulomb cap μ_r (SWEPT) | 0.0, 0.05, 0.10, 0.15, 0.20, 0.30 | — |
| Rolling pseudo-force stiffness k_roll | 2.5 × 10³ | N/m |
| Rolling pseudo-force damping γ_roll | 2.5 × 10⁻¹ | N·s/m |
| Density ρ | 2500 | kg/m³ |
| Radius R | 2.0 | mm |
| Mobile heap particles | 1200 | — |
| Confining-cylinder radius | 25 | mm |
| Reps per μ_r (manifest-recorded seeds) | 2 | — |
| Gravity g_z | −9.81 | m/s² |

E is softened to 10 MPa (a routine DEM practice) so the Rayleigh-criterion
timestep the solver auto-selects (≈ 2.6 × 10⁻⁵ s at R = 2 mm) is large enough
that each heap settles in a few seconds of wall-clock time. **μ_p is fixed** at
the measured glass value while **μ_r is swept**, so θ_r(μ_r) is isolated. Each
(μ_r, rep) case uses a distinct recorded `seed` on the inserter so the two reps
are independent random packs while remaining exactly regenerable (the inserter RNG
is otherwise deterministic).

### Rolling resistance — the `sds` spring–dashpot–slider model

Rolling resistance uses the **`sds`** (spring–dashpot–slider) model, the same one
LAMMPS's `pair_style granular … rolling sds k_roll γ_roll μ_roll` implements. It
uses a length-valued pseudo-displacement and pseudo-force:

```
v_roll = −r_eff·((ω_i−ω_j) × n)
F_roll = −k_roll·ξ_roll − γ_roll·v_roll,  capped at |F_roll| ≤ μ_roll·|F_n|
τ_roll = r_eff · (n × F_roll)
```

where ξ_roll is the accumulated rolling displacement (rescaled on slip), v_roll
is the rolling velocity, and r_eff is the reduced radius (the grain radius at a
wall). DIRT exposes this through `rolling_model = "sds"` with
`rolling_stiffness` (k_roll), `rolling_damping` (γ_roll), and `rolling_friction`
(μ_roll, the Coulomb cap) in `[[dem.materials]]`, and `dirt_wall` applies the same
sds rolling on the floor and confining walls.

**Dimensional correction:** the former torque-law values were 0.01 N·m/rad and
1e-6 N·m·s/rad. For the declared R=2 mm, the LAMMPS/DIRT pseudo-force law requires
division by R², giving **k_roll=2500 N/m** and **γ_roll=0.25 N·s/m**. This follows
`GranSubModRollingSDS::calculate_forces` and LAMMPS's documented pseudo-force
formulation; it was not selected to reach the 22–26° band. The stiffness is also
of the same order as the softened contact tangential stiffness (~2000 N/m). The
Coulomb cap **μ_roll** remains the only swept parameter.

### Base friction from a real frictional floor wall

The heap stands directly on a **frictional plane wall** at z = 0 (normal +z).
`dirt_wall` applies **Mindlin sliding (tangential) friction** on plane walls,
using the material's `friction` coefficient (μ) through `friction_ij` — exactly
the base friction the bottom layer needs so it cannot slide out and pancake the
heap into a thin monolayer. The same swept μ therefore governs both the
particle–particle contacts that set the pile's angle and the particle–floor
contacts that anchor its base.

This replaces an earlier workaround (a frozen rough particle bed standing in for
wall friction, from before `dirt_wall` had tangential friction): no second
material, no `[[group]]`/`[[freeze]]`, no base bed — just one frictional
`[[wall]]` plane. The confining cylinder wall now also carries friction, which is
harmless: it is deactivated at the lift before the heap forms.

## Parameter Sweep

- **Sliding friction μ_p**: 0.16 (FIXED, measured glass)
- **Rolling friction μ_r**: 0.0, 0.05, 0.10, 0.15, 0.20, 0.30 (SWEPT)
- **Replicates**: 2 independent random packs per μ_r. Their per-case inserter
  seeds are derived from a base seed or replayed from a manifest, giving a direct
  run-to-run spread for the reproducibility check without unrecoverable entropy.

In the lift-the-cylinder protocol the heap forms by a column *collapse* on the
frictional floor. With μ_p fixed, the intent is that raising μ_r arrests the
surface grains' rolling and steepens the cone.

### Current evidence status (2 manifest-recorded packs per μ_r)

The earlier table of six-point results was generated before the tangential-law
declaration correction above.  It is superseded and deliberately not retained as
current calibration evidence.  The regenerated 12-case campaign writes its
per-realization results to the ignored local evidence ledger and the checked
`sweep.py graph` command renders the measured-vs-band figure.  Its result must be
read from that command, because an unresolved flank is not converted into 0° and
partial campaigns are rejected.

**Current status: FAIL.** The regenerated precision-double campaign resolves
only 5 of 12 qualified deposits: μ_r=0.15 gives one 3.54° fit, μ_r=0.20 gives
3.41° and 0.00°, and μ_r=0.30 gives 0.00° for both replicates; the remaining
seven deposits have no resolvable flank.  It therefore fails the existing
campaign-integrity gate before monotonicity, spread, or glass-band closure can
be assessed.  The complete-state LAMMPS 22 Jul 2025 Update 4 replay of the
same qualified μ_r=0.30, replicate-0 pre-lift state measures 22.2854°, versus
DIRT's 0.0000°.  That is adverse implementation evidence, not a calibration
or experimental reference.  No calibrated μ_r is written into the campaign
closure.  The [22,26]° band, monotonicity allowance, replicate-spread gate,
geometry, seeds, and particle count remain unchanged.

### Independent reproduction receipt (2026-07-14)

An isolated rerun from the PR head regenerated all 12 cases from base seed
`20260706`, rebuilt the precision-double executable, and required each
solver-written fill/lift/heap-rest witness before fitting. It reproduced the
table above exactly: 8/12 unresolved flanks; resolved angles 3.4757° (`mu_r`
0.15), 1.4149° (0.20), and 0.6242°/5.8053° (0.30). The campaign ledger bound
each input, qualified snapshot, and result by SHA-256. This is a negative
reproduction of the DIRT protocol, not a calibration or experimental closure.

The former top-fed LAMMPS receipt is retained only as historical diagnostic
output; it is not parity evidence. New LAMMPS checks replay the DIRT pre-lift
snapshot with the declared SDS coefficients (`k_roll=2500 N/m`,
`gamma_roll=0.25 N s/m`). A replay is still not experimental validation or a
replacement reference.

The complete-state replay (2026-07-15, `mu_r=0.30`, replicate 0) is a **failed
implementation comparison**: DIRT's shared fitter returned 0.6242° after its
qualified lift, while LAMMPS 22 Jul 2025 Update 4 returned 16.8180° after
receiving the same solver-authored pre-lift positions, radii, linear velocities,
and angular velocities with the same nominal SDS coefficients. The earlier
position-only 17.1476° receipt is superseded because it changed the initial
condition by zeroing motion. The receipt records the source CSV, translated
LAMMPS data file, input, and final dump hashes. This large discrepancy is
evidence to investigate the post-lift contact/wall/integration path; it is
neither a parity pass nor a 22–26° calibration closure. The replay command now
refuses to run unless the DIRT pre-lift CSV, final CSV, and solver-written rest
qualification all exist, and the receipt hashes all three while reporting both
independently fitted angles. That is provenance and adversarial-software
evidence only; it does not make a DIRT-authored formation state into an
experimental reference.

### SDS contact-frame transport audit (2026-07-15)

An independent source-level audit of installed LAMMPS 22 Jul 2025 Update 4
found that its SDS rolling history is projected into a changing contact tangent
plane **and rescaled to retain its norm** (`GranSubMod::rotate_rescale_vec`).
DIRT previously projected that history but discarded the normal component's
magnitude. DIRT now performs the same projection-and-rescale transport for
particle contacts and wall contacts; a focused rotated-normal regression proves
the invariant. This is a constitutive-state correction, not a target,
tolerance, estimator, or parameter change.

A fresh isolated `mu_r=0.30`, replicate-0 run with that correction met DIRT's
existing fill and heap-rest witnesses at steps 12,600 and 26,800 respectively,
but its profile still had no resolvable flank. Replaying its complete
solver-written pre-lift state with the installed LAMMPS binary yielded 19.0475°
after LAMMPS independently met the same 1 cm/s rest threshold. The cross-code
discrepancy therefore remains; this single run is evidence of an unresolved
implementation/constitutive gap, not a calibration result, ensemble result, or
glass-band pass.

### Why μ_r does not close the gate (a physical limit, not a bug)

The `sds` rolling-resistance model **is** correctly implemented and applied — to
grain–grain contacts (`dirt_granular`) **and** to the floor plane
(`dirt_wall::wall_rolling_torque`), and the rolling torque is integrated into the
particle angular velocities. The problem is the *magnitude*: the rolling couple is
capped at `μ_r · |F_n| · r_eff` with `r_eff = R/2 ≈ 1 mm`, so it is an
intrinsically small contribution to the slope-holding torque. The lift-the-cylinder
**collapse is sliding-dominated**: at μ_p = 0.16 the surface grains *slide* out
(they do not roll), and once the slope exceeds ≈ atan(μ_p) ≈ 9° the deposit
avalanches flat. Direct probes confirm this — sweeping μ_r from 0 to 2.0 leaves the
apex unchanged (~10 mm), while raising μ_p (0.16 → 0.6) *does* lift it. The
deposit also spreads to the catch wall, so the angle is set by the collapse
energetics + sliding friction, not by μ_r.

To actually reach 22–26° one would need either a higher effective sliding friction,
a gentler (non-ballistic) deposition that does not mobilize the surface into a
sliding avalanche, or a rolling-resistance regime with a much larger couple than
`μ_r·|F_n|·R/2`. None of those is available within the fixed-μ_p lift-the-cylinder
protocol prescribed for this gate.

### Independent single-contact Coulomb control

Before launching the solver campaign, run:

```bash
python3 physical_feasibility.py
```

This independent, TOML-reading audit reports the Coulomb force balance for an
*isolated one-contact* grain on a slope: `mu_p >= tan(theta)`. The retained
22–26° target gives `tan(22°) = 0.4040`, while the declared fixed `mu_p = 0.16`
does not support that one-contact control. The SDS rolling term is a couple and
supplies no net tangential force in that control. A heap can recruit a
multi-contact force network and geometric interlocking, so this is deliberately
not a global repose bound, a campaign gate, a fitted result, or experimental
validation. It is an adversarial diagnostic printed before each run; a positive
control (`mu_p=0.5`) verifies that it does not simply encode the existing value.

The original gate is unchanged and remains unmet. A future study still needs a
complete campaign and an external, protocol-matched glass comparison; neither
the control nor a band match may be presented as canonical rolling-friction
calibration.

The "lift the cylinder" protocol, per case:
1. **fill** — 1200 mobile spheres are inserted inside a narrow 25 mm cylinder
   (a tall poured-column geometry), resting on the frictional floor wall, and
   settle into a packed column under gravity. After ten consecutive 200-step
   samples with its fastest particle below 2 mm/s, the cylinder wall is
   deactivated by name at runtime (the "lift").
2. **lift** — after ten consecutive 200-step samples below 2 mm/s, the column
   slumps across the frictional floor and relaxes into a
   cone. A wide outer cylinder (70 mm, beyond the heap toe) catches the few
   particles flung out during collapse so the count is conserved; it never
   touches the static heap. When the heap comes to rest (fastest particle
   < 2 mm/s for ten consecutive 200-step samples; a 150k-step cap after lift
   fails rather than fabricating a rest state), `main.rs` dumps every particle's final
   `(x, y, z, radius)`.

## Validation Criteria (the calibration gate)

Before interpreting a solver snapshot, run `python3 sweep.py estimator-check`.
It verifies the independent analysis invariant: the fitter recovers an analytic
24° cone and rejects a one-layer pancake as an unresolved flank. This regression
does not enter the physical pass decision or replace the 12 solver cases below.

| Check | Tolerance | Notes |
|-------|-----------|-------|
| Campaign completeness | exactly 6 μ_r × 2 manifest-recorded packs | no subset, duplicate, malformed, or non-finite row may reach the physical gate |
| θ_r monotonic in μ_r | mean may dip ≤ 2.5° between μ_r steps | stochastic slack |
| θ_r overall increase | θ_r(μ_r,max) > θ_r(μ_r,min) + 1° | rolling raises the heap |
| Glass-band diagnostic | some μ_r lands θ_r in [22°, 26°] | still insufficient to edit canonical glass: formation e differs |
| Reproducibility | per-μ_r std dev ≤ 5° (but > 0) | over the 2 manifest-recorded packs |

`graph` prints the per-μ_r table and a PASS/FAIL, exits non-zero on FAIL. A band
match is labelled a formation-protocol match and still fails canonical closure
until a protocol-matched canonical run and external experimental comparison exist.

## How to Run

Everything is driven by `sweep.py` (run from anywhere). With no argument it runs
all three stages in order.

```bash
# Everything: generate configs → build & run → validate & plot
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py --base-seed 20260706

# Or one stage at a time:
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py generate --base-seed 20260706
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py seed-check
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py start
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py graph
# Independent-code sentinel (requires lmp_serial/lmp; not a closure gate)
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py replay-generate --mu-r 0.30 --rep 0
# Run the generated zero-history DIRT config, then create a diagnostic receipt.
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py external --mu-r 0.30 --rep 0
```

`generate` writes both the per-case TOML files under `sweep/` and the exact seed
table under `data/seed_manifest.csv` by default. To reproduce a campaign exactly,
keep that CSV with the run record and regenerate from it:

```bash
python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py generate \
  --seed-manifest examples/SPH_glass_sphere_calibration/03_angle_of_repose/data/seed_manifest.csv
```

Rerunning `generate` with the same `--base-seed` or with the same
`--seed-manifest` rewrites byte-identical configs in the same checkout. Changing
the base seed produces a new independent campaign, and within one campaign every
replicate still has a distinct per-case seed.

`seed-check` is the lightweight reproducibility gate for the generation workflow:
it compares generated config SHA-256 fingerprints for a same-base rerun, a
manifest replay, and a changed-base campaign, then writes the evidence figure
shown below.

`graph` re-reads `data/repose_sweep.csv`, so you can re-validate and re-plot
without re-running the simulations. The LAMMPS cross-code leg is **opt-in**:
`start` runs only DIRT formation cases, while `external --mu-r 0.30` performs a
separate replay from a qualified DIRT pre-lift state when a supported LAMMPS
binary is installed. The replay is deliberately adversarial and returns a
non-zero `NOT ASSESSED` result even after it writes a qualified receipt; it is
not an overlay or a calibration pass.

### Cross-code replay (LAMMPS) — opt-in, not a closure gate

The bundled LAMMPS input is an executable audit template. Its SDS coefficients
match DIRT numerically; the prior `k_roll=100` mapping was not defensible
because LAMMPS's checked `GranSubModRollingSDS` uses the same length-valued
pseudo-force history and force law. Any future cross-code campaign must retain
equal coefficients and compare formation protocols explicitly.

### Single case (default config)

```bash
cargo run --release --example sphcal_angle_of_repose --no-default-features --features precision-double -- examples/SPH_glass_sphere_calibration/03_angle_of_repose/config.toml
```

This runs the representative case (μ_p = 0.16 fixed, μ_r = 0.15) and writes
`data/repose_results.csv` (the final particle positions).

## Generated evidence

`seed-check` and `graph` generate the reproducibility and measured-result plots
under `plots/`.  They are deliberately untracked: no prior run may be presented
as evidence for a fresh calibration.  A candidate closure must commit a newly
generated, hash-bound ledger and measured-vs-reference graph only after the
complete current campaign passes every retained gate.

## Assumptions

- **3D simulation**, monodisperse spheres (single radius).
- **Hertz–Mindlin** normal/tangential contact with viscoelastic damping (DIRT
  default), plus the `sds` (spring–dashpot–slider) rolling-resistance term
  (k_roll, γ_roll, and the swept μ_r).
- **Restitution e = 0.4 is a formation aid**, not the glass value; no result may
  transfer to canonical e = 0.926 glass without a protocol-matched external
  validation (see top).
- **Softened stiffness** (E = 10 MPa) for a tractable timestep — repose angle is
  governed by friction, not by absolute stiffness.
- **Frictional base from a real floor wall.** The heap stands on a frictional
  `[[wall]]` plane at z = 0; `dirt_wall` applies both Mindlin sliding friction
  (μ_p) and the `sds` rolling resistance (μ_r) on the floor.
- θ_r is fit on the **straight cone flank only**; on these shallow collapse
  deposits the fit can find no resolvable flank, which is recorded as an
  unmeasurable observable and fails the campaign rather than reporting 0°.
- The reference is **empirical** (the measured glass band [22, 26]°); the gate
  checks for a μ_r that lands θ_r in that band, and currently finds none.

## Validation status and authorship

This revision was prepared with AI assistance and independently checked against
the LAMMPS granular SDS source. It does **not** establish a calibrated glass
value. Closure requires (1) a deposition protocol that demonstrably yields a
resolved, arrested cone at fixed μ_p=0.16; (2) a protocol-matched experimental
glass-bead target; and (3) an adversarial DIRT/LAMMPS comparison using identical
SDS coefficients and formation conditions. No numerical result from the
pre-correction campaign is transferable.

## References

1. Y.C. Zhou, B.H. Xu, A.B. Yu, P. Zulli, "Rolling friction in the dynamic
   simulation of sandpile formation", *Physica A* 269 (1999) 536–553.
2. H.P. Zhu, Z.Y. Zhou, R.Y. Yang, A.B. Yu, "Discrete particle simulation of
   particulate systems: A review of major applications and findings",
   *Chemical Engineering Science* 63 (2008) 5728–5770.
3. J.M.N.T. Gray, "Particle segregation in dense granular flows",
   *Annu. Rev. Fluid Mech.* 50 (2018) 407–433 (heap/repose context).
