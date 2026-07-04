//! peri_dem_interop — peridynamics (pond-style) and DEM (dirt) as **two physics
//! tiers on one shared soil substrate**.
//!
//! This is *not* a coupling. There is no coupler, no field/particle transport,
//! no `grass_multi` sub-App. Instead, a single [`App`] schedules two force
//! systems over **one** soil [`Atom`] set, **one** neighbour list, and **one**
//! [`BondStore`]:
//!
//! - **Peridynamics tier (pond-style).** Each material point bonds, in the
//!   reference configuration, to every point within the horizon δ. Those bonds
//!   live in soil's [`BondStore`] (partner addressed by global tag). Every step
//!   the [`peri_force`] system sums the constant-micromodulus bond force
//!   `f = c·s·ŷ` over each point's surviving bonds; a bond whose stretch exceeds
//!   the critical stretch `s₀` **breaks** — and, crucially, is *removed* from the
//!   `BondStore`. The force law and the `s₀ = √(5 G₀ / 9Kδ)` relation are taken
//!   verbatim from POND's `pond_bond` / `pond_atom` (Silling & Askari, Comput.
//!   Struct. 83 (2005) 1526). See `NOTES` in the module tail for why breakage
//!   *removes* the bond here rather than marking it in place.
//!
//! - **DEM tier (dirt).** dirt's real [`HertzMindlinContactPlugin`] runs in the
//!   same `Force` phase over the same neighbour list. For every neighbour pair
//!   it consults the shared `BondStore` via `BondStore::are_excluded` and
//!   **skips any pair that is still peri-bonded** — so intact material feels only
//!   the peridynamic force, never a spurious contact. The instant a bond is
//!   removed by fracture, that pair is no longer excluded and the DEM tier takes
//!   over: the fragments interact through Hertz–Mindlin contact. This is exactly
//!   the mechanism dirt's own bonded-DEM (`dirt_bond`) uses to hand a broken
//!   bond back to contact — here it bridges the peri→DEM transition.
//!
//! # The scenario (closed system — the conservation gate)
//!
//! Two identical brittle bars are laid on the x-axis with a gap larger than the
//! horizon (so their peridynamic families never cross — they are two separate
//! specimens). Bar **A** is launched at the resting bar **B**. On impact the
//! unbonded A–B interface atoms interact through **DEM contact** (they were never
//! peri-bonded), driving a compressive stress wave into each bar. The wave
//! reflects off each free far end as **tension**; the bonds there exceed `s₀` and
//! break — a **spall** fragment separates from each bar (peri fracture; the
//! damage field → 1 on the spall planes). No walls, no gravity, no prescribed
//! motion: **every force is internal and pairwise equal-and-opposite**, so the
//! total mass and total momentum of the whole system are conserved across the
//! entire peri→DEM transition. That invariant is the hard acceptance gate and is
//! checked live by [`report`], which prints `CONSERVATION: PASS/FAIL` at the end.
//!
//! ```bash
//! cargo run --release --example peri_dem_interop \
//!     --no-default-features --features precision-double -- \
//!     examples/peri_dem_interop/config.toml
//! ```

use std::collections::HashMap;

use serde::Deserialize;

use dirt_atom::DemAtom;
use dirt_core::prelude::*;
use dirt_core::{dirt_atom, soil_core, soil_verlet};

// ─────────────────────────────────────────────────────────────────────────────
// Config
// ─────────────────────────────────────────────────────────────────────────────

/// The `[interop]` section: one brittle material, plus a list of bars.
#[derive(Deserialize, Clone, Default)]
struct InteropConfig {
    /// Point spacing Δx [m]. Cell volume per point is Δx³.
    spacing: f64,
    /// Peridynamic horizon δ [m] (typically ≈ 3·spacing).
    horizon: f64,
    /// Young's modulus E [Pa].
    youngs_mod: f64,
    /// Mass density ρ [kg/m³].
    density: f64,
    /// Fracture energy G₀ [J/m²] → s₀ = √(5 G₀ / 9Kδ). Ignored if
    /// `critical_stretch` is given.
    fracture_energy: Option<f64>,
    /// Critical bond stretch s₀ (overrides `fracture_energy`).
    critical_stretch: Option<f64>,
    /// DEM sphere radius as a fraction of the spacing (0.5 → touching lattice).
    #[serde(default = "half")]
    radius_scale: f64,
    /// The bars (each a separate brittle specimen).
    bars: Vec<BarDef>,
}

fn half() -> f64 {
    0.5
}

/// One brittle bar: a lattice block with a uniform initial velocity.
#[derive(Deserialize, Clone)]
struct BarDef {
    /// Reference-box lower corner [m].
    min: [f64; 3],
    /// Reference-box upper corner [m].
    max: [f64; 3],
    /// Uniform initial velocity [m/s].
    #[serde(default)]
    velocity: [f64; 3],
}

/// Resolved peridynamic constants, shared with the force system as a resource.
#[derive(Clone, Default)]
struct PeriParams {
    /// Constant bond micromodulus c = 18K/(πδ⁴) [N/m⁶].
    micromodulus: f64,
    /// Critical bond stretch s₀.
    critical_stretch: f64,
    /// Horizon δ [m].
    horizon: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-point peridynamic data (an example-local soil AtomData column)
// ─────────────────────────────────────────────────────────────────────────────

/// Per-point peridynamic extension, registered as a soil `AtomData` column so it
/// is permuted / migrated *with* the atoms by the substrate (mirrors POND's
/// `PdAtom`). `volume` is `#[forward]` so a ghost point would carry its Vⱼ.
#[derive(AtomData)]
struct PeriPoint {
    /// Cell volume Δx³ represented by the point [m³].
    #[forward]
    volume: Vec<f64>,
    /// Initial bond count (the damage denominator), set when the family is built.
    n0: Vec<f64>,
    /// Local damage: fraction of broken bonds, in [0, 1].
    damage: Vec<f64>,
}

impl Default for PeriPoint {
    fn default() -> Self {
        PeriPoint {
            volume: Vec::new(),
            n0: Vec::new(),
            damage: Vec::new(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Conservation bookkeeping (the acceptance gate)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Clone, Default)]
struct Conservation {
    initialized: bool,
    mass0: f64,
    p0: [f64; 3],
    max_mass_drift: f64,
    max_p_drift: f64,
    total_steps: usize,
    verdict_printed: bool,
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    let mut app = App::new();

    // Infrastructure (input, comm, domain, neighbour lists, groups, run, print).
    app.add_plugins(CorePlugins);

    // DEM tier: registers DemAtom + MaterialTable from [[dem.materials]], the
    // Hertz–Mindlin contact force, rotational dynamics, and translational Verlet.
    app.add_plugins(DemAtomPlugin);
    app.add_plugins(soil_verlet::VelocityVerletPlugin::new());
    app.add_plugins(HertzMindlinContactPlugin);
    app.add_plugins(RotationalDynamicsPlugin);

    // Shared substrate: the BondStore that both tiers ride. dirt's contact force
    // reads it (are_excluded); the peri tier writes/breaks bonds in it.
    app.add_plugins(soil_core::BondPlugin);

    // Peridynamic tier config + constants.
    let cfg = Config::load::<InteropConfig>(&mut app, "interop");
    let k = cfg.youngs_mod / (3.0 * (1.0 - 2.0 * 0.25)); // bond-based PD fixes ν=¼ → K=2E/3
    let c = 18.0 * k / (std::f64::consts::PI * cfg.horizon.powi(4));
    let s0 = match (cfg.critical_stretch, cfg.fracture_energy) {
        (Some(s0), _) => s0,
        (None, Some(g0)) => (5.0 * g0 / (9.0 * k * cfg.horizon)).sqrt(),
        (None, None) => f64::INFINITY,
    };
    app.add_resource(PeriParams {
        micromodulus: c,
        critical_stretch: s0,
        horizon: cfg.horizon,
    });

    // Register the per-point peri column.
    register_atom_data!(app, PeriPoint::default());
    app.add_resource(Conservation::default());

    // Setup: lattice → family → dt.
    app.add_setup_system(
        lattice_insert
            .after("domain_read_input")
            .label("interop_lattice"),
        ScheduleSetupSet::Setup,
    );
    // build_families runs in PostSetup, which the scheduler runs strictly after
    // the Setup phase where the lattice is laid down (no explicit .after needed).
    app.add_setup_system(
        build_families.label("interop_families"),
        ScheduleSetupSet::PostSetup,
    );
    app.add_setup_system(set_timestep, ScheduleSetupSet::PostSetup);

    // Peridynamic force + breakage, in the same Force phase as dirt's contact.
    app.add_update_system(peri_force, ParticleSimScheduleSet::Force);

    // Live conservation + fracture diagnostic.
    app.add_update_system(report, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

// ─────────────────────────────────────────────────────────────────────────────
// Setup systems
// ─────────────────────────────────────────────────────────────────────────────

/// Whether `pos` lies inside this rank's subdomain (single-process: the box).
#[inline]
fn owns_position(domain: &Domain, pos: &[f64; 3]) -> bool {
    (0..3).all(|d| pos[d] >= domain.sub_domain_low[d] && pos[d] < domain.sub_domain_high[d])
}

/// Lay down every bar as a uniform point cloud. Each point is simultaneously a
/// peridynamic point (PeriPoint volume) and a DEM sphere (DemAtom radius) on the
/// same soil Atom — the essence of the shared substrate.
fn lattice_insert(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    cfg: Res<InteropConfig>,
    domain: Res<Domain>,
    comm: Res<CommResource>,
) {
    let mut peri = registry.expect_mut::<PeriPoint>("lattice_insert");
    let mut dem = registry.expect_mut::<DemAtom>("lattice_insert");

    let dx = cfg.spacing;
    assert!(dx > 0.0, "interop.spacing must be > 0");
    let volume = dx * dx * dx;
    let mass = cfg.density * volume;
    let inv_mass = 1.0 / mass;
    let radius = cfg.radius_scale * dx;
    // Solid sphere moment of inertia I = 2/5 m r²  → inv_inertia.
    let inv_inertia = 1.0 / (0.4 * mass * radius * radius);

    let mut next_tag: u32 = atoms.natoms as u32 + 1;
    let mut placed: u64 = 0;

    for bar in &cfg.bars {
        let n: [i64; 3] =
            std::array::from_fn(|d| (((bar.max[d] - bar.min[d]) / dx).round() as i64).max(0));
        for ix in 0..n[0] {
            for iy in 0..n[1] {
                for iz in 0..n[2] {
                    let pos = [
                        bar.min[0] + (ix as f64 + 0.5) * dx,
                        bar.min[1] + (iy as f64 + 0.5) * dx,
                        bar.min[2] + (iz as f64 + 0.5) * dx,
                    ];
                    let tag = next_tag;
                    next_tag += 1;
                    if !owns_position(&domain, &pos) {
                        continue;
                    }

                    // Base soil columns. cutoff_radius = horizon so the one shared
                    // neighbour list captures the entire peri family (a superset
                    // of the shorter-range DEM contact pairs).
                    atoms.natoms += 1;
                    atoms.nlocal += 1;
                    atoms.tag.push(tag);
                    atoms.origin_index.push(0);
                    atoms.cutoff_radius.push(cfg.horizon as Real);
                    atoms.image.push([0, 0, 0]);
                    atoms.is_ghost.push(false);
                    atoms
                        .pos
                        .push([pos[0] as Real, pos[1] as Real, pos[2] as Real]);
                    atoms.vel.push([
                        bar.velocity[0] as Real,
                        bar.velocity[1] as Real,
                        bar.velocity[2] as Real,
                    ]);
                    atoms.force.push([0.0; 3]);
                    atoms.mass.push(mass as Real);
                    atoms.inv_mass.push(inv_mass as Real);
                    atoms.atom_type.push(0); // single material

                    // Peri column.
                    peri.volume.push(volume);
                    peri.n0.push(0.0);
                    peri.damage.push(0.0);

                    // DEM column.
                    dem.radius.push(radius);
                    dem.density.push(cfg.density);
                    dem.inv_inertia.push(inv_inertia);
                    dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
                    dem.omega.push([0.0; 3]);
                    dem.ang_mom.push([0.0; 3]);
                    dem.torque.push([0.0; 3]);
                    dem.body_id.push(0.0); // independent particles (not a rigid clump)

                    placed += 1;
                }
            }
        }
    }

    let total = comm.all_reduce_sum_f64(placed as f64) as u64;
    if comm.rank() == 0 {
        println!(
            "interop: placed {total} dual peri/DEM points across {} bar(s)",
            cfg.bars.len()
        );
    }
}

/// Build the peridynamic family in the reference (initial) configuration: bond
/// every point pair within the horizon. Bars are separated by more than δ, so
/// this produces **no** cross-bar bonds — the bars are independent specimens and
/// their eventual interaction is pure DEM. O(N²) over local points (single
/// process), matching POND's reference-config build.
fn build_families(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    params: Res<PeriParams>,
    comm: Res<CommResource>,
) {
    let mut bonds = registry.expect_mut::<BondStore>("build_families");
    let mut peri = registry.expect_mut::<PeriPoint>("build_families");
    let nlocal = atoms.nlocal as usize;

    while bonds.bonds.len() < nlocal {
        bonds.bonds.push(Vec::new());
    }
    let horizon = params.horizon;
    let mut nbonds = 0u64;

    for i in 0..nlocal {
        for j in (i + 1)..nlocal {
            let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
            let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
            let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist > 0.0 && dist <= horizon {
                bonds.bonds[i].push(BondEntry {
                    partner_tag: atoms.tag[j],
                    bond_type: 0,
                    r0: dist,
                });
                bonds.bonds[j].push(BondEntry {
                    partner_tag: atoms.tag[i],
                    bond_type: 0,
                    r0: dist,
                });
                nbonds += 1;
            }
        }
    }
    for i in 0..nlocal {
        peri.n0[i] = bonds.bonds[i].len() as f64;
    }

    let total = comm.all_reduce_sum_f64(nbonds as f64) as u64;
    if comm.rank() == 0 {
        println!("interop: formed {total} peridynamic bonds (reference configuration)");
    }
}

/// Copy the `[[run]] dt` into `Atom::dt` (the integrator + DEM contact read it).
fn set_timestep(mut atoms: ResMut<Atom>, run_config: Res<RunConfig>) {
    let dt = run_config.current_stage(0).dt;
    if dt > 0.0 {
        atoms.dt = dt;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Peridynamic force + breakage (removal-on-break)
// ─────────────────────────────────────────────────────────────────────────────

/// Bond-based peridynamic force with critical-stretch breakage.
///
/// For each local point `i`, sum `F_i = Σ_j c·s·(y_j−y_i)/|y_j−y_i|·V_i·V_j` over
/// its surviving bonds (constant-micromodulus bond-based PD; POND `pond_bond`).
/// A bond with stretch `s = (|y|−r₀)/r₀ > s₀` **fails**: it is dropped from the
/// point's list (so soil's `are_excluded` immediately stops excluding the pair
/// and dirt's contact tier can engage it). The damage is `1 − n/n₀`.
///
/// Momentum is conserved to round-off: bond `i–j` is evaluated once by `i` (force
/// `+F`) and once by `j` (force `−F`, same `s`, reversed `ŷ`), and both endpoints
/// take the identical break decision from the identical stretch — so a breaking
/// bond contributes zero to *both* endpoints, never a one-sided impulse.
fn peri_force(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>, params: Res<PeriParams>) {
    let mut bonds = registry.expect_mut::<BondStore>("peri_force");
    let mut peri = registry.expect_mut::<PeriPoint>("peri_force");
    let nlocal = atoms.nlocal as usize;
    let ntotal = atoms.len();
    if bonds.bonds.len() < nlocal {
        return;
    }

    // tag → current index (locals + any ghosts), so a partner resolves after a
    // spatial sort or migration.
    let mut tag_to_index: HashMap<u32, usize> = HashMap::with_capacity(ntotal);
    for idx in 0..ntotal {
        tag_to_index.insert(atoms.tag[idx], idx);
    }

    let c = params.micromodulus;
    let s0 = params.critical_stretch;

    for i in 0..nlocal {
        if bonds.bonds[i].is_empty() {
            continue;
        }
        let vi = peri.volume[i];
        let xi = [
            atoms.pos[i][0] as f64,
            atoms.pos[i][1] as f64,
            atoms.pos[i][2] as f64,
        ];

        let mut f = [0.0f64; 3];
        let src = std::mem::take(&mut bonds.bonds[i]);
        let mut kept: Vec<BondEntry> = Vec::with_capacity(src.len());

        for bond in src {
            let j = match tag_to_index.get(&bond.partner_tag) {
                Some(&j) => j,
                None => {
                    kept.push(bond); // partner absent this step; keep the bond
                    continue;
                }
            };
            let dxb = atoms.pos[j][0] as f64 - xi[0];
            let dyb = atoms.pos[j][1] as f64 - xi[1];
            let dzb = atoms.pos[j][2] as f64 - xi[2];
            let ylen = (dxb * dxb + dyb * dyb + dzb * dzb).sqrt();
            if ylen == 0.0 {
                kept.push(bond);
                continue;
            }
            let s = (ylen - bond.r0) / bond.r0;
            if s > s0 {
                // Bond fails: drop it (do NOT keep, apply no force). Both
                // endpoints make this identical decision → symmetric.
                continue;
            }
            let vj = peri.volume[j];
            let coeff = c * s * vi * vj / ylen;
            f[0] += coeff * dxb;
            f[1] += coeff * dyb;
            f[2] += coeff * dzb;
            kept.push(bond);
        }

        atoms.force[i][0] += f[0] as Accum;
        atoms.force[i][1] += f[1] as Accum;
        atoms.force[i][2] += f[2] as Accum;

        let n0 = peri.n0[i];
        peri.damage[i] = if n0 > 0.0 {
            1.0 - (kept.len() as f64) / n0
        } else {
            0.0
        };
        bonds.bonds[i] = kept;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Diagnostic + acceptance gate
// ─────────────────────────────────────────────────────────────────────────────

/// Every `thermo` steps: total mass, total momentum, kinetic energy, peak/summed
/// damage, surviving bonds, and the number of **active DEM contacts** (pairs that
/// overlap and are not peri-excluded — i.e. the DEM tier doing work between
/// fragments). Tracks the running mass/momentum drift and prints the final
/// `CONSERVATION` verdict.
#[allow(clippy::too_many_arguments)]
fn report(
    atoms: Res<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    params: Res<PeriParams>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    comm: Res<CommResource>,
    mut cons: ResMut<Conservation>,
) {
    let nlocal = atoms.nlocal as usize;

    // Conserved quantities.
    let mut mass = 0.0f64;
    let mut p = [0.0f64; 3];
    let mut ke = 0.0f64;
    for i in 0..nlocal {
        let m = atoms.mass[i] as f64;
        let v = [
            atoms.vel[i][0] as f64,
            atoms.vel[i][1] as f64,
            atoms.vel[i][2] as f64,
        ];
        mass += m;
        p[0] += m * v[0];
        p[1] += m * v[1];
        p[2] += m * v[2];
        ke += 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    let mass = comm.all_reduce_sum_f64(mass);
    let p = [
        comm.all_reduce_sum_f64(p[0]),
        comm.all_reduce_sum_f64(p[1]),
        comm.all_reduce_sum_f64(p[2]),
    ];
    let ke = comm.all_reduce_sum_f64(ke);
    let _ = params;

    // Drift tracking — cheap, done every step so the max drift is exact.
    if !cons.initialized {
        cons.initialized = true;
        cons.mass0 = mass;
        cons.p0 = p;
        cons.total_steps = run_config.current_stage(0).steps as usize;
    }
    let mass_drift = (mass - cons.mass0).abs();
    let p_drift =
        ((p[0] - cons.p0[0]).powi(2) + (p[1] - cons.p0[1]).powi(2) + (p[2] - cons.p0[2]).powi(2))
            .sqrt();
    cons.max_mass_drift = cons.max_mass_drift.max(mass_drift);
    cons.max_p_drift = cons.max_p_drift.max(p_drift);

    let step = run_state.total_cycle;
    const REPORT_EVERY: usize = 200;
    let is_final = step + 1 >= cons.total_steps;
    let do_diag = step % REPORT_EVERY == 0 || (is_final && !cons.verdict_printed);
    if !do_diag {
        return;
    }

    // Fracture + DEM-contact diagnostics (only on report steps — the exclusion
    // scan is O(pairs × family) and need not run every step).
    let bonds = registry.expect::<BondStore>("report");
    let peri = registry.expect::<PeriPoint>("report");
    let dem = registry.expect::<DemAtom>("report");
    let mut dmax = 0.0f64;
    let mut dsum = 0.0f64;
    let mut nbonds = 0u64;
    for i in 0..nlocal {
        dmax = dmax.max(peri.damage[i]);
        dsum += peri.damage[i];
        nbonds += bonds.bonds[i].len() as u64;
    }
    // Count active DEM contacts on the shared neighbour list: overlapping pairs
    // that are NOT peri-excluded (precisely the pair set the DEM tier forces). At
    // t=0 every in-horizon pair is peri-bonded → this is 0; it rises only as the
    // bars touch and as bonds break (the peri→DEM handoff, made visible).
    let mut dem_contacts = 0u64;
    for (i, j) in neighbor.pairs(nlocal) {
        if bonds.are_excluded(i, j, &atoms.tag) {
            continue;
        }
        let dxb = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
        let dyb = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
        let dzb = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
        let dist = (dxb * dxb + dyb * dyb + dzb * dzb).sqrt();
        if dist < (dem.radius[i] + dem.radius[j]) {
            dem_contacts += 1;
        }
    }
    let nbonds = comm.all_reduce_sum_f64(nbonds as f64) as u64;
    let dsum = comm.all_reduce_sum_f64(dsum);
    let dem_contacts = comm.all_reduce_sum_f64(dem_contacts as f64) as u64;

    if comm.rank() == 0 && step % REPORT_EVERY == 0 {
        println!(
            "step {step:6}  KE={ke:.4e} J  |p|={:.6e}  bonds={nbonds:7}  \
             dmax={dmax:.3}  Σdmg={dsum:8.2}  DEMcontacts={dem_contacts:5}",
            (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
        );
    }

    // Final verdict (printed once, at/after the last step).
    if comm.rank() == 0 && !cons.verdict_printed && is_final {
        cons.verdict_printed = true;
        // Momentum scale for a relative tolerance: |p0| (bar A's momentum) plus a
        // KE-based floor so a p0=0 symmetric case still has a scale.
        let p0_mag = (cons.p0[0].powi(2) + cons.p0[1].powi(2) + cons.p0[2].powi(2)).sqrt();
        let p_scale = p0_mag.max((2.0 * cons.mass0 * ke).sqrt()).max(1e-30);
        let rel_p = cons.max_p_drift / p_scale;
        let rel_m = if cons.mass0 != 0.0 {
            cons.max_mass_drift / cons.mass0
        } else {
            0.0
        };
        const TOL: f64 = 1e-9;
        let pass = rel_p < TOL && rel_m < TOL;
        println!("──────────────────────────────────────────────────────────────");
        println!("mass0                = {:.6e} kg", cons.mass0);
        println!("|p0|                 = {p0_mag:.6e} kg·m/s");
        println!(
            "max mass drift       = {:.3e}  (rel {rel_m:.3e})",
            cons.max_mass_drift
        );
        println!(
            "max momentum drift   = {:.3e}  (rel {rel_p:.3e})",
            cons.max_p_drift
        );
        println!("tolerance (relative) = {TOL:.1e}");
        println!("surviving bonds      = {nbonds}   peak damage = {dmax:.3}");
        println!("CONSERVATION: {}", if pass { "PASS" } else { "FAIL" });
        println!("──────────────────────────────────────────────────────────────");
    }
}

// NOTES — why breakage *removes* the bond (vs POND marking it in place):
// POND's `pond_bond` flags a failed bond with `bond_type = BOND_BROKEN` and
// counts the flags for damage, because a pure-peridynamics run never needs the
// broken pair to do anything else. Here the *same* pair must transition to the
// DEM tier, and soil's `BondStore::are_excluded` excludes a pair whenever a
// `BondEntry` for the partner exists — regardless of its `bond_type` (dirt uses
// `bond_type` as a genuine bond-type index, so it cannot be overloaded as a
// broken flag). Removing the entry is therefore the correct, substrate-level
// signal that the pair has left the peridynamic tier — and it is exactly what
// dirt's own bonded DEM (`dirt_bond`) does on breakage. Damage is recovered from
// the surviving bond count against `n₀`, so no in-place flag is needed.
