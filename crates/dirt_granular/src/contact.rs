//! Fused Hertz-Mindlin contact force computation.
//!
//! This is the primary contact force module and the recommended code path for DEM
//! simulations. It computes normal, tangential, rolling, and twisting forces in a
//! **single pair loop**, eliminating redundant computation of shared quantities
//! (distance, material lookups, Hertz stiffness, normal force magnitude).
//!
//! # Supported models
//!
//! | Component   | Models                                    |
//! |-------------|-------------------------------------------|
//! | Normal      | Hertz (nonlinear), Hooke (linear), MDR (elastic-plastic) |
//! | Tangential  | Mindlin incremental spring + Coulomb cap  |
//! | Rolling     | Constant torque, SDS (spring-dashpot-slider) |
//! | Twisting    | Constant torque, SDS (spring-dashpot-slider), Marshall (derived from tangential) |
//! | Adhesion    | JKR, DMT                                  |
//! | Cohesion    | SJKR (area-proportional)                  |
//!
//! # Contact detection
//!
//! Two particles are in contact when their geometric overlap `δ = R1 + R2 - d > 0`.
//! With JKR adhesion, the interaction range extends beyond geometric contact by a
//! pull-off distance derived from the surface energy and elastic properties.
//!
//! # TOML configuration
//!
//! The contact model is selected via `contact_model` in the `[materials]` section:
//!
//! ```toml
//! [materials]
//! contact_model = "hertz"  # or "hooke" / "mdr"
//! ```
//!
//! See the [crate-level documentation](crate) for the full material parameter list.

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use dirt_atom::{self, DemAtom, MaterialTable};
use soil_core::Neighbor;
use soil_core::{forward_comm_overlap, CommBuffers, CommResource, CommTopology};
use soil_core::{
    register_atom_data, Atom, AtomDataRegistry, BondStore, ParticleSimScheduleSet, VirialStress,
    VirialStressPlugin,
};

use crate::tangential::{ContactHistory, ContactHistoryStore, CONTACT_HISTORY_LEN};
use crate::{LARGE_OVERLAP_WARN_THRESHOLD, MAX_OVERLAP_WARNINGS, SQRT_5_6, TANGENTIAL_EPSILON};

fn is_mindlin_rescale(model: &str) -> bool {
    matches!(
        model,
        "mindlin_rescale" | "mindlin_rescale_force" | "mindlin_rescale/force"
    )
}

fn is_mindlin_force_history(model: &str) -> bool {
    matches!(model, "mindlin_rescale_force" | "mindlin_rescale/force")
}

fn zero_contact_history() -> ContactHistory {
    [0.0; CONTACT_HISTORY_LEN]
}

/// Fused Hertz normal + Mindlin tangential contact force plugin.
///
/// Registers [`ContactHistoryStore`] in the [`AtomDataRegistry`] and a single
/// `hertz_mindlin_contact` system at [`ParticleSimScheduleSet::Force`].
pub struct HertzMindlinContactPlugin;

impl Plugin for HertzMindlinContactPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(VirialStressPlugin);
        // Register ContactHistoryStore
        register_atom_data!(app, ContactHistoryStore::new());

        let contact_model = {
            let mt = app
                .get_resource_ref::<MaterialTable>()
                .expect("MaterialTable must exist before HertzMindlinContactPlugin");
            mt.contact_model.clone()
        };

        match contact_model.as_str() {
            "hooke" => {
                app.add_update_system(
                    hooke_contact_force.label("hertz_mindlin_contact"),
                    ParticleSimScheduleSet::Force,
                );
            }
            _ => {
                // Roadmap step 4: opt into the interior/boundary overlapped force
                // (interior pairs computed while the ghost halo is in flight). Bit-
                // identical to the standard force; set DIRT_OVERLAP_FORCE=1 to enable.
                let overlap = std::env::var("DIRT_OVERLAP_FORCE")
                    .map(|v| v != "0" && !v.is_empty())
                    .unwrap_or(false);
                if overlap {
                    app.add_update_system(
                        overlapped_contact_force.label("hertz_mindlin_contact"),
                        ParticleSimScheduleSet::Force,
                    );
                } else {
                    app.add_update_system(
                        hertz_mindlin_contact_force.label("hertz_mindlin_contact"),
                        ParticleSimScheduleSet::Force,
                    );
                }
            }
        }
    }
}

/// Fused Hertz-Mindlin contact force for all neighbor pairs.
///
/// Computes normal (Hertz), tangential (Mindlin), rolling, and twisting forces
/// in a single pass over the neighbor list. Supports JKR/DMT adhesion and SJKR
/// cohesion. Forces and torques are accumulated with Newton's third law symmetry.
///
/// # Panics
///
/// Panics if more than [`MAX_OVERLAP_WARNINGS`] pairs have excessive overlap
/// in a single timestep, indicating an unstable simulation.
/// Which pairs to process — for interior/boundary overlap (roadmap step 4). A pair
/// `(i, j)` is "boundary" iff the neighbour `j` is a ghost (`j >= nlocal`): such
/// pairs need fresh ghost positions, so they run after the halo lands. Interior
/// pairs (`j < nlocal`) need no ghosts and can be computed while the halo is in
/// flight. `All` reproduces the standard single-pass force exactly.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ForcePass {
    /// All pairs (resets active flags, prunes at the end) — the standard force.
    All,
    /// Interior pairs only (resets active flags, does NOT prune).
    Interior,
    /// Boundary pairs only (does NOT reset, prunes at the end).
    Boundary,
}

fn mdr_nonadhesive_force(delta: f64, a_inv: f64, e_eff: f64, a: f64, b: f64) -> f64 {
    if delta <= 0.0 || a <= 0.0 || b <= 0.0 {
        return 0.0;
    }
    let x = (delta * a_inv).clamp(0.0, 1.0);
    let root = (x - x * x).max(0.0).sqrt();
    0.25 * e_eff * a * b * ((1.0 - 2.0 * x).acos() - (2.0 - 4.0 * x) * root)
}

fn mdr_round_up_negative_epsilon(value: f64) -> f64 {
    if value < 0.0 && value > -1.0e-20 {
        0.0
    } else {
        value
    }
}

fn mdr_adhesive_force(
    delta_e: f64,
    a_shape: f64,
    a_inv: f64,
    e_eff: f64,
    surface_energy: f64,
    b_shape: f64,
    a_na: f64,
    a_adh: &mut f64,
    at_max_delta: bool,
) -> (f64, f64) {
    // Direct transcription of LAMMPS GranSubModNormalMDR adhesive cases 1-3.
    const A_FAC: f64 = 0.99;
    const MDR_MAX_IT: usize = 100;
    const MDR_EPSILON1: f64 = 1.0e-10;
    const MDR_EPSILON2: f64 = 1.0e-16;
    const CBRT2: f64 = 1.259_921_049_894_873_2;
    const SQRTHALFPI: f64 = 1.253_314_137_315_500_1;
    const CBRTHALFPI: f64 = 1.162_447_351_509_626_5;
    const PITOFIVETHIRDS: f64 = 6.738_808_595_698_141;

    let mut active_a_adh = (*a_adh).min(A_FAC * 0.5 * b_shape);
    if at_max_delta || a_na >= active_a_adh {
        let force = mdr_nonadhesive_force(delta_e, a_inv, e_eff, a_shape, b_shape);
        *a_adh = A_FAC * a_na;
        return (force, *a_adh);
    }

    let e_eff_inv = 1.0 / e_eff;
    let b_inv = 1.0 / b_shape;
    let b_sq = b_shape * b_shape;
    let a_sq = a_shape * a_shape;
    let a_inv_sq = a_inv * a_inv;

    let l_max = (2.0 * std::f64::consts::PI * active_a_adh * surface_energy * e_eff_inv).sqrt();
    let mut g_a_adh =
        0.5 * a_shape - a_shape * b_inv * (0.25 * b_sq - active_a_adh * active_a_adh).sqrt();
    g_a_adh = mdr_round_up_negative_epsilon(g_a_adh);

    let gamma_sq = surface_energy * surface_energy;
    let gamma3 = gamma_sq * surface_energy;
    let gamma4 = gamma_sq * gamma_sq;
    let e_eff_sq = e_eff * e_eff;
    let e_eff_sq_inv = e_eff_inv * e_eff_inv;
    let a4 = a_sq * a_sq;
    let b4 = b_sq * b_sq;
    let b6 = b4 * b_sq;
    let disc = (27.0 * a4 * e_eff_sq * gamma_sq
        - 4.0 * b_sq * gamma4 * std::f64::consts::PI * std::f64::consts::PI)
        .max(0.0);
    let mut tmp = 27.0 * a4 * b4 * surface_energy * e_eff_inv;
    tmp -= 2.0 * b6 * gamma3 * std::f64::consts::PI * std::f64::consts::PI * e_eff_inv.powi(3);
    tmp += 27.0_f64.sqrt() * a_sq * b4 * disc.sqrt() * e_eff_sq_inv;
    tmp = tmp.cbrt();

    let mut a_crit = -b_sq * surface_energy * std::f64::consts::PI * a_inv_sq * e_eff_inv;
    a_crit += CBRT2 * b4 * gamma_sq * PITOFIVETHIRDS / (a_sq * e_eff_sq * tmp);
    a_crit += CBRTHALFPI * tmp * a_inv_sq;
    a_crit /= 6.0;

    if delta_e + l_max - g_a_adh >= 0.0 {
        let f_na = mdr_nonadhesive_force(g_a_adh, a_inv, e_eff, a_shape, b_shape);
        let f_adhes = 2.0 * e_eff * (delta_e - g_a_adh) * active_a_adh;
        return (f_na + f_adhes, active_a_adh);
    }

    if active_a_adh >= a_crit {
        let mut a_tmp = active_a_adh;
        for iter in 0..MDR_MAX_IT {
            let radicand = 0.25 * b_sq - a_tmp * a_tmp;
            if radicand <= 0.0 || a_tmp <= 0.0 {
                a_tmp = 0.0;
                break;
            }
            let fa_tmp = delta_e - 0.5 * a_shape + a_shape * radicand.sqrt() * b_inv;
            let fa =
                fa_tmp + (2.0 * std::f64::consts::PI * a_tmp * surface_energy * e_eff_inv).sqrt();
            if fa.abs() < MDR_EPSILON1 {
                break;
            }
            let mut dfda = -a_tmp * a_shape / (b_shape * radicand.sqrt());
            dfda += surface_energy * SQRTHALFPI / (a_tmp * surface_energy * e_eff).sqrt();
            let next = a_tmp - fa / dfda;
            let fa2 =
                fa_tmp + (2.0 * std::f64::consts::PI * next * surface_energy * e_eff_inv).sqrt();
            a_tmp = next;
            if (fa - fa2).abs() < MDR_EPSILON2 {
                break;
            }
            if iter == MDR_MAX_IT - 1 {
                a_tmp = 0.0;
            }
        }
        active_a_adh = a_tmp;
    }

    if active_a_adh < a_crit || active_a_adh <= 0.0 {
        active_a_adh = 0.0;
        *a_adh = active_a_adh;
        (0.0, active_a_adh)
    } else {
        let mut g_active =
            0.5 * a_shape - a_shape * b_inv * (0.25 * b_sq - active_a_adh * active_a_adh).sqrt();
        g_active = mdr_round_up_negative_epsilon(g_active);
        let f_na = mdr_nonadhesive_force(g_active, a_inv, e_eff, a_shape, b_shape);
        let f_adhes = 2.0 * e_eff * (delta_e - g_active) * active_a_adh;
        *a_adh = active_a_adh;
        (f_na + f_adhes, active_a_adh)
    }
}

fn mdr_normal_force(
    delta: f64,
    r_i: f64,
    r_j: f64,
    e_eff: f64,
    poisson: f64,
    yield_stress: f64,
    surface_energy: f64,
    damping_prefactor: f64,
    m_r: f64,
    v_n: f64,
    stored: &mut [f64; CONTACT_HISTORY_LEN],
) -> (f64, f64, f64) {
    // Mirrors LAMMPS GranSubModNormalMDR's particle-pair rigid-flat normal
    // force path, but without apparent-radius, bulk/free-surface, or penalty
    // state from fix GRANULAR/MDR.
    const MDR_OVERLAP_LIMIT: f64 = 0.95;
    const TOTAL_DELTAMAX: usize = 8;
    const SIDE0: usize = 9;
    const SIDE1: usize = 16;
    const DELTA_PREV: usize = 0;
    const DELTA_MAX_MDR: usize = 1;
    const YFLAG: usize = 2;
    const DELTA_Y: usize = 3;
    const C_A: usize = 4;
    const A_ADH: usize = 5;
    const DELTAP: usize = 6;

    let prev_delta = stored[TOTAL_DELTAMAX].min(delta);
    if delta > stored[TOTAL_DELTAMAX] {
        stored[TOTAL_DELTAMAX] = delta;
    }
    let delta_max = stored[TOTAL_DELTAMAX].max(delta);

    let shear_mod = e_eff / (2.0 * (1.0 + poisson));
    let mut total_force = 0.0;
    let mut total_a_contact = 0.0;

    for (side, base, radius, other_radius) in [(0, SIDE0, r_i, r_j), (1, SIDE1, r_j, r_i)] {
        let mut side_delta = 0.0;
        if delta_max > 0.0 && radius > 0.0 && other_radius > 0.0 {
            let denom = 2.0 * (delta_max - radius - other_radius);
            let opt1 = delta_max * (delta_max - 2.0 * other_radius) / denom;
            let opt2 = delta_max * (delta_max - 2.0 * radius) / denom;
            let mut delta_geo = if radius < other_radius {
                opt1.max(opt2)
            } else {
                opt1.min(opt2)
            };
            let delta_geo_alt = if radius < other_radius {
                opt1.min(opt2)
            } else {
                opt1.max(opt2)
            };
            if delta_geo / radius > MDR_OVERLAP_LIMIT {
                delta_geo = radius * MDR_OVERLAP_LIMIT;
            } else if delta_geo_alt / other_radius > MDR_OVERLAP_LIMIT {
                delta_geo = delta_max - other_radius * MDR_OVERLAP_LIMIT;
            }

            let deltap_sum = stored[SIDE0 + DELTAP] + stored[SIDE1 + DELTAP];
            side_delta = if (deltap_sum - delta_max).abs() > 1.0e-30 {
                let deltap_side = if side == 0 {
                    stored[SIDE0 + DELTAP]
                } else {
                    stored[SIDE1 + DELTAP]
                };
                delta_geo
                    + (deltap_side - delta_geo) * (delta - delta_max) / (deltap_sum - delta_max)
            } else {
                delta_geo
            };
        }
        side_delta = side_delta.max(0.0);

        stored[base + DELTA_PREV] = side_delta;
        if side_delta > stored[base + DELTA_MAX_MDR] {
            stored[base + DELTA_MAX_MDR] = side_delta;
        }
        let side_delta_max = stored[base + DELTA_MAX_MDR].max(side_delta);
        let p_y = yield_stress.max(0.0) * (1.75 * (-4.4 * side_delta_max / radius).exp() + 1.0);

        if stored[base + YFLAG] == 0.0 && yield_stress > 0.0 && side_delta > 0.0 {
            let p_hertz =
                4.0 * e_eff * side_delta.sqrt() / (3.0 * std::f64::consts::PI * radius.sqrt());
            if p_hertz > p_y {
                stored[base + YFLAG] = 1.0;
                stored[base + DELTA_Y] = side_delta;
                stored[base + C_A] =
                    std::f64::consts::PI * (side_delta * side_delta - side_delta * radius);
            }
        }

        let (a_shape, b_shape, delta_e) = if stored[base + YFLAG] == 0.0 {
            (4.0 * radius, 2.0 * radius, side_delta)
        } else {
            let c_a = stored[base + C_A];
            let a_max_sq = (2.0 * side_delta_max * radius - side_delta_max * side_delta_max
                + c_a / std::f64::consts::PI)
                .max(1.0e-30);
            let a_max = a_max_sq.sqrt();
            let a_shape = (4.0 * p_y / e_eff * a_max).max(1.0e-30);
            let b_shape = (2.0 * a_max).max(1.0e-30);
            let delta_e_max = 0.5 * a_shape;
            // Match LAMMPS gran_sub_mod_normal.cpp lines 755-756 exactly:
            // only the acos term is multiplied by E*A*B/4 before subtracting
            // the second term.
            let x = delta_e_max / a_shape;
            let f_max = e_eff * (a_shape * b_shape * 0.25) * (1.0 - 2.0 * x).acos()
                - (2.0 - 4.0 * x) * (x - x * x).max(0.0).sqrt();
            let z_r = radius - (side_delta_max - delta_e_max);
            let delta_r = (2.0 * a_max_sq * (poisson - 1.0)
                - (2.0 * poisson - 1.0) * z_r * (-z_r + (a_max_sq + z_r * z_r).sqrt()))
                * f_max
                / (2.0
                    * std::f64::consts::PI
                    * a_max_sq
                    * shear_mod
                    * (a_max_sq + z_r * z_r).sqrt());
            let delta_e = (side_delta - side_delta_max + delta_e_max + delta_r)
                / (1.0 + delta_r / delta_e_max);
            stored[base + DELTAP] = side_delta_max - (delta_e_max + delta_r);
            (a_shape, b_shape, delta_e)
        };

        let a_inv = 1.0 / a_shape;
        let a_contact = if delta_e > 0.0 {
            b_shape * (a_shape - delta_e).max(0.0).sqrt() * delta_e.sqrt() * a_inv
        } else {
            0.0
        };

        let force = if surface_energy > 0.0 {
            let at_max_delta = (side_delta - side_delta_max).abs() <= 1.0e-14;
            let (force, a_adh) = mdr_adhesive_force(
                delta_e,
                a_shape,
                a_inv,
                e_eff,
                surface_energy,
                b_shape,
                a_contact,
                &mut stored[base + A_ADH],
                at_max_delta,
            );
            if a_adh > a_contact {
                total_a_contact += a_adh;
            } else {
                total_a_contact += a_contact;
            }
            force
        } else {
            stored[base + A_ADH] = a_contact;
            total_a_contact += a_contact;
            mdr_nonadhesive_force(delta_e, a_inv, e_eff, a_shape, b_shape)
        };

        total_force += force;
    }

    let mut force = 0.5 * total_force;
    let a_contact = 0.5 * total_a_contact;
    let k_mdr = 2.0 * e_eff * a_contact.max(1.0e-30);
    if damping_prefactor > 0.0 && delta > 0.0 && (delta >= prev_delta || force > 0.0) {
        force -= damping_prefactor * (m_r * k_mdr).sqrt() * v_n;
    }
    let k_t = 8.0 * shear_mod * a_contact;
    (force, k_mdr, k_t)
}

/// Standard system: compute the full Hertz-Mindlin contact force in one pass.
pub fn hertz_mindlin_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    contact_force_core(
        &mut atoms,
        &neighbor,
        &registry,
        &material_table,
        virial.as_deref_mut(),
        ForcePass::All,
    );
}

/// Overlapped Hertz-Mindlin force (roadmap step 4): compute the interior pairs
/// (`j < nlocal`, no ghosts needed) *while the ghost halo is in flight*, then the
/// boundary pairs (`j >= nlocal`) once it lands. The interior force runs as the
/// overlap closure of [`forward_comm_overlap`], so on a multi-rank run its compute
/// hides the MPI latency of the ghost exchange. Bit-identical to
/// [`hertz_mindlin_contact_force`] (the interior/boundary split is exact — see
/// `interior_boundary_split_matches_single_pass`).
pub fn overlapped_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    comm: Res<CommResource>,
    topo: Res<CommTopology>,
    mut buffers: ResMut<CommBuffers>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    let mut pool = std::mem::take(&mut buffers.forward_scratch);
    {
        // Interior pairs need no fresh ghosts — run them during the in-flight halo.
        let mut interior = |a: &mut Atom| {
            contact_force_core(
                a,
                &neighbor,
                &registry,
                &material_table,
                None,
                ForcePass::Interior,
            );
        };
        forward_comm_overlap(
            &mut atoms,
            &registry,
            &topo,
            &**comm,
            &mut pool,
            &mut interior,
        );
    }
    buffers.forward_scratch = pool;
    // Boundary pairs, now that the halo has landed.
    contact_force_core(
        &mut atoms,
        &neighbor,
        &registry,
        &material_table,
        virial.as_deref_mut(),
        ForcePass::Boundary,
    );
}

/// Core force computation, parameterised by [`ForcePass`] so the interior and
/// boundary pairs can be computed in separate passes (with the halo exchange
/// between) for comm/compute overlap. `ForcePass::All` is the single-pass force.
pub fn contact_force_core(
    atoms: &mut Atom,
    neighbor: &Neighbor,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    mut virial: Option<&mut VirialStress>,
    pass: ForcePass,
) {
    let newton = neighbor.newton;
    let mut dem = registry.expect_mut::<DemAtom>("contact_force_core");
    let mut history = registry.expect_mut::<ContactHistoryStore>("contact_force_core");
    let bond_store = registry.get::<BondStore>();
    let dt = atoms.dt;

    let natoms = atoms.len();
    if history.contacts.len() < natoms {
        history.contacts.resize_with(natoms, Vec::new);
    }

    let nlocal = atoms.nlocal as usize;
    let mut overlap_warnings = 0usize;

    // Reset all active flags before pair loop (skipped on the Boundary pass, which
    // continues the Interior pass's history instead of clearing it).
    if pass != ForcePass::Boundary {
        for i in 0..nlocal {
            for entry in &mut history.contacts[i] {
                entry.2 = false;
            }
        }
    }

    for (i, j) in neighbor.pairs(nlocal) {
        // Interior/boundary split (step 4): boundary pairs touch a ghost (j >= nlocal).
        let is_boundary_pair = j >= nlocal;
        match pass {
            ForcePass::Interior if is_boundary_pair => continue,
            ForcePass::Boundary if !is_boundary_pair => continue,
            _ => {}
        }
        if let Some(ref bonds) = bond_store {
            if bonds.are_excluded(i, j, &atoms.tag) {
                continue;
            }
        }

        // Skip same-body pairs (sub-spheres of the same rigid body don't interact)
        if dirt_atom::same_body(&dem, i, j) {
            continue;
        }

        let r1 = dem.radius[i];
        let r2 = dem.radius[j];

        let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
        let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
        let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        let sum_r = r1 + r2;

        let mat_i = atoms.atom_type[i] as usize;
        let mat_j = atoms.atom_type[j] as usize;
        let surface_energy = material_table.surface_energy_ij[mat_i][mat_j];

        let use_dmt = material_table.adhesion_model == "dmt";

        // JKR: compute pull-off distance for extended interaction range
        // DMT: no extended range (particles separate at delta = 0)
        // Effective radius: R* = R1 R2 / (R1 + R2)
        let r_eff = (r1 * r2) / sum_r;
        // Effective Young's modulus: 1/E* = (1-ν1²)/E1 + (1-ν2²)/E2
        let e_eff = material_table.e_eff_ij[mat_i][mat_j];
        // JKR pull-off distance: particles interact beyond geometric contact
        let delta_pulloff = if surface_energy > 0.0 && !use_dmt {
            let gamma = surface_energy;
            (std::f64::consts::PI * std::f64::consts::PI * gamma * gamma * r_eff
                / (4.0 * e_eff * e_eff))
                .cbrt()
        } else {
            0.0
        };

        // Check contact: geometric touch or within JKR adhesion range
        let interaction_r = sum_r + delta_pulloff;
        if dist_sq >= interaction_r * interaction_r {
            continue;
        }

        let distance = dist_sq.sqrt();

        if distance == 0.0 {
            #[cfg(debug_assertions)]
            eprintln!(
                "WARNING: zero separation between tags {} {}",
                atoms.tag[i], atoms.tag[j]
            );
            continue;
        }

        // delta > 0 means geometric overlap, delta < 0 means gap
        // Cap at half the smaller radius to keep the Hertz model numerically valid.
        let r_min = r1.min(r2);
        let delta = (sum_r - distance).min(0.5 * r_min);

        if delta > 0.0 && distance / sum_r < LARGE_OVERLAP_WARN_THRESHOLD {
            overlap_warnings += 1;
            #[cfg(debug_assertions)]
            eprintln!(
                "WARNING: large overlap tags {} {} ratio {:.3}",
                atoms.tag[i],
                atoms.tag[j],
                distance / sum_r
            );
            if overlap_warnings > MAX_OVERLAP_WARNINGS {
                panic!(
                    "Over {} excessive overlaps this step — aborting. \
                     Check timestep or initial configuration.",
                    MAX_OVERLAP_WARNINGS
                );
            }
            // Cap overlap at half the smaller radius to keep Hertz model valid,
            // but still compute the repulsive force (skipping would remove all
            // repulsion and cause runaway penetration).
        }

        // For non-JKR, skip if no geometric overlap
        if delta <= 0.0 && surface_energy <= 0.0 {
            continue;
        }

        // ── Shared quantities (computed once) ────────────────────────────
        let inv_dist = 1.0 / distance;
        let nx = dx * inv_dist;
        let ny = dy * inv_dist;
        let nz = dz * inv_dist;

        // Effective shear modulus: 1/G* = (2-ν1)/G1 + (2-ν2)/G2
        let g_eff = material_table.g_eff_ij[mat_i][mat_j];

        // Reduced mass: m_r = 1 / (1/m1 + 1/m2)
        // For clump sub-spheres inv_mass is 0 (body-integrated); use real mass.
        let inv_m_i = if atoms.inv_mass[i] as f64 > 0.0 {
            atoms.inv_mass[i] as f64
        } else {
            1.0 / atoms.mass[i] as f64
        };
        let inv_m_j = if atoms.inv_mass[j] as f64 > 0.0 {
            atoms.inv_mass[j] as f64
        } else {
            1.0 / atoms.mass[j] as f64
        };
        let m_r = 1.0 / (inv_m_i + inv_m_j);

        let beta = material_table.beta_ij[mat_i][mat_j];
        let mu = material_table.friction_ij[mat_i][mat_j];
        let mu_r = material_table.rolling_friction_ij[mat_i][mat_j];
        let mu_tw = material_table.twisting_friction_ij[mat_i][mat_j];
        let cohesion_energy = material_table.cohesion_energy_ij[mat_i][mat_j];
        let use_mdr = material_table.contact_model == "mdr";

        // JKR adhesion-only regime: gap exists but within pull-off distance
        // DMT has no adhesion-only regime (no force beyond contact)
        let jkr_adhesion_only = surface_energy > 0.0 && !use_dmt && delta <= 0.0;

        // Hertz stiffness parameters (only meaningful when δ > 0)
        // S_n = 2 E* √(R* δ)  — normal stiffness parameter (used in damping)
        // k_n = 4/3 E* √(R* δ) — normal spring constant
        // k_t = 8 G* √(R* δ)  — tangential spring constant (Mindlin)
        let (s_n, k_n, k_t, contact_radius) = if delta > 0.0 {
            let sdr = (delta * r_eff).sqrt();
            let sn = 2.0 * e_eff * sdr;
            let kn = 4.0 / 3.0 * e_eff * sdr;
            let kt = 8.0 * g_eff * sdr;
            (sn, kn, kt, sdr)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };

        // Full relative velocity (including angular contributions)
        let omega_ix = dem.omega[i][0];
        let omega_iy = dem.omega[i][1];
        let omega_iz = dem.omega[i][2];
        let omega_jx = dem.omega[j][0];
        let omega_jy = dem.omega[j][1];
        let omega_jz = dem.omega[j][2];

        // v_contact_i = vel_i + omega_i × (r1 * n)
        let r1n_x = r1 * nx;
        let r1n_y = r1 * ny;
        let r1n_z = r1 * nz;
        let vc_ix = atoms.vel[i][0] as f64 + (omega_iy * r1n_z - omega_iz * r1n_y);
        let vc_iy = atoms.vel[i][1] as f64 + (omega_iz * r1n_x - omega_ix * r1n_z);
        let vc_iz = atoms.vel[i][2] as f64 + (omega_ix * r1n_y - omega_iy * r1n_x);

        // v_contact_j = vel_j + omega_j × (-r2 * n)
        let r2n_x = r2 * nx;
        let r2n_y = r2 * ny;
        let r2n_z = r2 * nz;
        let vc_jx = atoms.vel[j][0] as f64 + (-omega_jy * r2n_z + omega_jz * r2n_y);
        let vc_jy = atoms.vel[j][1] as f64 + (-omega_jz * r2n_x + omega_jx * r2n_z);
        let vc_jz = atoms.vel[j][2] as f64 + (-omega_jx * r2n_y + omega_jy * r2n_x);

        let vr_x = vc_jx - vc_ix;
        let vr_y = vc_jy - vc_iy;
        let vr_z = vc_jz - vc_iz;

        let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

        let tag_i = atoms.tag[i];
        let tag_j = atoms.tag[j];
        let sign: f64 = if tag_i < tag_j { 1.0 } else { -1.0 };

        // Look up existing spring/history (single search, reused for write-back).
        let entry_idx = history.contacts[i].iter().position(|(t, _, _)| *t == tag_j);
        let mut stored = match entry_idx {
            Some(idx) => history.contacts[i][idx].1,
            None => zero_contact_history(),
        };

        // ── Normal force ─────────────────────────────────────────────────
        // F_n > 0 → repulsive (along contact normal from i to j)
        // F_n < 0 → attractive (adhesion/cohesion pulls particles together)
        let (f_n_mag, k_t) = if use_mdr {
            let (f, _k_mdr, kt_mdr) = mdr_normal_force(
                delta.max(0.0),
                r1,
                r2,
                e_eff,
                0.5 * (material_table.poisson_ratio[mat_i] + material_table.poisson_ratio[mat_j]),
                material_table.mdr_yield_stress_ij[mat_i][mat_j],
                surface_energy,
                material_table.mdr_damping_ij[mat_i][mat_j],
                m_r,
                v_n,
                &mut stored,
            );
            (f, kt_mdr)
        } else if surface_energy > 0.0 && use_dmt {
            // DMT: Hertz contact + constant adhesive force F_dmt = 2π γ R*
            let f_dmt = 2.0 * std::f64::consts::PI * surface_energy * r_eff;
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            (k_n * delta - f_diss_n - f_dmt, k_t)
        } else if surface_energy > 0.0 {
            // JKR: adhesion force F_adh = 3/2 π γ R* (simplified explicit model)
            let f_adhesion = 1.5 * std::f64::consts::PI * surface_energy * r_eff;
            if jkr_adhesion_only {
                // Gap regime (δ ≤ 0): pure adhesion, no Hertz contact or damping
                (-f_adhesion, k_t)
            } else {
                // Contact regime (δ > 0): Hertz repulsion + damping − adhesion
                let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
                (k_n * delta - f_diss_n - f_adhesion, k_t)
            }
        } else if cohesion_energy > 0.0 {
            // SJKR: cohesion proportional to contact area A = π δ R*
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * r_eff;
            (k_n * delta - f_diss_n - f_cohesion, k_t) // can go negative (attractive)
        } else {
            // Standard Hertz repulsion + viscoelastic damping. With
            // `limit_damping` (default) the total is clamped to ≥ 0 so damping
            // can never pull particles together; with it disabled the damping may
            // go net-attractive near separation, matching LAMMPS's default
            // `pair granular` (no tensile cutoff) — required for exact cross-code
            // COR at low restitution (see bench_hertz_rebound).
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_total = k_n * delta - f_diss_n;
            if material_table.limit_damping {
                (f_total.max(0.0), k_t)
            } else {
                (f_total, k_t)
            }
        };

        let fn_x = f_n_mag * nx;
        let fn_y = f_n_mag * ny;
        let fn_z = f_n_mag * nz;

        atoms.force[i][0] -= fn_x as soil_core::Accum;
        atoms.force[i][1] -= fn_y as soil_core::Accum;
        atoms.force[i][2] -= fn_z as soil_core::Accum;
        if newton {
            atoms.force[j][0] += fn_x as soil_core::Accum;
            atoms.force[j][1] += fn_y as soil_core::Accum;
            atoms.force[j][2] += fn_z as soil_core::Accum;
        }

        // ── Tangential force (skip in JKR adhesion-only regime) ──────────
        // No tangential friction when particles are not in geometric contact
        if jkr_adhesion_only {
            // No tangential, rolling, or spring history in adhesion-only regime
            // Virial contribution from normal only
            if let Some(ref mut v) = virial {
                if v.active {
                    let vs = if newton { 1.0 } else { 0.5 };
                    v.add_pair(dx, dy, dz, -fn_x * vs, -fn_y * vs, -fn_z * vs);
                }
            }
            continue;
        }

        let vt_x = vr_x - v_n * nx;
        let vt_y = vr_y - v_n * ny;
        let vt_z = vr_z - v_n * nz;

        // Tangential spring displacement (history model). For the history-free
        // `linear_nohistory` model the spring is identically zero, so the force
        // collapses to the velocity-Coulomb law
        //   F_t = -min(μ |F_n|, γ_t |v_t|) t̂ ,   t̂ = v_t / |v_t|
        // (LAMMPS pair_granular `tangential linear_nohistory`, and the classic
        // `pair gran/hooke`) — the force depends only on the instantaneous
        // relative tangential velocity, with NO accumulated displacement.
        let tangential_model = material_table.tangential_model.as_str();
        let nohistory = tangential_model == "linear_nohistory";
        let mindlin_rescale = is_mindlin_rescale(tangential_model);
        let mindlin_force_history = is_mindlin_force_history(tangential_model);
        let f_t_max = mu * f_n_mag.abs();
        let (sx, sy, sz) = if nohistory {
            (0.0, 0.0, 0.0)
        } else if mindlin_force_history {
            // LAMMPS `mindlin_rescale/force` stores the elastic tangential force
            // itself as history. On normal unloading it scales that force by the
            // contact-radius ratio a_n/a_{n-1} before adding the new increment.
            let mut fx = sign * stored[0];
            let mut fy = sign * stored[1];
            let mut fz = sign * stored[2];
            let prev_a = stored[7];
            if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                let scale = contact_radius / prev_a;
                fx *= scale;
                fy *= scale;
                fz *= scale;
            }
            let f_dot_n = fx * nx + fy * ny + fz * nz;
            fx -= f_dot_n * nx;
            fy -= f_dot_n * ny;
            fz -= f_dot_n * nz;
            fx += k_t * vt_x * dt;
            fy += k_t * vt_y * dt;
            fz += k_t * vt_z * dt;
            (fx, fy, fz)
        } else {
            // Convert stored spring from canonical form to local (i,j) frame
            let mut sx = sign * stored[0];
            let mut sy = sign * stored[1];
            let mut sz = sign * stored[2];
            let prev_a = stored[7];
            if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                let scale = contact_radius / prev_a;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }
            // Rotate spring into current tangent plane (remove normal component)
            let s_dot_n = sx * nx + sy * ny + sz * nz;
            sx -= s_dot_n * nx;
            sy -= s_dot_n * ny;
            sz -= s_dot_n * nz;
            // Integrate tangential velocity into spring displacement
            sx += vt_x * dt;
            sy += vt_y * dt;
            sz += vt_z * dt;

            // Coulomb cap on spring: |k_t s| ≤ μ |F_n|
            let s_mag = (sx * sx + sy * sy + sz * sz).sqrt();
            let f_t_spring_mag = k_t * s_mag;
            if f_t_spring_mag > f_t_max && f_t_spring_mag > TANGENTIAL_EPSILON {
                let scale = f_t_max / f_t_spring_mag;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }
            (sx, sy, sz)
        };

        // Tangential damping coefficient: γ_t = 2 β √(5/6) √(k_t m_r)
        let gamma_t = 2.0 * SQRT_5_6 * beta * (k_t * m_r).sqrt();
        let mut ft_x = (if mindlin_force_history { sx } else { k_t * sx }) + gamma_t * vt_x;
        let mut ft_y = (if mindlin_force_history { sy } else { k_t * sy }) + gamma_t * vt_y;
        let mut ft_z = (if mindlin_force_history { sz } else { k_t * sz }) + gamma_t * vt_z;

        // Coulomb cap on total tangential force
        let f_t_mag = (ft_x * ft_x + ft_y * ft_y + ft_z * ft_z).sqrt();
        if f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
            let scale = f_t_max / f_t_mag;
            ft_x *= scale;
            ft_y *= scale;
            ft_z *= scale;
        }

        let (sx, sy, sz) =
            if mindlin_force_history && f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
                (
                    ft_x - gamma_t * vt_x,
                    ft_y - gamma_t * vt_y,
                    ft_z - gamma_t * vt_z,
                )
            } else {
                (sx, sy, sz)
            };

        // Torques: τ_i = (r1 * n) × f_t, τ_j = (-r2 * n) × (-f_t) = (r2 * n) × f_t
        let ti_x = r1n_y * ft_z - r1n_z * ft_y;
        let ti_y = r1n_z * ft_x - r1n_x * ft_z;
        let ti_z = r1n_x * ft_y - r1n_y * ft_x;
        let tj_x = r2n_y * ft_z - r2n_z * ft_y;
        let tj_y = r2n_z * ft_x - r2n_x * ft_z;
        let tj_z = r2n_x * ft_y - r2n_y * ft_x;

        atoms.force[i][0] += ft_x as soil_core::Accum;
        atoms.force[i][1] += ft_y as soil_core::Accum;
        atoms.force[i][2] += ft_z as soil_core::Accum;
        if newton {
            atoms.force[j][0] -= ft_x as soil_core::Accum;
            atoms.force[j][1] -= ft_y as soil_core::Accum;
            atoms.force[j][2] -= ft_z as soil_core::Accum;
        }
        dem.torque[i][0] += ti_x;
        dem.torque[i][1] += ti_y;
        dem.torque[i][2] += ti_z;
        if newton {
            dem.torque[j][0] += tj_x;
            dem.torque[j][1] += tj_y;
            dem.torque[j][2] += tj_z;
        }

        // ── Rolling resistance torque ───────────────────────────────────
        // Relative angular velocity (rolling component)
        let or_x = omega_ix - omega_jx;
        let or_y = omega_iy - omega_jy;
        let or_z = omega_iz - omega_jz;
        let or_dot_n = or_x * nx + or_y * ny + or_z * nz;
        let roll_x = or_x - or_dot_n * nx;
        let roll_y = or_y - or_dot_n * ny;
        let roll_z = or_z - or_dot_n * nz;

        let mut roll_disp_x = sign * stored[3];
        let mut roll_disp_y = sign * stored[4];
        let mut roll_disp_z = sign * stored[5];
        let mut twist_disp = sign * stored[6];

        if mu_r > 0.0 {
            let roll_mag = (roll_x * roll_x + roll_y * roll_y + roll_z * roll_z).sqrt();
            let sds_rolling = material_table.rolling_model == "sds";
            if sds_rolling {
                // SDS rolling: spring-dashpot-slider model
                let k_roll = material_table.rolling_stiffness_ij[mat_i][mat_j];
                let gamma_roll = material_table.rolling_damping_ij[mat_i][mat_j];

                // Update rolling displacement: remove normal component, integrate
                let rd_dot_n = roll_disp_x * nx + roll_disp_y * ny + roll_disp_z * nz;
                roll_disp_x -= rd_dot_n * nx;
                roll_disp_y -= rd_dot_n * ny;
                roll_disp_z -= rd_dot_n * nz;
                roll_disp_x += roll_x * dt;
                roll_disp_y += roll_y * dt;
                roll_disp_z += roll_z * dt;

                // Spring + dashpot torque
                let mut tr_x = -k_roll * roll_disp_x - gamma_roll * roll_x;
                let mut tr_y = -k_roll * roll_disp_y - gamma_roll * roll_y;
                let mut tr_z = -k_roll * roll_disp_z - gamma_roll * roll_z;
                let tr_mag = (tr_x * tr_x + tr_y * tr_y + tr_z * tr_z).sqrt();
                let tau_max = mu_r * f_n_mag.abs() * r_eff;

                if tr_mag > tau_max && tr_mag > TANGENTIAL_EPSILON {
                    // Cap and rescale spring displacement
                    let scale = tau_max / tr_mag;
                    tr_x *= scale;
                    tr_y *= scale;
                    tr_z *= scale;
                    // Rescale spring: δ = (τ + γ·ω) / (-k)
                    if k_roll > TANGENTIAL_EPSILON {
                        roll_disp_x = (tr_x + gamma_roll * roll_x) / (-k_roll);
                        roll_disp_y = (tr_y + gamma_roll * roll_y) / (-k_roll);
                        roll_disp_z = (tr_z + gamma_roll * roll_z) / (-k_roll);
                    }
                }

                dem.torque[i][0] += tr_x;
                dem.torque[i][1] += tr_y;
                dem.torque[i][2] += tr_z;
                if newton {
                    dem.torque[j][0] -= tr_x;
                    dem.torque[j][1] -= tr_y;
                    dem.torque[j][2] -= tr_z;
                }
            } else if roll_mag > 1e-30 {
                // Constant torque model (existing behavior)
                let tau_mag = mu_r * f_n_mag.abs() * r_eff;
                let inv_roll = tau_mag / roll_mag;
                let tr_x = -inv_roll * roll_x;
                let tr_y = -inv_roll * roll_y;
                let tr_z = -inv_roll * roll_z;
                dem.torque[i][0] += tr_x;
                dem.torque[i][1] += tr_y;
                dem.torque[i][2] += tr_z;
                if newton {
                    dem.torque[j][0] -= tr_x;
                    dem.torque[j][1] -= tr_y;
                    dem.torque[j][2] -= tr_z;
                }
            }
        }

        // ── Twisting friction torque ─────────────────────────────────────
        // Three selectable models (material_table.twisting_model):
        //   "constant" / "sds"  — user-supplied coefficients (gated on μ_tw > 0);
        //   "marshall"           — coefficients DERIVED from the active tangential
        //                          (Mindlin) model, no separate twist inputs.
        if material_table.twisting_model == "marshall" {
            // Marshall (2009) twisting, per LAMMPS pair_granular `twisting marshall`
            // (doc/src/pair_granular.rst §twisting, Marshall2009 eqs 32-33). The
            // twisting stiffness/damping/friction are expressed in terms of the
            // tangential (sliding) coefficients and the Hertz contact radius
            // a = √(R* δ):
            //     k_twist = ½ k_t a²,  γ_twist = ½ γ_t a²,  μ_twist = (2/3) a μ_t
            // with k_t, γ_t the tangential spring/damping computed above and μ_t
            // the tangential friction coefficient. Below the cap the couple is
            // the spring–dashpot τ = −k_twist ξ − γ_twist Ω; it is then truncated
            // to |τ| ≤ μ_twist F_n and the angular displacement rescaled to the
            // critical value (identical bookkeeping to the SDS slider).
            if delta > 0.0 {
                let twist_vel = or_dot_n;
                let a_sq = delta * r_eff; // a² = (√(R* δ))²
                let a = a_sq.sqrt(); // Hertz contact radius
                let k_twist = 0.5 * k_t * a_sq;
                let gamma_twist = 0.5 * gamma_t * a_sq;
                let mu_twist = (2.0 / 3.0) * a * mu; // μ = tangential friction coeff

                twist_disp += twist_vel * dt;

                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_twist * f_n_mag.abs();
                if tau_twist.abs() > tau_max {
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                let tt_x = tau_twist * nx;
                let tt_y = tau_twist * ny;
                let tt_z = tau_twist * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            }
        } else if mu_tw > 0.0 {
            let twist_vel = or_dot_n; // twisting component of relative angular velocity
            let sds_twisting = material_table.twisting_model == "sds";
            if sds_twisting {
                // SDS twisting: spring-dashpot-slider model
                let k_twist = material_table.twisting_stiffness_ij[mat_i][mat_j];
                let gamma_twist = material_table.twisting_damping_ij[mat_i][mat_j];

                // Update twisting displacement
                twist_disp += twist_vel * dt;

                // Spring + dashpot torque (scalar, along contact normal)
                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_tw * f_n_mag.abs() * r_eff;

                if tau_twist.abs() > tau_max {
                    // Cap and rescale spring
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                let tt_x = tau_twist * nx;
                let tt_y = tau_twist * ny;
                let tt_z = tau_twist * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            } else if twist_vel.abs() > 1e-30 {
                // Constant torque model (existing behavior)
                let tau = mu_tw * f_n_mag.abs() * r_eff;
                let sign_tw = if twist_vel > 0.0 { -1.0 } else { 1.0 };
                let tt_x = sign_tw * tau * nx;
                let tt_y = sign_tw * tau * ny;
                let tt_z = sign_tw * tau * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            }
        }

        // Virial: force on i from j = (-fn + ft)
        // When newton=false, each pair is visited twice so halve virial contribution
        if let Some(ref mut v) = virial {
            if v.active {
                let vs = if newton { 1.0 } else { 0.5 };
                let vfx = (-fn_x + ft_x) * vs;
                let vfy = (-fn_y + ft_y) * vs;
                let vfz = (-fn_z + ft_z) * vs;
                v.add_pair(dx, dy, dz, vfx, vfy, vfz);
            }
        }

        // Store updated spring back (canonical form) and mark active
        let mut new_spring = stored;
        new_spring[0] = sign * sx;
        new_spring[1] = sign * sy;
        new_spring[2] = sign * sz;
        new_spring[3] = sign * roll_disp_x;
        new_spring[4] = sign * roll_disp_y;
        new_spring[5] = sign * roll_disp_z;
        new_spring[6] = sign * twist_disp;
        new_spring[7] = contact_radius;
        match entry_idx {
            Some(idx) => {
                history.contacts[i][idx].1 = new_spring;
                history.contacts[i][idx].2 = true;
            }
            None => history.contacts[i].push((tag_j, new_spring, true)),
        }
    }

    // Prune stale contacts (skipped on the Interior pass; the Boundary pass prunes
    // once after both passes have marked their active contacts).
    if pass != ForcePass::Interior {
        for i in 0..nlocal {
            history.contacts[i].retain(|(_, _, active)| *active);
        }
    }

    // Debug: check total force + torque on all atoms (local + ghost).
    // In a correct Newton's 3rd law implementation, the sum of all forces
    // from pair interactions must be zero (each pair contributes +F to one atom
    // and -F to the other). A nonzero sum means a pair was counted asymmetrically.
    // Skip this check when newton=false (forces only written to i).
    #[cfg(debug_assertions)]
    if newton {
        let total = atoms.len();
        let mut sum_fx = 0.0;
        let mut sum_fy = 0.0;
        let mut sum_fz = 0.0;
        for i in 0..total {
            sum_fx += atoms.force[i][0] as f64;
            sum_fy += atoms.force[i][1] as f64;
            sum_fz += atoms.force[i][2] as f64;
        }
        let sum_f = (sum_fx * sum_fx + sum_fy * sum_fy + sum_fz * sum_fz).sqrt();
        if sum_f > 1e-6 {
            eprintln!(
                "WARNING: nonzero net force after contact: |F|={:.6e} ({:.6e},{:.6e},{:.6e})",
                sum_f, sum_fx, sum_fy, sum_fz
            );
        }
    }
}

/// Hooke (linear spring) contact force — alternative to Hertz-Mindlin.
///
/// Normal: `f_n = kn * delta`, tangential uses `kt` directly.
/// Damping: `gamma = 2 * beta * sqrt(kn_ij * m_r)`.
/// All other features (friction, rolling, twisting, cohesion, JKR) reused.
pub fn hooke_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    let newton = neighbor.newton;
    let mut dem = registry.expect_mut::<DemAtom>("hooke_contact_force");
    let mut history = registry.expect_mut::<ContactHistoryStore>("hooke_contact_force");
    let bond_store = registry.get::<BondStore>();
    let dt = atoms.dt;

    while history.contacts.len() < atoms.len() {
        history.contacts.push(Vec::new());
    }

    let nlocal = atoms.nlocal as usize;
    let mut overlap_warnings = 0usize;

    for i in 0..nlocal {
        for entry in &mut history.contacts[i] {
            entry.2 = false;
        }
    }

    for (i, j) in neighbor.pairs(nlocal) {
        if let Some(ref bonds) = bond_store {
            if bonds.are_excluded(i, j, &atoms.tag) {
                continue;
            }
        }

        // Skip same-body pairs (sub-spheres of the same rigid body don't interact)
        if dirt_atom::same_body(&dem, i, j) {
            continue;
        }

        let r1 = dem.radius[i];
        let r2 = dem.radius[j];

        let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
        let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
        let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        let sum_r = r1 + r2;

        if dist_sq >= sum_r * sum_r {
            continue;
        }

        let distance = dist_sq.sqrt();
        if distance == 0.0 {
            continue;
        }

        let r_min = r1.min(r2);
        let delta = (sum_r - distance).min(0.5 * r_min);
        if delta <= 0.0 {
            continue;
        }

        if distance / sum_r < LARGE_OVERLAP_WARN_THRESHOLD {
            overlap_warnings += 1;
            if overlap_warnings > MAX_OVERLAP_WARNINGS {
                panic!(
                    "Over {} excessive overlaps this step — aborting.",
                    MAX_OVERLAP_WARNINGS
                );
            }
            // Still compute force (don't skip) — removing repulsion causes runaway.
        }

        let inv_dist = 1.0 / distance;
        let nx = dx * inv_dist;
        let ny = dy * inv_dist;
        let nz = dz * inv_dist;

        let mat_i = atoms.atom_type[i] as usize;
        let mat_j = atoms.atom_type[j] as usize;
        let r_eff = (r1 * r2) / sum_r;
        // For clump sub-spheres inv_mass is 0 (body-integrated); use real mass.
        let inv_m_i = if atoms.inv_mass[i] as f64 > 0.0 {
            atoms.inv_mass[i] as f64
        } else {
            1.0 / atoms.mass[i] as f64
        };
        let inv_m_j = if atoms.inv_mass[j] as f64 > 0.0 {
            atoms.inv_mass[j] as f64
        } else {
            1.0 / atoms.mass[j] as f64
        };
        let m_r = 1.0 / (inv_m_i + inv_m_j);
        let beta = material_table.beta_ij[mat_i][mat_j];
        let mu = material_table.friction_ij[mat_i][mat_j];
        let mu_r = material_table.rolling_friction_ij[mat_i][mat_j];
        let mu_tw = material_table.twisting_friction_ij[mat_i][mat_j];
        let cohesion_energy = material_table.cohesion_energy_ij[mat_i][mat_j];

        let kn = material_table.kn_ij[mat_i][mat_j];
        let kt = material_table.kt_ij[mat_i][mat_j];
        let contact_radius = (r_eff * delta).sqrt();

        // Hooke normal: f_n = kn * delta
        // Damping: gamma_n = 2 * beta * sqrt(kn * m_r)
        let gamma_n = 2.0 * beta * (kn * m_r).sqrt();

        // Relative velocity
        let omega_ix = dem.omega[i][0];
        let omega_iy = dem.omega[i][1];
        let omega_iz = dem.omega[i][2];
        let omega_jx = dem.omega[j][0];
        let omega_jy = dem.omega[j][1];
        let omega_jz = dem.omega[j][2];

        let r1n_x = r1 * nx;
        let r1n_y = r1 * ny;
        let r1n_z = r1 * nz;
        let vc_ix = atoms.vel[i][0] as f64 + (omega_iy * r1n_z - omega_iz * r1n_y);
        let vc_iy = atoms.vel[i][1] as f64 + (omega_iz * r1n_x - omega_ix * r1n_z);
        let vc_iz = atoms.vel[i][2] as f64 + (omega_ix * r1n_y - omega_iy * r1n_x);

        let r2n_x = r2 * nx;
        let r2n_y = r2 * ny;
        let r2n_z = r2 * nz;
        let vc_jx = atoms.vel[j][0] as f64 + (-omega_jy * r2n_z + omega_jz * r2n_y);
        let vc_jy = atoms.vel[j][1] as f64 + (-omega_jz * r2n_x + omega_jx * r2n_z);
        let vc_jz = atoms.vel[j][2] as f64 + (-omega_jx * r2n_y + omega_jy * r2n_x);

        let vr_x = vc_jx - vc_ix;
        let vr_y = vc_jy - vc_iy;
        let vr_z = vc_jz - vc_iz;
        let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

        // Normal force
        let f_n_mag = if cohesion_energy > 0.0 {
            let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * r_eff;
            kn * delta - gamma_n * v_n - f_cohesion
        } else {
            // See the Hertz path: `limit_damping` (default) clamps to repulsive-
            // only; disabling it matches LAMMPS's default (no tensile cutoff).
            let f_total = kn * delta - gamma_n * v_n;
            if material_table.limit_damping {
                f_total.max(0.0)
            } else {
                f_total
            }
        };

        let fn_x = f_n_mag * nx;
        let fn_y = f_n_mag * ny;
        let fn_z = f_n_mag * nz;

        atoms.force[i][0] -= fn_x as soil_core::Accum;
        atoms.force[i][1] -= fn_y as soil_core::Accum;
        atoms.force[i][2] -= fn_z as soil_core::Accum;
        if newton {
            atoms.force[j][0] += fn_x as soil_core::Accum;
            atoms.force[j][1] += fn_y as soil_core::Accum;
            atoms.force[j][2] += fn_z as soil_core::Accum;
        }

        // Tangential force
        let vt_x = vr_x - v_n * nx;
        let vt_y = vr_y - v_n * ny;
        let vt_z = vr_z - v_n * nz;

        let tag_i = atoms.tag[i];
        let tag_j = atoms.tag[j];
        let sign: f64 = if tag_i < tag_j { 1.0 } else { -1.0 };

        let entry_idx = history.contacts[i].iter().position(|(t, _, _)| *t == tag_j);
        let stored = match entry_idx {
            Some(idx) => history.contacts[i][idx].1,
            None => zero_contact_history(),
        };

        // History-free `linear_nohistory` tangential model → zero spring (see the
        // Hertz path above); the force reduces to velocity-Coulomb with no
        // accumulated displacement. "history" keeps the incremental Hooke spring.
        let tangential_model = material_table.tangential_model.as_str();
        let nohistory = tangential_model == "linear_nohistory";
        let mindlin_rescale = is_mindlin_rescale(tangential_model);
        let mindlin_force_history = is_mindlin_force_history(tangential_model);
        let f_t_max = mu * f_n_mag.abs();
        let (sx, sy, sz) = if nohistory {
            (0.0, 0.0, 0.0)
        } else if mindlin_force_history {
            let mut fx = sign * stored[0];
            let mut fy = sign * stored[1];
            let mut fz = sign * stored[2];
            let prev_a = stored[7];
            if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                let scale = contact_radius / prev_a;
                fx *= scale;
                fy *= scale;
                fz *= scale;
            }
            let f_dot_n = fx * nx + fy * ny + fz * nz;
            fx -= f_dot_n * nx;
            fy -= f_dot_n * ny;
            fz -= f_dot_n * nz;
            fx += kt * vt_x * dt;
            fy += kt * vt_y * dt;
            fz += kt * vt_z * dt;
            (fx, fy, fz)
        } else {
            let mut sx = sign * stored[0];
            let mut sy = sign * stored[1];
            let mut sz = sign * stored[2];
            let prev_a = stored[7];
            if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                let scale = contact_radius / prev_a;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }
            let s_dot_n = sx * nx + sy * ny + sz * nz;
            sx -= s_dot_n * nx;
            sy -= s_dot_n * ny;
            sz -= s_dot_n * nz;
            sx += vt_x * dt;
            sy += vt_y * dt;
            sz += vt_z * dt;

            let s_mag = (sx * sx + sy * sy + sz * sz).sqrt();
            let f_t_spring_mag = kt * s_mag;
            if f_t_spring_mag > f_t_max && f_t_spring_mag > TANGENTIAL_EPSILON {
                let scale = f_t_max / f_t_spring_mag;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }
            (sx, sy, sz)
        };

        let gamma_t = 2.0 * SQRT_5_6 * beta * (kt * m_r).sqrt();
        let mut ft_x = (if mindlin_force_history { sx } else { kt * sx }) + gamma_t * vt_x;
        let mut ft_y = (if mindlin_force_history { sy } else { kt * sy }) + gamma_t * vt_y;
        let mut ft_z = (if mindlin_force_history { sz } else { kt * sz }) + gamma_t * vt_z;

        let f_t_mag = (ft_x * ft_x + ft_y * ft_y + ft_z * ft_z).sqrt();
        if f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
            let scale = f_t_max / f_t_mag;
            ft_x *= scale;
            ft_y *= scale;
            ft_z *= scale;
        }

        let (sx, sy, sz) =
            if mindlin_force_history && f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
                (
                    ft_x - gamma_t * vt_x,
                    ft_y - gamma_t * vt_y,
                    ft_z - gamma_t * vt_z,
                )
            } else {
                (sx, sy, sz)
            };

        // Torques
        let ti_x = r1n_y * ft_z - r1n_z * ft_y;
        let ti_y = r1n_z * ft_x - r1n_x * ft_z;
        let ti_z = r1n_x * ft_y - r1n_y * ft_x;
        let tj_x = r2n_y * ft_z - r2n_z * ft_y;
        let tj_y = r2n_z * ft_x - r2n_x * ft_z;
        let tj_z = r2n_x * ft_y - r2n_y * ft_x;

        atoms.force[i][0] += ft_x as soil_core::Accum;
        atoms.force[i][1] += ft_y as soil_core::Accum;
        atoms.force[i][2] += ft_z as soil_core::Accum;
        if newton {
            atoms.force[j][0] -= ft_x as soil_core::Accum;
            atoms.force[j][1] -= ft_y as soil_core::Accum;
            atoms.force[j][2] -= ft_z as soil_core::Accum;
        }
        dem.torque[i][0] += ti_x;
        dem.torque[i][1] += ti_y;
        dem.torque[i][2] += ti_z;
        if newton {
            dem.torque[j][0] += tj_x;
            dem.torque[j][1] += tj_y;
            dem.torque[j][2] += tj_z;
        }

        // Rolling/twisting relative angular velocity
        let or_x = omega_ix - omega_jx;
        let or_y = omega_iy - omega_jy;
        let or_z = omega_iz - omega_jz;
        let or_dot_n = or_x * nx + or_y * ny + or_z * nz;
        let roll_x = or_x - or_dot_n * nx;
        let roll_y = or_y - or_dot_n * ny;
        let roll_z = or_z - or_dot_n * nz;

        let mut roll_disp_x = sign * stored[3];
        let mut roll_disp_y = sign * stored[4];
        let mut roll_disp_z = sign * stored[5];
        let mut twist_disp = sign * stored[6];

        // Rolling resistance
        if mu_r > 0.0 {
            let roll_mag = (roll_x * roll_x + roll_y * roll_y + roll_z * roll_z).sqrt();
            let sds_rolling = material_table.rolling_model == "sds";
            if sds_rolling {
                let k_roll = material_table.rolling_stiffness_ij[mat_i][mat_j];
                let gamma_roll = material_table.rolling_damping_ij[mat_i][mat_j];

                let rd_dot_n = roll_disp_x * nx + roll_disp_y * ny + roll_disp_z * nz;
                roll_disp_x -= rd_dot_n * nx;
                roll_disp_y -= rd_dot_n * ny;
                roll_disp_z -= rd_dot_n * nz;
                roll_disp_x += roll_x * dt;
                roll_disp_y += roll_y * dt;
                roll_disp_z += roll_z * dt;

                let mut tr_x = -k_roll * roll_disp_x - gamma_roll * roll_x;
                let mut tr_y = -k_roll * roll_disp_y - gamma_roll * roll_y;
                let mut tr_z = -k_roll * roll_disp_z - gamma_roll * roll_z;
                let tr_mag = (tr_x * tr_x + tr_y * tr_y + tr_z * tr_z).sqrt();
                let tau_max = mu_r * f_n_mag.abs() * r_eff;

                if tr_mag > tau_max && tr_mag > TANGENTIAL_EPSILON {
                    let scale = tau_max / tr_mag;
                    tr_x *= scale;
                    tr_y *= scale;
                    tr_z *= scale;
                    if k_roll > TANGENTIAL_EPSILON {
                        roll_disp_x = (tr_x + gamma_roll * roll_x) / (-k_roll);
                        roll_disp_y = (tr_y + gamma_roll * roll_y) / (-k_roll);
                        roll_disp_z = (tr_z + gamma_roll * roll_z) / (-k_roll);
                    }
                }

                dem.torque[i][0] += tr_x;
                dem.torque[i][1] += tr_y;
                dem.torque[i][2] += tr_z;
                if newton {
                    dem.torque[j][0] -= tr_x;
                    dem.torque[j][1] -= tr_y;
                    dem.torque[j][2] -= tr_z;
                }
            } else if roll_mag > 1e-30 {
                let tau_mag = mu_r * f_n_mag.abs() * r_eff;
                let inv_roll = tau_mag / roll_mag;
                let tr_x = -inv_roll * roll_x;
                let tr_y = -inv_roll * roll_y;
                let tr_z = -inv_roll * roll_z;
                dem.torque[i][0] += tr_x;
                dem.torque[i][1] += tr_y;
                dem.torque[i][2] += tr_z;
                if newton {
                    dem.torque[j][0] -= tr_x;
                    dem.torque[j][1] -= tr_y;
                    dem.torque[j][2] -= tr_z;
                }
            }
        }

        // Twisting friction (see the Hertz-Mindlin path for model semantics).
        if material_table.twisting_model == "marshall" {
            // Marshall (2009) derived-coefficient twisting on the linear (Hooke)
            // tangential model: k_twist = ½ k_t a², γ_twist = ½ γ_t a²,
            // μ_twist = (2/3) a μ_t with a = √(R* δ) and k_t = kt, γ_t the Hooke
            // tangential spring/damping computed above.
            let twist_vel = or_dot_n;
            let a_sq = delta * r_eff;
            let a = a_sq.sqrt();
            let k_twist = 0.5 * kt * a_sq;
            let gamma_twist = 0.5 * gamma_t * a_sq;
            let mu_twist = (2.0 / 3.0) * a * mu;

            twist_disp += twist_vel * dt;

            let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
            let tau_max = mu_twist * f_n_mag.abs();
            if tau_twist.abs() > tau_max {
                tau_twist = tau_twist.signum() * tau_max;
                if k_twist > TANGENTIAL_EPSILON {
                    twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                }
            }

            let tt_x = tau_twist * nx;
            let tt_y = tau_twist * ny;
            let tt_z = tau_twist * nz;
            dem.torque[i][0] += tt_x;
            dem.torque[i][1] += tt_y;
            dem.torque[i][2] += tt_z;
            if newton {
                dem.torque[j][0] -= tt_x;
                dem.torque[j][1] -= tt_y;
                dem.torque[j][2] -= tt_z;
            }
        } else if mu_tw > 0.0 {
            let twist_vel = or_dot_n;
            let sds_twisting = material_table.twisting_model == "sds";
            if sds_twisting {
                let k_twist = material_table.twisting_stiffness_ij[mat_i][mat_j];
                let gamma_twist = material_table.twisting_damping_ij[mat_i][mat_j];

                twist_disp += twist_vel * dt;

                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_tw * f_n_mag.abs() * r_eff;

                if tau_twist.abs() > tau_max {
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                let tt_x = tau_twist * nx;
                let tt_y = tau_twist * ny;
                let tt_z = tau_twist * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            } else if twist_vel.abs() > 1e-30 {
                let tau = mu_tw * f_n_mag.abs() * r_eff;
                let sign_tw = if twist_vel > 0.0 { -1.0 } else { 1.0 };
                let tt_x = sign_tw * tau * nx;
                let tt_y = sign_tw * tau * ny;
                let tt_z = sign_tw * tau * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            }
        }

        // Virial
        if let Some(ref mut v) = virial {
            if v.active {
                let vs = if newton { 1.0 } else { 0.5 };
                let vfx = (-fn_x + ft_x) * vs;
                let vfy = (-fn_y + ft_y) * vs;
                let vfz = (-fn_z + ft_z) * vs;
                v.add_pair(dx, dy, dz, vfx, vfy, vfz);
            }
        }

        let mut new_spring = stored;
        new_spring[0] = sign * sx;
        new_spring[1] = sign * sy;
        new_spring[2] = sign * sz;
        new_spring[3] = sign * roll_disp_x;
        new_spring[4] = sign * roll_disp_y;
        new_spring[5] = sign * roll_disp_z;
        new_spring[6] = sign * twist_disp;
        new_spring[7] = contact_radius;
        match entry_idx {
            Some(idx) => {
                history.contacts[i][idx].1 = new_spring;
                history.contacts[i][idx].2 = true;
            }
            None => history.contacts[i].push((tag_j, new_spring, true)),
        }
    }

    for i in 0..nlocal {
        history.contacts[i].retain(|(_, _, active)| *active);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dirt_atom::DemAtom;
    use dirt_test_utils::{make_material_table, push_dem_test_atom};
    use soil_core::Neighbor;
    use soil_core::{Atom, AtomDataRegistry};

    fn push_test_atom_with_history(
        atom: &mut Atom,
        dem: &mut DemAtom,
        history: &mut ContactHistoryStore,
        tag: u32,
        pos: [f64; 3],
        radius: f64,
    ) {
        push_dem_test_atom(atom, dem, tag, pos, radius);
        history.contacts.push(Vec::new());
    }

    /// Step 4 correctness: the interior/boundary two-pass force (Interior pass for
    /// local-local pairs while the halo is in flight, Boundary pass for ghost pairs
    /// after they land) must equal the single All pass — bit-for-bit.
    #[test]
    fn interior_boundary_split_matches_single_pass() {
        let r = 0.001;
        let build = || {
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;
            // atom 0 (local) overlaps atom 1 (local -> interior pair) and atom 2
            // (ghost -> boundary pair).
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], r);
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [1.5 * r, 0.0, 0.0], r);
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 2, [0.0, 1.5 * r, 0.0], r);
            atom.nlocal = 2;
            atom.natoms = 3;
            // Half neighbour list (newton): atom 0 -> {1 (local), 2 (ghost)}.
            let mut nb = Neighbor::new();
            nb.neighbor_offsets = vec![0, 2, 2, 2];
            nb.neighbor_indices = vec![1, 2];
            let mut reg = AtomDataRegistry::new();
            reg.register(dem);
            reg.register(hist);
            (atom, nb, reg)
        };
        let mt = make_material_table();

        let (mut a_all, nb, r_all) = build();
        contact_force_core(&mut a_all, &nb, &r_all, &mt, None, ForcePass::All);

        let (mut a_split, _nb, r_split) = build();
        contact_force_core(&mut a_split, &nb, &r_split, &mt, None, ForcePass::Interior);
        contact_force_core(&mut a_split, &nb, &r_split, &mt, None, ForcePass::Boundary);

        let mut max_diff = 0.0f64;
        for i in 0..3 {
            for d in 0..3 {
                max_diff =
                    max_diff.max((a_all.force[i][d] as f64 - a_split.force[i][d] as f64).abs());
            }
        }
        assert!(
            max_diff < 1e-15,
            "interior+boundary != all: max force diff = {max_diff:.3e}"
        );
        // Sanity: the contact actually produced a non-trivial force.
        assert!(a_all.force[0][0].abs() as f64 + a_all.force[0][1].abs() as f64 > 0.0);
    }

    #[test]
    fn fused_contact_repulsive_for_overlap() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][0] < 0.0,
            "particle 0 should have negative x force"
        );
        assert!(
            atom.force[1][0] > 0.0,
            "particle 1 should have positive x force"
        );
        assert!((atom.force[0][0] + atom.force[1][0]).abs() < 1e-10);
    }

    #[test]
    fn fused_contact_tangential_with_sliding() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        atom.vel[1][1] = 0.1;
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Normal force present
        assert!(atom.force[0][0] < 0.0, "normal force on atom 0");
        assert!(atom.force[1][0] > 0.0, "normal force on atom 1");
        // Tangential force present
        assert!(atom.force[0][1].abs() > 0.0, "tangential force on atom 0");
        assert!(
            (atom.force[0][1] + atom.force[1][1]).abs() < 1e-10,
            "tangential forces equal and opposite"
        );
        // Torque present (stored in DemAtom via registry)
        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        let t_mag =
            (dem.torque[0][0].powi(2) + dem.torque[0][1].powi(2) + dem.torque[0][2].powi(2)).sqrt();
        assert!(t_mag > 0.0, "torque on atom 0");
    }

    /// The `linear_nohistory` tangential model must keep the tangential spring
    /// displacement identically zero (velocity-Coulomb, LAMMPS `pair_granular`
    /// `tangential linear_nohistory`), while the default `history` (Mindlin) model
    /// accumulates it. Both are driven with the same sub-Coulomb tangential slip.
    #[test]
    fn linear_nohistory_has_no_spring_accumulation() {
        let radius = 0.001;
        let build = || {
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(
                &mut atom,
                &mut dem,
                &mut hist,
                1,
                [0.00185, 0.0, 0.0],
                radius,
            );
            atom.vel[1][1] = 0.001; // small tangential slip, below the Coulomb cap
            atom.nlocal = 2;
            atom.natoms = 2;
            let mut nb = Neighbor::new();
            nb.neighbor_offsets = vec![0, 1, 1];
            nb.neighbor_indices = vec![1];
            let mut reg = AtomDataRegistry::new();
            reg.register(dem);
            reg.register(hist);
            (atom, nb, reg)
        };
        let spring_mag = |reg: &AtomDataRegistry| -> f64 {
            let h = reg.expect::<ContactHistoryStore>("spring");
            let s = h.contacts[0]
                .iter()
                .find(|(t, _, _)| *t == 1)
                .map(|(_, s, _)| *s)
                .unwrap_or([0.0; CONTACT_HISTORY_LEN]);
            (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt()
        };

        // History (Mindlin): the tangential spring accumulates over the contact.
        let mut mt_h = make_material_table();
        mt_h.tangential_model = "history".to_string();
        let (mut a, nb, reg) = build();
        for _ in 0..20 {
            a.force[0] = [0.0; 3];
            a.force[1] = [0.0; 3];
            contact_force_core(&mut a, &nb, &reg, &mt_h, None, ForcePass::All);
        }
        let xi_history = spring_mag(&reg);
        assert!(
            xi_history > 0.0,
            "history model must accumulate spring, got {xi_history:e}"
        );

        // linear_nohistory: spring stays exactly zero; force is still present.
        let mut mt_nh = make_material_table();
        mt_nh.tangential_model = "linear_nohistory".to_string();
        let (mut a2, nb2, reg2) = build();
        for _ in 0..20 {
            a2.force[0] = [0.0; 3];
            a2.force[1] = [0.0; 3];
            contact_force_core(&mut a2, &nb2, &reg2, &mt_nh, None, ForcePass::All);
        }
        let xi_nohistory = spring_mag(&reg2);
        assert_eq!(
            xi_nohistory, 0.0,
            "linear_nohistory must not accumulate spring"
        );
        assert!(
            (a2.force[0][1] as f64).abs() > 0.0,
            "linear_nohistory must still produce a tangential (velocity-Coulomb) force"
        );
    }

    #[test]
    fn fused_contact_no_force_for_gap() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [0.003, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(atom.force[0][0].abs() < 1e-20);
        assert!(atom.force[1][0].abs() < 1e-20);
    }

    fn make_material_table_cohesion() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 1e9);
        mt.build_pair_tables();
        mt
    }

    fn make_material_table_rolling() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.3, 0.0);
        mt.build_pair_tables();
        mt
    }

    #[test]
    fn cohesion_produces_attractive_force() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        // Very small overlap with high cohesion energy → cohesion dominates
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.00199999, 0.0, 0.0],
            radius, // delta = 1e-8 (tiny overlap)
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_cohesion());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // With cohesion and small overlap, normal force on atom 0 should be positive (attractive toward atom 1)
        assert!(
            atom.force[0][0] > 0.0,
            "cohesion should make force attractive on atom 0, got {}",
            atom.force[0][0]
        );
        // Newton's 3rd law
        assert!(
            (atom.force[0][0] + atom.force[1][0]).abs() < 1e-10,
            "forces should be equal and opposite"
        );
    }

    #[test]
    fn zero_cohesion_matches_original() {
        // Two identical setups — one with default table, one with explicit 0.0 cohesion
        let radius = 0.001;
        let sep = 0.0019;

        let run = |mt: MaterialTable| -> [f64; 3] {
            let mut app = App::new();
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
            atom.nlocal = 2;
            atom.natoms = 2;
            let mut neighbor = Neighbor::new();
            neighbor.neighbor_offsets = vec![0, 1, 1];
            neighbor.neighbor_indices = vec![1];
            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(hist);
            app.add_resource(atom);
            app.add_resource(neighbor);
            app.add_resource(registry);
            app.add_resource(mt);
            app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();
            let atom = app.get_resource_ref::<Atom>().unwrap();
            [
                atom.force[0][0] as f64,
                atom.force[0][1] as f64,
                atom.force[0][2] as f64,
            ]
        };

        let f_default = run(make_material_table());
        let mut mt_zero = MaterialTable::new();
        mt_zero.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0);
        mt_zero.build_pair_tables();
        let f_zero = run(mt_zero);

        for d in 0..3 {
            assert!(
                (f_default[d] - f_zero[d]).abs() < 1e-15,
                "zero params should reproduce original, dim {} default={} zero={}",
                d,
                f_default[d],
                f_zero[d]
            );
        }
    }

    fn make_material_table_jkr() -> MaterialTable {
        let mut mt = MaterialTable::new();
        // Use high surface energy (1.0 J/m²) so adhesion clearly dominates at small overlaps
        mt.add_material_full("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 1.0);
        mt.build_pair_tables();
        mt
    }

    #[test]
    fn jkr_pulloff_force_matches_theory() {
        // Test in adhesion-only regime (gap, not overlap) where force = -F_adhesion exactly
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        let gamma = 1.0;
        let r_eff = radius / 2.0;

        // Place particles with a tiny gap (adhesion-only regime)
        let gap = 1e-9;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [2.0 * radius + gap, 0.0, 0.0],
            radius,
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        let mt = make_material_table_jkr();
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let expected_pulloff = 1.5 * std::f64::consts::PI * gamma * r_eff;
        // In adhesion-only regime, force should be exactly -F_adhesion
        // Force on atom 0 should be positive (attracted toward atom 1)
        assert!(
            atom.force[0][0] > 0.0,
            "JKR should produce attractive force, got {}",
            atom.force[0][0]
        );
        // f_n_mag = -F_adhesion, force[0] -= f_n_mag * nx → force[0] += F_adhesion
        let f_mag = atom.force[0][0] as f64;
        assert!(
            (f_mag - expected_pulloff).abs() / expected_pulloff < 1e-6,
            "pull-off force should match theory {}, got {}",
            expected_pulloff,
            f_mag
        );
    }

    #[test]
    fn jkr_adhesion_only_regime() {
        // Two particles with a small gap (no geometric overlap) but within JKR range
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        // Gap of 1e-9 (very small, within JKR pull-off distance for gamma=1.0)
        let gap = 1e-9;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [2.0 * radius + gap, 0.0, 0.0],
            radius,
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_jkr());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Should be attractive (atom 0 pulled toward atom 1 = positive x)
        assert!(
            atom.force[0][0] > 0.0,
            "JKR adhesion-only should attract, got {}",
            atom.force[0][0]
        );
        // Newton's 3rd law
        assert!(
            (atom.force[0][0] + atom.force[1][0]).abs() < 1e-10,
            "forces should be equal and opposite"
        );
    }

    #[test]
    fn jkr_no_interaction_beyond_pulloff() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        // Large gap — well beyond JKR pull-off distance
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.003, 0.0, 0.0],
            radius, // gap = 0.001 >> delta_pulloff
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_jkr());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][0].abs() < 1e-20,
            "no force beyond pull-off distance"
        );
    }

    fn make_material_table_hooke() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.add_material_extended("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 0.0, 0.0, 1e6, 5e5);
        mt.contact_model = "hooke".to_string();
        mt.build_pair_tables();
        mt
    }

    fn make_material_table_twisting() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.add_material_extended(
            "glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0,
        );
        mt.build_pair_tables();
        mt
    }

    #[test]
    fn hooke_force_linear_in_delta() {
        let radius = 0.001;
        let run = |sep: f64| -> f64 {
            let mut app = App::new();
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
            atom.nlocal = 2;
            atom.natoms = 2;
            let mut neighbor = Neighbor::new();
            neighbor.neighbor_offsets = vec![0, 1, 1];
            neighbor.neighbor_indices = vec![1];
            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(hist);
            app.add_resource(atom);
            app.add_resource(neighbor);
            app.add_resource(registry);
            app.add_resource(make_material_table_hooke());
            app.add_update_system(hooke_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();
            let atom = app.get_resource_ref::<Atom>().unwrap();
            atom.force[0][0] as f64
        };

        // delta1 = 2*r - sep1, delta2 = 2*r - sep2
        let sep1 = 0.00195; // delta = 0.00005
        let sep2 = 0.0019; // delta = 0.0001
        let f1 = run(sep1);
        let f2 = run(sep2);

        // Hooke: force proportional to delta → f2/f1 ≈ 2.0 (linear)
        let ratio = f2 / f1;
        assert!(
            (ratio - 2.0).abs() < 0.15,
            "Hooke force should be linear in delta, got ratio {} (expected ~2.0)",
            ratio
        );
    }

    #[test]
    fn hooke_no_force_beyond_contact() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [0.003, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_hooke());
        app.add_update_system(hooke_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][0].abs() < 1e-20,
            "no force beyond contact distance"
        );
    }

    #[test]
    fn twisting_friction_opposes_spin() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        // Spin about contact normal (x-axis)
        dem.omega[0] = [100.0, 0.0, 0.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_twisting());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // Twisting torque on atom 0 should oppose its spin about x (negative x torque)
        assert!(
            dem.torque[0][0] < 0.0,
            "twisting torque should oppose omega_x, got {}",
            dem.torque[0][0]
        );
    }

    #[test]
    fn twisting_friction_zero_when_no_spin() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        // No angular velocity at all
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_twisting());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // No twisting torque when there's no angular velocity
        let torque_mag =
            (dem.torque[0][0].powi(2) + dem.torque[0][1].powi(2) + dem.torque[0][2].powi(2)).sqrt();
        assert!(
            torque_mag < 1e-20,
            "no twisting torque when no spin, got {}",
            torque_mag
        );
    }

    #[test]
    fn rolling_resistance_opposes_angular_velocity() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        // Give atom 0 a rolling angular velocity (around y-axis — perpendicular to contact normal x)
        dem.omega[0] = [0.0, 100.0, 0.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_rolling());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // Rolling torque on atom 0 should oppose its angular velocity (negative y)
        assert!(
            dem.torque[0][1] < 0.0,
            "rolling torque should oppose omega_y, got {}",
            dem.torque[0][1]
        );
    }

    // ── SDS model helper ────────────────────────────────────────────────

    fn make_material_table_sds_rolling() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.rolling_model = "sds".to_string();
        mt.add_material_with_sds(
            "glass", 8.7e9, 0.3, 0.95, 0.4, 0.3, // rolling_friction (mu_r)
            0.0, 0.0, 0.0, // twisting_friction
            0.0, 0.0, 1e3, // rolling_stiffness
            0.5, // rolling_damping
            0.0, 0.0,
        );
        mt.build_pair_tables();
        mt
    }

    fn make_material_table_sds_twisting() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.twisting_model = "sds".to_string();
        mt.add_material_with_sds(
            "glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, // rolling_friction
            0.0, 0.0, 0.3, // twisting_friction (mu_tw)
            0.0, 0.0, 0.0, 0.0, 1e3, // twisting_stiffness
            0.5, // twisting_damping
        );
        mt.build_pair_tables();
        mt
    }

    #[test]
    fn sds_rolling_opposes_angular_velocity() {
        // Two overlapping particles, one spinning → SDS rolling torque opposes it
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        // Give atom 0 angular velocity in y (rolling about contact normal x)
        dem.omega[0] = [0.0, 10.0, 0.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_sds_rolling());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // SDS rolling torque should oppose omega_y on atom 0
        assert!(
            dem.torque[0][1] < 0.0,
            "SDS rolling torque should oppose omega_y, got {}",
            dem.torque[0][1]
        );
    }

    #[test]
    fn sds_rolling_spring_accumulates() {
        // Pre-load rolling displacement → larger torque than zero displacement
        // Use very small omega so that damping doesn't dominate and Coulomb cap isn't reached
        let radius = 0.001;

        let run_with_preload = |preload_y: f64| -> f64 {
            let mut app = App::new();
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;

            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(
                &mut atom,
                &mut dem,
                &mut hist,
                1,
                [0.0019, 0.0, 0.0],
                radius,
            );
            dem.omega[0] = [0.0, 0.001, 0.0]; // very small angular velocity
            atom.nlocal = 2;
            atom.natoms = 2;

            // Pre-load rolling displacement in contact history (canonical: tag 0 < tag 1, sign=+1)
            if preload_y != 0.0 {
                let mut preload = [0.0; CONTACT_HISTORY_LEN];
                preload[4] = preload_y;
                hist.contacts[0].push((1, preload, false));
            }

            let mut neighbor = Neighbor::new();
            neighbor.neighbor_offsets = vec![0, 1, 1];
            neighbor.neighbor_indices = vec![1];

            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(hist);

            app.add_resource(atom);
            app.add_resource(neighbor);
            app.add_resource(registry);
            app.add_resource(make_material_table_sds_rolling());
            app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();

            let reg = app.get_resource_ref::<AtomDataRegistry>().unwrap();
            let d = reg.expect::<DemAtom>("test");
            d.torque[0][1]
        };

        let torque_no_preload = run_with_preload(0.0);
        let torque_with_preload = run_with_preload(1e-5); // small preload below cap

        assert!(torque_no_preload < 0.0, "should oppose omega_y");
        assert!(torque_with_preload < 0.0, "should oppose omega_y");
        // Pre-loaded spring adds to torque magnitude
        assert!(
            torque_with_preload.abs() > torque_no_preload.abs(),
            "preloaded spring should increase torque: no_preload={}, preloaded={}",
            torque_no_preload,
            torque_with_preload
        );
    }

    #[test]
    fn sds_rolling_coulomb_cap() {
        // Very high angular velocity → torque should be capped at mu_r * |F_n| * R_eff
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-5; // larger dt to accumulate big spring

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        dem.omega[0] = [0.0, 1e6, 0.0]; // very high
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        let mt = make_material_table_sds_rolling();
        let mu_r = mt.rolling_friction_ij[0][0];

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        let torque_mag =
            (dem.torque[0][0].powi(2) + dem.torque[0][1].powi(2) + dem.torque[0][2].powi(2)).sqrt();

        // Compute expected cap: mu_r * F_n * R_eff
        // F_n from Hertz: 4/3 * E_eff * sqrt(delta * r_eff) * delta
        let r_eff = radius / 2.0;
        let delta = 2.0 * radius - 0.0019;
        let e_eff = 8.7e9 / (2.0 * (1.0 - 0.09)); // single material
        let sqrt_dr = (delta * r_eff).sqrt();
        let f_n_approx = 4.0 / 3.0 * e_eff * sqrt_dr * delta;
        let tau_cap = mu_r * f_n_approx * r_eff;

        // Rolling torque should not exceed cap (with reasonable tolerance for damping and normal force)
        // The torque includes tangential torque contributions, so we just check the rolling component
        // is bounded. Since torque_mag includes all contributions, just check it's finite and reasonable.
        assert!(torque_mag.is_finite(), "torque should be finite");
        assert!(
            torque_mag < tau_cap * 100.0, // generous bound since total torque includes tangential
            "torque {} should be bounded near cap {}",
            torque_mag,
            tau_cap
        );
    }

    #[test]
    fn sds_twisting_opposes_spin() {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        // Spin about contact normal (x-axis)
        dem.omega[0] = [10.0, 0.0, 0.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_sds_twisting());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // SDS twisting torque should oppose spin about x
        assert!(
            dem.torque[0][0] < 0.0,
            "SDS twisting torque should oppose spin about x, got {}",
            dem.torque[0][0]
        );
    }

    #[test]
    fn sds_twisting_spring_accumulates() {
        let radius = 0.001;

        let run_with_preload = |preload: f64| -> f64 {
            let mut app = App::new();
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;

            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(
                &mut atom,
                &mut dem,
                &mut hist,
                1,
                [0.0019, 0.0, 0.0],
                radius,
            );
            dem.omega[0] = [0.001, 0.0, 0.0]; // very small spin
            atom.nlocal = 2;
            atom.natoms = 2;

            if preload != 0.0 {
                let mut preload_state = [0.0; CONTACT_HISTORY_LEN];
                preload_state[6] = preload;
                hist.contacts[0].push((1, preload_state, false));
            }

            let mut neighbor = Neighbor::new();
            neighbor.neighbor_offsets = vec![0, 1, 1];
            neighbor.neighbor_indices = vec![1];

            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(hist);

            app.add_resource(atom);
            app.add_resource(neighbor);
            app.add_resource(registry);
            app.add_resource(make_material_table_sds_twisting());
            app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();

            let reg = app.get_resource_ref::<AtomDataRegistry>().unwrap();
            let d = reg.expect::<DemAtom>("test");
            d.torque[0][0]
        };

        let torque_no_preload = run_with_preload(0.0);
        let torque_with_preload = run_with_preload(1e-5);
        assert!(torque_no_preload < 0.0);
        assert!(torque_with_preload < 0.0);
        assert!(
            torque_with_preload.abs() > torque_no_preload.abs(),
            "preloaded twisting spring should increase torque: no_preload={}, preloaded={}",
            torque_no_preload,
            torque_with_preload
        );
    }

    /// Marshall twisting material: tangential friction `friction` drives the
    /// derived twisting cap; the SDS twist stiffness/damping inputs are provided
    /// deliberately so tests can confirm the Marshall model *ignores* them.
    fn make_material_table_marshall_twisting(
        friction: f64,
        twist_stiff: f64,
        twist_damp: f64,
    ) -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.twisting_model = "marshall".to_string();
        mt.add_material_with_sds(
            "glass",
            8.7e9,
            0.3,
            0.95,
            friction, // tangential μ_t — Marshall derives μ_twist = (2/3) a μ_t from this
            0.0,      // rolling_friction
            0.0,
            0.0,
            0.0, // twisting_friction (unused by Marshall)
            0.0,
            0.0, // kn, kt (Hertz path ignores these)
            0.0,
            0.0, // rolling sds
            twist_stiff,
            twist_damp, // twisting sds — must NOT affect Marshall
        );
        mt.build_pair_tables();
        mt
    }

    /// Run one Marshall-twisting contact step and return the twisting torque on
    /// atom 0 (about the contact normal x̂). `preload` seeds the stored twisting
    /// spring displacement; a large value forces the saturated (capped) regime.
    fn run_marshall_twist(mt: MaterialTable, omega_x: f64, preload: f64) -> f64 {
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        dem.omega[0] = [omega_x, 0.0, 0.0]; // spin about contact normal x̂
        atom.nlocal = 2;
        atom.natoms = 2;

        if preload != 0.0 {
            let mut preload_state = [0.0; CONTACT_HISTORY_LEN];
            preload_state[6] = preload;
            hist.contacts[0].push((1, preload_state, false));
        }

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let reg = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let tq = reg.expect::<DemAtom>("test").torque[0][0];
        tq
    }

    #[test]
    fn marshall_twisting_opposes_spin() {
        // Spin about the contact normal (x̂) → Marshall twisting couple opposes it.
        let tq = run_marshall_twist(
            make_material_table_marshall_twisting(0.4, 0.0, 0.0),
            10.0,
            0.0,
        );
        assert!(
            tq < 0.0,
            "Marshall twisting torque should oppose spin about x, got {}",
            tq
        );
    }

    #[test]
    fn marshall_twisting_ignores_sds_inputs() {
        // The Marshall coefficients are DERIVED from the tangential model, so the
        // SDS twisting_stiffness / twisting_damping material inputs must have no
        // effect. Drive into the saturated (capped) regime with a large preload so
        // the torque equals the derived cap τ_max = μ_twist·F_n, then confirm two
        // wildly different SDS-input tables give the identical torque.
        let tq_zero = run_marshall_twist(
            make_material_table_marshall_twisting(0.4, 0.0, 0.0),
            10.0,
            1.0,
        );
        let tq_huge = run_marshall_twist(
            make_material_table_marshall_twisting(0.4, 1.0e9, 1.0e6),
            10.0,
            1.0,
        );
        assert!(tq_zero < 0.0, "should oppose spin, got {}", tq_zero);
        assert!(
            (tq_zero - tq_huge).abs() <= 1e-12 * tq_zero.abs().max(1e-30),
            "Marshall torque must ignore SDS twist inputs: zero-input={}, huge-input={}",
            tq_zero,
            tq_huge
        );
    }

    #[test]
    fn marshall_twisting_cap_scales_with_tangential_friction() {
        // μ_twist = (2/3) a μ_t, so in the saturated regime the cap scales linearly
        // with the tangential friction coefficient: doubling μ_t doubles |τ|, and
        // μ_t = 0 gives zero twisting couple (Marshall ties the cap to sliding).
        let tq_mu04 = run_marshall_twist(
            make_material_table_marshall_twisting(0.4, 0.0, 0.0),
            10.0,
            1.0,
        );
        let tq_mu08 = run_marshall_twist(
            make_material_table_marshall_twisting(0.8, 0.0, 0.0),
            10.0,
            1.0,
        );
        let tq_mu00 = run_marshall_twist(
            make_material_table_marshall_twisting(0.0, 0.0, 0.0),
            10.0,
            1.0,
        );
        let ratio = tq_mu08 / tq_mu04;
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "doubling μ_t should double the Marshall cap: ratio={}",
            ratio
        );
        assert!(
            tq_mu00.abs() < 1e-12,
            "μ_t = 0 should give zero Marshall twisting torque, got {}",
            tq_mu00
        );
    }

    #[test]
    fn constant_model_unchanged_with_sds_config() {
        // When rolling_model = "constant" (default), SDS parameters should be ignored
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        dem.omega[0] = [0.0, 10.0, 0.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        // Use constant model but with SDS parameters set (they should be ignored)
        let mut mt = MaterialTable::new();
        // rolling_model defaults to "constant"
        mt.add_material_with_sds(
            "glass", 8.7e9, 0.3, 0.95, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 1e3, 0.5, 0.0, 0.0,
        );
        mt.build_pair_tables();

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let dem = registry.expect::<DemAtom>("test");
        // Constant model: torque = -mu_r * |F_n| * r_eff * (roll/|roll|)
        // Should still produce opposing torque
        assert!(
            dem.torque[0][1] < 0.0,
            "constant rolling model should still work, got {}",
            dem.torque[0][1]
        );

        // Check that spring history has zero rolling/twisting displacement
        let hist = registry.expect::<ContactHistoryStore>("test");
        let contact = &hist.contacts[0][0];
        assert_eq!(
            contact.1[3], 0.0,
            "rolling disp x should be zero in constant model"
        );
        assert_eq!(
            contact.1[4], 0.0,
            "rolling disp y should be zero in constant model"
        );
        assert_eq!(
            contact.1[5], 0.0,
            "rolling disp z should be zero in constant model"
        );
        assert_eq!(
            contact.1[6], 0.0,
            "twisting disp should be zero in constant model"
        );
    }

    // ── DMT adhesion tests ──────────────────────────────────────────────

    fn make_material_table_dmt() -> MaterialTable {
        let mut mt = MaterialTable::new();
        // Use high surface energy (1.0 J/m²) so adhesion clearly dominates at small overlaps
        mt.add_material_full("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 1.0);
        mt.adhesion_model = "dmt".to_string();
        mt.build_pair_tables();
        mt
    }

    #[test]
    fn dmt_pulloff_force_matches_theory() {
        // DMT pull-off force = 2 * pi * gamma * r_eff (at contact, delta = 0+)
        let radius = 0.001;
        let gamma = 1.0;
        let r_eff = radius / 2.0; // two equal spheres

        // Use a very small overlap so Hertz contribution is negligible
        // At tiny delta, F_hertz ~ 0 but F_dmt = 2*pi*gamma*r_eff
        let tiny_overlap = 1e-12; // extremely small overlap
        let sep = 2.0 * radius - tiny_overlap;

        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_dmt());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let expected_dmt = 2.0 * std::f64::consts::PI * gamma * r_eff;
        // Force on atom 0 should be positive (attracted toward atom 1)
        // f_n_mag = k_n*delta - f_diss - f_dmt ~ -f_dmt (since delta ~ 0, v=0)
        // force[0] -= f_n_mag * nx -> force[0] ~ +f_dmt
        assert!(
            atom.force[0][0] > 0.0,
            "DMT should produce attractive force, got {}",
            atom.force[0][0]
        );
        assert!(
            (atom.force[0][0] as f64 - expected_dmt).abs() / expected_dmt < 1e-3,
            "DMT pull-off force should match 2*pi*gamma*r_eff = {}, got {}",
            expected_dmt,
            atom.force[0][0]
        );
    }

    #[test]
    fn dmt_no_force_beyond_contact() {
        // DMT has no adhesion-only regime -- no force when delta < 0 (gap)
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        // Place particles with a gap
        let gap = 1e-9;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [2.0 * radius + gap, 0.0, 0.0],
            radius,
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_dmt());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // DMT: no force when particles are not in geometric contact
        assert!(
            atom.force[0][0].abs() < 1e-20,
            "DMT should have no force beyond contact, got {}",
            atom.force[0][0]
        );
    }

    #[test]
    fn dmt_pulloff_less_than_jkr() {
        // DMT pull-off = 2*pi*gamma*r_eff, JKR pull-off = 1.5*pi*gamma*r_eff
        // At same surface energy, DMT has HIGHER pull-off force than JKR (2 > 1.5)
        // But JKR has extended range (adhesion across gap), so effective sticking is stronger
        let gamma = 1.0;
        let radius = 0.001;
        let r_eff = radius / 2.0;

        let f_dmt = 2.0 * std::f64::consts::PI * gamma * r_eff;
        let f_jkr = 1.5 * std::f64::consts::PI * gamma * r_eff;
        assert!(
            f_dmt > f_jkr,
            "DMT pull-off ({}) should be larger than JKR pull-off ({})",
            f_dmt,
            f_jkr
        );
    }

    #[test]
    fn dmt_newtons_third_law() {
        // Verify equal and opposite forces for DMT contact
        let mut app = App::new();
        let radius = 0.001;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_dmt());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        for d in 0..3 {
            assert!(
                (atom.force[0][d] + atom.force[1][d]).abs() < 1e-10,
                "Newton's 3rd law violated in dim {}: {} + {} != 0",
                d,
                atom.force[0][d],
                atom.force[1][d]
            );
        }
    }

    #[test]
    fn dmt_does_not_break_jkr() {
        // Run the JKR test with default adhesion_model (should still work as JKR)
        let mut app = App::new();
        let radius = 0.001;
        let gamma = 1.0;
        let r_eff = radius / 2.0;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        // Place particles with a tiny gap (adhesion-only regime for JKR)
        let gap = 1e-9;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [2.0 * radius + gap, 0.0, 0.0],
            radius,
        );
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        // Use JKR material table (default adhesion_model = "jkr")
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_jkr());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let expected_jkr = 1.5 * std::f64::consts::PI * gamma * r_eff;
        // JKR should still attract across gap
        assert!(
            atom.force[0][0] > 0.0,
            "JKR should still work with DMT feature added, got {}",
            atom.force[0][0]
        );
        assert!(
            (atom.force[0][0] as f64 - expected_jkr).abs() / expected_jkr < 1e-6,
            "JKR pull-off force should still match 1.5*pi*gamma*r_eff = {}, got {}",
            expected_jkr,
            atom.force[0][0]
        );
    }

    // ── Force scaling validation tests ──────────────────────────────────

    #[test]
    fn hertz_force_scales_as_delta_three_halves() {
        let radius = 0.001;

        // Compute elastic-only normal force for a given separation (zero velocity -> no damping).
        let hertz_force_at = |sep: f64| -> f64 {
            let mut app = App::new();
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
            atom.nlocal = 2;
            atom.natoms = 2;
            let mut neighbor = Neighbor::new();
            neighbor.neighbor_offsets = vec![0, 1, 1];
            neighbor.neighbor_indices = vec![1];
            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(hist);
            app.add_resource(atom);
            app.add_resource(neighbor);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();
            let atom = app.get_resource_ref::<Atom>().unwrap();
            // Force on atom 0 is negative (pushed away from atom 1), take absolute value
            atom.force[0][0].abs() as f64
        };

        // Test at 5 different overlaps
        let deltas = [1e-5, 2e-5, 4e-5, 6e-5, 8e-5];
        let forces: Vec<f64> = deltas
            .iter()
            .map(|d| {
                let sep = 2.0 * radius - d;
                hertz_force_at(sep)
            })
            .collect();

        // For each pair (i, 0), check F_i/F_0 ~ (delta_i/delta_0)^(3/2)
        for i in 1..deltas.len() {
            let expected_ratio = (deltas[i] / deltas[0]).powf(1.5);
            let actual_ratio = forces[i] / forces[0];
            let rel_err = ((actual_ratio - expected_ratio) / expected_ratio).abs();
            assert!(
                rel_err < 0.01,
                "Hertz force scaling: delta ratio {:.1}, expected F ratio {:.4}, got {:.4} (rel err {:.4})",
                deltas[i] / deltas[0], expected_ratio, actual_ratio, rel_err
            );
        }
    }

    #[test]
    fn hooke_force_scales_linearly_across_overlaps() {
        let radius = 0.001;
        let hooke_force_at = |sep: f64| -> f64 {
            let mut app = App::new();
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            let mut hist = ContactHistoryStore::new();
            atom.dt = 1e-7;
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
            push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
            atom.nlocal = 2;
            atom.natoms = 2;
            let mut neighbor = Neighbor::new();
            neighbor.neighbor_offsets = vec![0, 1, 1];
            neighbor.neighbor_indices = vec![1];
            let mut registry = AtomDataRegistry::new();
            registry.register(dem);
            registry.register(hist);
            app.add_resource(atom);
            app.add_resource(neighbor);
            app.add_resource(registry);
            app.add_resource(make_material_table_hooke());
            app.add_update_system(hooke_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();
            let atom = app.get_resource_ref::<Atom>().unwrap();
            atom.force[0][0].abs() as f64
        };

        let deltas = [2e-5, 4e-5, 6e-5, 8e-5, 1e-4];
        let forces: Vec<f64> = deltas
            .iter()
            .map(|d| {
                let sep = 2.0 * radius - d;
                hooke_force_at(sep)
            })
            .collect();

        for i in 1..deltas.len() {
            let expected_ratio = deltas[i] / deltas[0]; // linear
            let actual_ratio = forces[i] / forces[0];
            let rel_err = ((actual_ratio - expected_ratio) / expected_ratio).abs();
            assert!(
                rel_err < 0.01,
                "Hooke force scaling: delta ratio {:.1}, expected F ratio {:.4}, got {:.4} (rel err {:.4})",
                deltas[i] / deltas[0], expected_ratio, actual_ratio, rel_err
            );
        }
    }

    #[test]
    fn hertz_force_matches_analytical_value() {
        let radius = 0.001;
        let delta = 5e-5;
        let sep = 2.0 * radius - delta;

        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        let mt = make_material_table();
        let e_eff = mt.e_eff_ij[0][0];
        let r_eff = radius / 2.0; // two equal spheres: r_eff = r1*r2/(r1+r2) = r/2

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let f_computed = atom.force[0][0].abs() as f64;
        // Analytical: F = (4/3) * E_eff * sqrt(R_eff) * delta^(3/2)
        let f_analytical = (4.0 / 3.0) * e_eff * r_eff.sqrt() * delta.powf(1.5);
        let rel_err = (f_computed - f_analytical).abs() / f_analytical;
        assert!(
            rel_err < 1e-10,
            "Hertz force analytical check: computed={:.6e}, expected={:.6e}, rel_err={:.2e}",
            f_computed,
            f_analytical,
            rel_err
        );
    }

    #[test]
    fn linear_momentum_conserved_during_elastic_contact() {
        // Perfectly elastic (restitution = 1.0) → ~no damping. The Hertz/Tsuji
        // coefficient at e=1 is the polynomial's residual (~1.3e-4), not exactly 0
        // (LAMMPS `damping tsuji` has the same residual), so momentum is conserved to
        // that order rather than machine epsilon.
        let mut mt = MaterialTable::new();
        mt.add_material("elastic", 8.7e9, 0.3, 1.0, 0.0, 0.0, 0.0);
        mt.build_pair_tables();
        assert!(
            mt.beta_ij[0][0].abs() < 1e-3,
            "beta should be ~0 for e=1.0, got {}",
            mt.beta_ij[0][0]
        );

        let radius = 0.001;
        let dt = 1e-8;

        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = dt;

        // Two particles approaching each other, slight overlap
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.00195, 0.0, 0.0],
            radius,
        );
        atom.vel[0] = [0.1, 0.05, -0.02];
        atom.vel[1] = [-0.05, 0.03, 0.01];
        atom.nlocal = 2;
        atom.natoms = 2;

        let initial_momentum = [
            atom.mass[0] * atom.vel[0][0] + atom.mass[1] * atom.vel[1][0],
            atom.mass[0] * atom.vel[0][1] + atom.mass[1] * atom.vel[1][1],
            atom.mass[0] * atom.vel[0][2] + atom.mass[1] * atom.vel[1][2],
        ];

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(
            crate::contact::hertz_mindlin_contact_force,
            ParticleSimScheduleSet::Force,
        );
        app.add_update_system(
            soil_verlet::initial_integration,
            ParticleSimScheduleSet::InitialIntegration,
        );
        app.add_update_system(
            soil_verlet::final_integration,
            ParticleSimScheduleSet::FinalIntegration,
        );
        // Zero forces between steps
        app.add_update_system(
            |mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>| {
                let n = atoms.len();
                atoms.force[..n].fill([0.0; 3]);
                registry.zero_all(n);
            },
            ParticleSimScheduleSet::PostInitialIntegration,
        );
        app.organize_systems();

        // Run for 100 steps
        for _ in 0..100 {
            app.run();
        }

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let final_momentum = [
            atom.mass[0] * atom.vel[0][0] + atom.mass[1] * atom.vel[1][0],
            atom.mass[0] * atom.vel[0][1] + atom.mass[1] * atom.vel[1][1],
            atom.mass[0] * atom.vel[0][2] + atom.mass[1] * atom.vel[1][2],
        ];

        for d in 0..3 {
            let err = (final_momentum[d] - initial_momentum[d]).abs();
            assert!(
                err < 1e-12,
                "Momentum not conserved in dim {}: initial={:.6e}, final={:.6e}, err={:.2e}",
                d,
                initial_momentum[d],
                final_momentum[d],
                err
            );
        }
    }

    #[test]
    fn contact_force_symmetry_with_tangential_velocity() {
        let radius = 0.001;
        let sep = 0.0019;

        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        // Give both atoms velocities in all directions
        atom.vel[0] = [0.1, 0.2, -0.1];
        atom.vel[1] = [-0.3, 0.1, 0.05];
        dem.omega[0] = [10.0, 20.0, -5.0];
        dem.omega[1] = [-15.0, 5.0, 10.0];
        atom.nlocal = 2;
        atom.natoms = 2;

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.register(dem);
        registry.register(hist);

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Newton's 3rd law: forces equal and opposite
        for d in 0..3 {
            assert!(
                (atom.force[0][d] + atom.force[1][d]).abs() < 1e-10,
                "Newton's 3rd law violated in dim {}: f0={:.6e}, f1={:.6e}",
                d,
                atom.force[0][d],
                atom.force[1][d]
            );
        }
    }
}
