//! The bond force kernel and its per-rebuild pair cache.
//!
//! [`bond_force`] walks a flat list of `(owner, partner)` pairs that is built
//! once per neighbor rebuild ([`BondPairCache`]), instead of rebuilding a
//! tag → index hash map and searching each atom's history list every step.
//! Local and ghost indices only change when SOIL migrates, sorts, or rebuilds
//! borders, which all happen on a rebuild step, so the cache is rebuilt then,
//! after any bond breaks, and whenever the atom counts change.
//!
//! Two force laws are available (see [`BondFormulation`]):
//!
//! * [`BondFormulation::Dirt`] keeps DIRT's historical law.
//! * [`BondFormulation::BpmGuo`] follows LAMMPS `bond_style bpm/guo`
//!   (`bond_bpm_guo.cpp`, EMS2026 revision) step for step.

use std::f64::consts::PI;

use grass_scheduler::prelude::*;
use soil_core::{
    Atom, BondEntry, BondStore, Domain, Neighbor, Optional, ParticlesWith, VirialStress, Write,
};

use crate::breakage;
use crate::plasticity;
use crate::{
    BondBreakage, BondConfig, BondFormulation, BondHistoryEntry, BondHistoryStore, BondMetrics,
    BondPlasticity, BondRadiusRule,
};
use dirt_atom::DemAtom;

const NONE: u32 = u32::MAX;

/// One bond as the kernel sees it: owner `i` (local, lower tag), partner `j`
/// (local or ghost), the bond type, the index of its history in
/// `BondHistoryStore::history[i]`, plus the per-bond constants.
#[derive(Clone, Copy, Debug)]
struct BondPair {
    i: u32,
    j: u32,
    h: u32,
    bond_type: u32,
    r0: f64,
    r_b: f64,
    area: f64,
    jpol: f64,
    iben: f64,
    k_n: f64,
    k_t: f64,
    k_tor: f64,
    k_bend: f64,
    gamma_n: f64,
    gamma_t: f64,
    gamma_tor: f64,
    gamma_bend: f64,
}

/// Flat per-rebuild list of the bonds this rank computes.
#[derive(Default)]
pub struct BondPairCache {
    pairs: Vec<BondPair>,
    /// Dense tag → local/ghost index map (`u32::MAX` = absent). Local atoms win
    /// over ghost images of the same tag.
    slot: Vec<u32>,
    valid: bool,
    nlocal: usize,
    ntotal: usize,
    /// Bonds whose partner was not visible when the cache was built.
    missing: usize,
    /// Bond-force calls so far (drives `breakage_start_step`).
    steps: u64,
    /// Half-bonds dropped because the owning rank had already broken them.
    pub pruned_orphans: u64,
}

impl BondPairCache {
    /// Force a rebuild at the next force call (e.g. after bonds were added).
    pub fn invalidate(&mut self) {
        self.valid = false;
    }

    #[inline]
    fn lookup(&self, tag: u32) -> Option<usize> {
        match self.slot.get(tag as usize) {
            Some(&s) if s != NONE => Some(s as usize),
            _ => None,
        }
    }

    fn build_slot_map(&mut self, atoms: &Atom) {
        let n = atoms.len();
        let max_tag = atoms.tag[..n].iter().copied().max().unwrap_or(0) as usize;
        self.slot.clear();
        self.slot.resize(max_tag + 1, NONE);
        for idx in 0..n {
            let s = &mut self.slot[atoms.tag[idx] as usize];
            if *s == NONE {
                *s = idx as u32;
            }
        }
    }
}

/// Minimum-image separation `pos[j] − pos[i]` (triclinic-aware).
#[inline]
fn min_image(atoms: &Atom, domain: &Domain, i: usize, j: usize) -> [f64; 3] {
    let mut d = [
        atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64,
        atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64,
        atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64,
    ];
    let pflags = domain.periodic_flags();
    if !(pflags[0] || pflags[1] || pflags[2]) {
        return d;
    }
    let box_size = domain.size;
    if domain.triclinic {
        // Map to fractional (lamda) coordinates through H⁻¹ (H upper-triangular
        // with tilt factors [xy, xz, yz]), wrap periodic components, map back.
        let [xy, xz, yz] = domain.tilt;
        let lz = if box_size[2] != 0.0 {
            d[2] / box_size[2]
        } else {
            0.0
        };
        let ly = if box_size[1] != 0.0 {
            (d[1] - yz * lz) / box_size[1]
        } else {
            0.0
        };
        let lx = if box_size[0] != 0.0 {
            (d[0] - xy * ly - xz * lz) / box_size[0]
        } else {
            0.0
        };
        let mut lam = [lx, ly, lz];
        for k in 0..3 {
            if pflags[k] {
                lam[k] -= lam[k].round();
            }
        }
        d[0] = box_size[0] * lam[0] + xy * lam[1] + xz * lam[2];
        d[1] = box_size[1] * lam[1] + yz * lam[2];
        d[2] = box_size[2] * lam[2];
    } else {
        for k in 0..3 {
            if pflags[k] {
                d[k] -= box_size[k] * (d[k] / box_size[k]).round();
            }
        }
    }
    d
}

/// Per-bond failure thresholds, drawn exactly as `init_bond_history` does.
pub(crate) fn sample_thresholds(
    breakage: &BondBreakage,
    bond_type: u32,
    tag_a: u32,
    tag_b: u32,
    r0: f64,
) -> [f64; 4] {
    let crit = breakage.resolve(bond_type);
    let u = if crit.uses_guo_hash() {
        [breakage::guo_weibull_u(tag_a, tag_b, breakage.seed); 4]
    } else {
        breakage::per_bond_uniform_samples(tag_a, tag_b, breakage.seed)
    };
    crit.sample(r0, u).t
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// LAMMPS `reproject()`: rotate `v` onto the axis `n` (`along`) or into the
/// plane normal to it, keeping its magnitude. A zero projection leaves `v`.
#[inline]
fn reproject(v: &mut [f64; 3], n: [f64; 3], along: bool) {
    let len = norm(*v);
    let p = dot(*v, n);
    let u = if along {
        [p * n[0], p * n[1], p * n[2]]
    } else {
        [v[0] - p * n[0], v[1] - p * n[1], v[2] - p * n[2]]
    };
    let ulen = norm(u);
    if ulen > 0.0 {
        let s = len / ulen;
        *v = [s * u[0], s * u[1], s * u[2]];
    }
}

/// Rebuild the pair cache: map tags to indices, resolve deferred `r0` for
/// bonds that crossed a rank cut at load time, seed any missing history,
/// drop half-bonds that another rank already broke, and precompute each
/// owned bond's stiffness and damping constants.
#[allow(clippy::too_many_arguments)]
fn rebuild_cache(
    cache: &mut BondPairCache,
    atoms: &Atom,
    dem: &DemAtom,
    bonds: &mut BondStore,
    hist: &mut BondHistoryStore,
    cfg: &BondConfig,
    breakage: &BondBreakage,
    domain: &Domain,
) {
    cache.build_slot_map(atoms);
    let nlocal = atoms.nlocal() as usize;
    cache.nlocal = nlocal;
    cache.ntotal = atoms.len();
    cache.pairs.clear();
    cache.missing = 0;
    while hist.history.len() < bonds.bonds.len() {
        hist.history.push(Vec::new());
    }

    let guo = cfg.formulation == BondFormulation::BpmGuo;
    let ratio = cfg.bond_radius_ratio;
    let radius_rule = cfg.bond_radius.unwrap_or(if guo {
        BondRadiusRule::Mean
    } else {
        BondRadiusRule::Min
    });
    let e_mod = cfg.youngs_modulus;
    let g_mod = cfg.shear_modulus;
    let e_bend_mod = cfg.bending_modulus.or(e_mod);
    let g_tor_mod = cfg.twist_modulus.or(g_mod);

    for i in 0..nlocal.min(bonds.bonds.len()) {
        let tag_i = atoms.tag[i];

        // A bond the lower-tag owner already broke on its rank survives here
        // only as this atom's half. The owner's row (local, or a ghost copy
        // refreshed at the last border exchange) no longer lists us: drop it.
        let before = bonds.bonds[i].len();
        {
            let (rows_before, rows_rest) = bonds.bonds.split_at_mut(i);
            let (row_i, rows_after) = rows_rest.split_first_mut().expect("row i exists");
            row_i.retain(|b: &BondEntry| {
                if b.partner_tag >= tag_i {
                    return true;
                }
                let Some(p) = cache.lookup(b.partner_tag) else {
                    return true;
                };
                let row_p = if p < i {
                    &rows_before[p]
                } else if p > i {
                    match rows_after.get(p - i - 1) {
                        Some(r) => r,
                        None => return true,
                    }
                } else {
                    return true;
                };
                row_p.iter().any(|e| e.partner_tag == tag_i)
            });
        }
        if bonds.bonds[i].len() != before {
            cache.pruned_orphans += (before - bonds.bonds[i].len()) as u64;
            let keep: Vec<u32> = bonds.bonds[i].iter().map(|b| b.partner_tag).collect();
            hist.history[i].retain(|h| keep.contains(&h.partner_tag));
        }

        for b_idx in 0..bonds.bonds[i].len() {
            let partner = bonds.bonds[i][b_idx].partner_tag;
            let Some(j) = cache.lookup(partner) else {
                cache.missing += 1;
                continue;
            };
            // Resolve a deferred reference length now that the partner is
            // visible (both halves of a cross-rank bond see the same positions,
            // so they agree on r0).
            if bonds.bonds[i][b_idx].r0 <= 0.0 {
                let d = min_image(atoms, domain, i, j);
                let r0 = norm(d);
                if r0 <= 0.0 {
                    continue;
                }
                bonds.bonds[i][b_idx].r0 = r0;
            }
            if tag_i >= partner {
                continue; // the lower tag computes the bond
            }
            let bond = bonds.bonds[i][b_idx].clone();
            let h_idx = match hist.history[i]
                .iter()
                .position(|h| h.partner_tag == partner)
            {
                Some(h) => h,
                None => {
                    hist.history[i].push(BondHistoryEntry {
                        partner_tag: partner,
                        thresholds: sample_thresholds(
                            breakage,
                            bond.bond_type,
                            tag_i,
                            partner,
                            bond.r0,
                        ),
                        ..Default::default()
                    });
                    hist.history[i].len() - 1
                }
            };

            let r_b = ratio
                * match radius_rule {
                    BondRadiusRule::Min => dem.radius[i].min(dem.radius[j]),
                    BondRadiusRule::Mean => 0.5 * (dem.radius[i] + dem.radius[j]),
                };
            let area = PI * r_b * r_b;
            let jpol = 0.5 * PI * r_b.powi(4);
            let iben = 0.5 * jpol;
            let len = bond.r0;
            let k_n = e_mod.map_or(cfg.normal_stiffness, |e| e * area / len);
            let k_t = g_mod.map_or(cfg.shear_stiffness, |g| g * area / len);
            let k_tor = g_tor_mod.map_or(cfg.twist_stiffness, |g| g * jpol / len);
            let k_bend = e_bend_mod.map_or(cfg.bending_stiffness, |e| e * iben / len);

            // Damping masses. DIRT: reduced mass and reduced moment of inertia of
            // the pair. bpm/guo: half the owner's mass, and half the owner's
            // sphere inertia taken about the bond radius (bond_inertia()).
            let (m_damp, moi_damp) = if guo {
                let m_i = atoms.mass[i] as f64;
                (0.5 * m_i, 0.5 * 0.4 * m_i * r_b * r_b)
            } else {
                let m_i = atoms.mass[i] as f64;
                let m_j = atoms.mass[j] as f64;
                let m_red = if m_i + m_j > 0.0 {
                    m_i * m_j / (m_i + m_j)
                } else {
                    0.0
                };
                let inv = |x: f64| if x > 0.0 { 1.0 / x } else { 0.0 };
                let moi_i = inv(dem.inv_inertia[i]);
                let moi_j = inv(dem.inv_inertia[j]);
                let moi_red = if moi_i + moi_j > 0.0 {
                    moi_i * moi_j / (moi_i + moi_j)
                } else {
                    0.0
                };
                (m_red, moi_red)
            };
            let crit = |beta: f64, m: f64, k: f64| 2.0 * beta * (m * k.max(0.0)).sqrt();
            let gamma_n = cfg
                .normal_damping
                .unwrap_or_else(|| crit(cfg.beta_normal, m_damp, k_n));
            let gamma_t = cfg
                .shear_damping
                .unwrap_or_else(|| crit(cfg.beta_shear, m_damp, k_t));
            let gamma_tor = cfg
                .twist_damping
                .unwrap_or_else(|| crit(cfg.beta_twist, moi_damp, k_tor));
            let gamma_bend = cfg
                .bending_damping
                .unwrap_or_else(|| crit(cfg.beta_bending, moi_damp, k_bend));

            cache.pairs.push(BondPair {
                i: i as u32,
                j: j as u32,
                bond_type: bond.bond_type,
                h: h_idx as u32,
                r0: bond.r0,
                r_b,
                area,
                jpol,
                iben,
                k_n,
                k_t,
                k_tor,
                k_bend,
                gamma_n,
                gamma_t,
                gamma_tor,
                gamma_bend,
            });
        }
    }
    cache.valid = true;
}

/// Loads applied by one bond on its owner `i` (the partner gets the force
/// negated; see the torque rule in each formulation).
struct Applied {
    force: [f64; 3],
    torque_i: [f64; 3],
    torque_j: [f64; 3],
    strain: f64,
}

/// DIRT's historical law, unchanged: shear displacement with its axial part
/// removed each step, a single total relative rotation split into twist and
/// bending about the current axis, lever arm `r₀/2`, and the configured
/// plasticity return maps.
#[allow(clippy::too_many_arguments)]
#[inline]
fn dirt_bond(
    p: &BondPair,
    d: [f64; 3],
    dist: f64,
    atoms: &Atom,
    dem: &DemAtom,
    h: &mut BondHistoryEntry,
    plasticity: &plasticity::BondPlasticityModel,
    crit: &dyn breakage::BreakageCriterion,
    breakage_active: bool,
    dt: f64,
) -> Option<Applied> {
    let (i, j) = (p.i as usize, p.j as usize);
    let nhat = [d[0] / dist, d[1] / dist, d[2] / dist];
    let half_l = 0.5 * p.r0;
    let r1 = [half_l * nhat[0], half_l * nhat[1], half_l * nhat[2]];
    let w_i = dem.omega[i];
    let w_j = dem.omega[j];
    let wr_i = cross(w_i, r1);
    let wr_j = cross(w_j, r1);
    let v_rel = [
        (atoms.vel[j][0] as f64 - wr_j[0]) - (atoms.vel[i][0] as f64 + wr_i[0]),
        (atoms.vel[j][1] as f64 - wr_j[1]) - (atoms.vel[i][1] as f64 + wr_i[1]),
        (atoms.vel[j][2] as f64 - wr_j[2]) - (atoms.vel[i][2] as f64 + wr_i[2]),
    ];
    let v_n_s = dot(v_rel, nhat);
    let v_t = [
        v_rel[0] - v_n_s * nhat[0],
        v_rel[1] - v_n_s * nhat[1],
        v_rel[2] - v_n_s * nhat[2],
    ];
    let delta = dist - p.r0;

    let (f_n_cons, eps_p, eps_max) = match plasticity.axial.as_ref() {
        Some(model) => {
            let eps = if p.r0 > 0.0 { delta / p.r0 } else { 0.0 };
            plasticity::update_axial(eps, h.eps_p_axial, h.eps_max_axial, p.k_n, p.r0, model)
        }
        None => (p.k_n * delta, h.eps_p_axial, h.eps_max_axial),
    };
    h.eps_p_axial = eps_p;
    h.eps_max_axial = eps_max;
    let f_n_mag = f_n_cons + p.gamma_n * v_n_s;

    let s_n = dot(h.delta_t, nhat);
    for k in 0..3 {
        h.delta_t[k] -= s_n * nhat[k];
        h.delta_t[k] += v_t[k] * dt;
    }
    let ds = h.delta_t;
    let f_t = [
        p.k_t * ds[0] + p.gamma_t * v_t[0],
        p.k_t * ds[1] + p.gamma_t * v_t[1],
        p.k_t * ds[2] + p.gamma_t * v_t[2],
    ];

    let w_rel = [w_j[0] - w_i[0], w_j[1] - w_i[1], w_j[2] - w_i[2]];
    let w_rel_n_s = dot(w_rel, nhat);
    let w_n = [
        w_rel_n_s * nhat[0],
        w_rel_n_s * nhat[1],
        w_rel_n_s * nhat[2],
    ];
    let w_t = [w_rel[0] - w_n[0], w_rel[1] - w_n[1], w_rel[2] - w_n[2]];
    for k in 0..3 {
        h.delta_theta[k] += w_rel[k] * dt;
    }
    let dth = h.delta_theta;
    let dth_n_s = dot(dth, nhat);
    let dth_twist = [dth_n_s * nhat[0], dth_n_s * nhat[1], dth_n_s * nhat[2]];
    let dth_bend = [
        dth[0] - dth_twist[0],
        dth[1] - dth_twist[1],
        dth[2] - dth_twist[2],
    ];
    let m_tor = [
        p.k_tor * dth_twist[0] + p.gamma_tor * w_n[0],
        p.k_tor * dth_twist[1] + p.gamma_tor * w_n[1],
        p.k_tor * dth_twist[2] + p.gamma_tor * w_n[2],
    ];
    let (m_bend_cons, theta_p, theta_max) = match plasticity.bending.as_ref() {
        Some(model) => plasticity::update_bending(
            dth_bend,
            h.theta_p_bend,
            h.theta_max_bend,
            p.k_bend,
            model,
            p.r_b,
            p.r0,
        ),
        None => (
            [
                p.k_bend * dth_bend[0],
                p.k_bend * dth_bend[1],
                p.k_bend * dth_bend[2],
            ],
            h.theta_p_bend,
            h.theta_max_bend,
        ),
    };
    h.theta_p_bend = theta_p;
    h.theta_max_bend = theta_max;
    let m_bend = [
        m_bend_cons[0] + p.gamma_bend * w_t[0],
        m_bend_cons[1] + p.gamma_bend * w_t[1],
        m_bend_cons[2] + p.gamma_bend * w_t[2],
    ];

    if breakage_active {
        let l_now = dist.max(f64::MIN_POSITIVE);
        let geom = breakage::BondGeom {
            r_b: p.r_b,
            area: p.area,
            iben: p.iben,
            jpol: p.jpol,
            l0: p.r0,
        };
        let loads = breakage::BondLoads {
            f_n: f_n_mag,
            f_t_mag: norm(f_t),
            m_bend_mag: norm(m_bend),
            m_tor_mag: norm(m_tor),
        };
        let kin = breakage::BondKinematics {
            eps_axial: delta / p.r0,
            gamma_shear: norm(ds) / l_now,
            kappa_bend: norm(dth_bend) / l_now,
            kappa_tor: dth_n_s.abs() / l_now,
        };
        let thr = breakage::BondThresholds { t: h.thresholds };
        if crit.check(&geom, &loads, &kin, &thr).is_some() {
            return None;
        }
    }

    let tau_shear = cross(r1, f_t);
    let m_total = [
        m_tor[0] + m_bend[0],
        m_tor[1] + m_bend[1],
        m_tor[2] + m_bend[2],
    ];
    Some(Applied {
        force: [
            f_n_mag * nhat[0] + f_t[0],
            f_n_mag * nhat[1] + f_t[1],
            f_n_mag * nhat[2] + f_t[2],
        ],
        torque_i: [
            tau_shear[0] + m_total[0],
            tau_shear[1] + m_total[1],
            tau_shear[2] + m_total[2],
        ],
        torque_j: [
            tau_shear[0] - m_total[0],
            tau_shear[1] - m_total[1],
            tau_shear[2] - m_total[2],
        ],
        strain: delta / p.r0,
    })
}

/// LAMMPS `bond_style bpm/guo` (EMS2026), step for step:
///
/// 1. Re-align the stored shear force, twist and bending moments with the
///    current axis, keeping their magnitudes (`reproject()`).
/// 2. Relative velocities at the bond midpoint, lever arm `L/2` with the
///    *current* length.
/// 3. Accumulate the geometric bend angle, shear displacement and twist angle
///    (never reprojected; they feed the strain break criteria and the bend
///    plasticity peak).
/// 4. Plasticity: bend peak, then the tension-only axial rule.
/// 5. Axial force from the total stretch `K_n (L − L₀ − offset)`.
/// 6. Add this step's shear, twist and bend increments (bend scaled by the
///    plastic slope at a new peak angle).
/// 7. Break check on the elastic loads; a broken bond applies nothing.
/// 8. Damping (owner mass) and application.
#[allow(clippy::too_many_arguments)]
#[inline]
fn guo_bond(
    p: &BondPair,
    d: [f64; 3],
    dist: f64,
    atoms: &Atom,
    dem: &DemAtom,
    h: &mut BondHistoryEntry,
    plasticity: &plasticity::BondPlasticityModel,
    crit: &dyn breakage::BreakageCriterion,
    breakage_active: bool,
    dt: f64,
) -> Option<Applied> {
    let (i, j) = (p.i as usize, p.j as usize);
    let n12 = [d[0] / dist, d[1] / dist, d[2] / dist];

    // 1. The shear history is stored as a displacement (force / K_t); since
    //    K_t is constant per bond this reprojects exactly like the force.
    reproject(&mut h.delta_t, n12, false);
    reproject(&mut h.twist_spring, n12, true);
    reproject(&mut h.bend_spring, n12, false);

    // 2.
    let r1 = [
        0.5 * dist * n12[0],
        0.5 * dist * n12[1],
        0.5 * dist * n12[2],
    ];
    let r2 = [-r1[0], -r1[1], -r1[2]];
    let w_i = dem.omega[i];
    let w_j = dem.omega[j];
    let v1w = cross(w_i, r1);
    let v2w = cross(w_j, r2);
    let vcdiff = [
        (atoms.vel[j][0] as f64 + v2w[0]) - (atoms.vel[i][0] as f64 + v1w[0]),
        (atoms.vel[j][1] as f64 + v2w[1]) - (atoms.vel[i][1] as f64 + v1w[1]),
        (atoms.vel[j][2] as f64 + v2w[2]) - (atoms.vel[i][2] as f64 + v1w[2]),
    ];
    let vn12 = dot(vcdiff, n12);
    let vn = [vn12 * n12[0], vn12 * n12[1], vn12 * n12[2]];
    let vt = [vcdiff[0] - vn[0], vcdiff[1] - vn[1], vcdiff[2] - vn[2]];
    let w_rel = [w_j[0] - w_i[0], w_j[1] - w_i[1], w_j[2] - w_i[2]];
    let wn12 = dot(w_rel, n12);
    let wn = [wn12 * n12[0], wn12 * n12[1], wn12 * n12[2]];
    let wt = [w_rel[0] - wn[0], w_rel[1] - wn[1], w_rel[2] - wn[2]];

    // 3. Geometric accumulators (delta_theta holds the bend angle here).
    for k in 0..3 {
        h.delta_theta[k] += wt[k] * dt;
        h.shear_acc[k] += vt[k] * dt;
        h.twist_acc[k] += wn[k] * dt;
    }
    let theta_bond = norm(h.delta_theta);

    // 4. Bend plasticity marks a new peak angle; axial plasticity updates the
    //    tension-only offset.
    let mut bend_scale = 1.0;
    if let Some(model) = plasticity.bending.as_ref() {
        if theta_bond > h.theta_max_bend {
            h.theta_max_bend = theta_bond;
            bend_scale = model.slope_multiplier(theta_bond, p.k_bend, p.r_b, p.r0);
        }
    }
    let ovlp = dist - p.r0;
    let eps = ovlp / p.r0;
    // 5.
    let f_n_cons = match plasticity.axial.as_ref() {
        Some(model) => {
            let (f, eps_p, eps_max) =
                plasticity::update_axial(eps, h.eps_p_axial, h.eps_max_axial, p.k_n, p.r0, model);
            h.eps_p_axial = eps_p;
            h.eps_max_axial = eps_max;
            f
        }
        None => p.k_n * ovlp,
    };

    // 6.
    for k in 0..3 {
        h.delta_t[k] += vt[k] * dt;
        h.twist_spring[k] += wn[k] * dt;
        h.bend_spring[k] += bend_scale * wt[k] * dt;
    }
    let ftan = [
        p.k_t * h.delta_t[0],
        p.k_t * h.delta_t[1],
        p.k_t * h.delta_t[2],
    ];
    let twist = [
        p.k_tor * h.twist_spring[0],
        p.k_tor * h.twist_spring[1],
        p.k_tor * h.twist_spring[2],
    ];
    let bend = [
        p.k_bend * h.bend_spring[0],
        p.k_bend * h.bend_spring[1],
        p.k_bend * h.bend_spring[2],
    ];

    // 7.
    if breakage_active {
        let geom = breakage::BondGeom {
            r_b: p.r_b,
            area: p.area,
            iben: p.iben,
            jpol: p.jpol,
            l0: p.r0,
        };
        let loads = breakage::BondLoads {
            f_n: f_n_cons,
            f_t_mag: norm(ftan),
            m_bend_mag: norm(bend),
            m_tor_mag: norm(twist),
        };
        let kin = breakage::BondKinematics {
            eps_axial: eps,
            gamma_shear: norm(h.shear_acc) / p.r0,
            kappa_bend: theta_bond / p.r0,
            kappa_tor: norm(h.twist_acc) / p.r0,
        };
        let thr = breakage::BondThresholds { t: h.thresholds };
        if crit.check(&geom, &loads, &kin, &thr).is_some() {
            return None;
        }
    }

    // 8.
    let ft_tot = [
        ftan[0] + p.gamma_t * vt[0],
        ftan[1] + p.gamma_t * vt[1],
        ftan[2] + p.gamma_t * vt[2],
    ];
    let m_el = [
        twist[0] + bend[0] + p.gamma_tor * wn[0] + p.gamma_bend * wt[0],
        twist[1] + bend[1] + p.gamma_tor * wn[1] + p.gamma_bend * wt[1],
        twist[2] + bend[2] + p.gamma_tor * wn[2] + p.gamma_bend * wt[2],
    ];
    let torq1 = cross(r1, ft_tot);
    // torq2 = r2 × (−ft_tot) = r1 × ft_tot.
    Some(Applied {
        force: [
            f_n_cons * n12[0] + ft_tot[0] + p.gamma_n * vn[0],
            f_n_cons * n12[1] + ft_tot[1] + p.gamma_n * vn[1],
            f_n_cons * n12[2] + ft_tot[2] + p.gamma_n * vn[2],
        ],
        torque_i: [m_el[0] + torq1[0], m_el[1] + torq1[1], m_el[2] + torq1[2]],
        torque_j: [torq1[0] - m_el[0], torq1[1] - m_el[1], torq1[2] - m_el[2]],
        strain: eps,
    })
}

/// Computes BPM bond forces and moments for every bond whose lower-tag atom
/// is local. Forces on ghost partners reach their owners through SOIL's
/// reverse communication.
#[allow(clippy::too_many_arguments)]
pub fn bond_force(
    mut atoms: ResMut<Atom>,
    particles: ParticlesWith<
        '_,
        (
            Write<DemAtom>,
            Write<BondHistoryStore>,
            Optional<Write<BondStore>>,
        ),
    >,
    bond_config: Res<BondConfig>,
    breakage: Res<BondBreakage>,
    plasticity: Res<BondPlasticity>,
    mut metrics: ResMut<BondMetrics>,
    mut cache: ResMut<BondPairCache>,
    neighbor: Option<Res<Neighbor>>,
    mut virial: Option<ResMut<VirialStress>>,
    domain: Res<Domain>,
) {
    let step = cache.steps;
    cache.steps += 1;
    let break_min_step = bond_config.breakage_start_step.max(env_break_min_step());
    let breakage_active = step >= break_min_step;
    let log_breaks = std::env::var_os("YARN_BREAK_LOG").is_some();

    particles.with(|(mut dem, mut hist, bonds)| {
        let Some(mut bonds) = bonds else {
            return;
        };
        let nlocal = atoms.nlocal() as usize;
        if bonds.bonds.len() < nlocal || bond_config.bond_radius_ratio <= 0.0 {
            return;
        }

        let rebuilt_neighbors = neighbor.as_ref().map_or(true, |n| n.steps_since_build == 0);
        if !cache.valid
            || rebuilt_neighbors
            || cache.nlocal != nlocal
            || cache.ntotal != atoms.len()
        {
            rebuild_cache(
                &mut cache,
                &atoms,
                &dem,
                &mut bonds,
                &mut hist,
                &bond_config,
                &breakage,
                &domain,
            );
        }
        metrics.missing_partner_skips += cache.missing;

        let dt = atoms.dt;
        let guo = bond_config.formulation == BondFormulation::BpmGuo;
        let model = &plasticity.model;
        let mut broken: Vec<(usize, usize)> = Vec::new();

        let mut crit_type = u32::MAX;
        let mut crit = breakage.resolve(crit_type);
        for p in cache.pairs.iter() {
            let (i, j) = (p.i as usize, p.j as usize);
            let d = min_image(&atoms, &domain, i, j);
            let dist = norm(d);
            if dist < 1e-20 {
                continue;
            }
            // Bonds come in long runs of one type; resolving through the
            // per-type map once per run keeps SipHash out of the loop.
            if p.bond_type != crit_type {
                crit_type = p.bond_type;
                crit = breakage.resolve(crit_type);
            }
            let h = &mut hist.history[i][p.h as usize];
            let applied = if guo {
                guo_bond(
                    p,
                    d,
                    dist,
                    &atoms,
                    &dem,
                    h,
                    model,
                    crit,
                    breakage_active,
                    dt,
                )
            } else {
                dirt_bond(
                    p,
                    d,
                    dist,
                    &atoms,
                    &dem,
                    h,
                    model,
                    crit,
                    breakage_active,
                    dt,
                )
            };
            let Some(a) = applied else {
                broken.push((i, j));
                continue;
            };
            for k in 0..3 {
                atoms.force[i][k] += a.force[k] as soil_core::Accum;
                atoms.force[j][k] -= a.force[k] as soil_core::Accum;
                dem.torque[i][k] += a.torque_i[k];
                dem.torque[j][k] += a.torque_j[k];
            }
            if let Some(v) = virial.as_mut() {
                if v.active {
                    v.add_pair(d[0], d[1], d[2], a.force[0], a.force[1], a.force[2]);
                }
            }
            metrics.strain_sum += a.strain;
            metrics.bond_count += 1;
        }

        if !broken.is_empty() {
            for &(i, j) in &broken {
                let (tag_i, tag_j) = (atoms.tag[i], atoms.tag[j]);
                if log_breaks {
                    let zmid = 0.5 * (atoms.pos[i][2] as f64 + atoms.pos[j][2] as f64);
                    eprintln!("[break] step={step} a={tag_i} b={tag_j} zmid={zmid:.1}");
                }
                for (row, partner) in [(i, tag_j), (j, tag_i)] {
                    if row < bonds.bonds.len() {
                        bonds.bonds[row].retain(|b| b.partner_tag != partner);
                    }
                    if row < hist.history.len() {
                        hist.history[row].retain(|h| h.partner_tag != partner);
                    }
                }
            }
            metrics.bonds_broken_this_step += broken.len();
            metrics.total_bonds_broken += broken.len();
            // Row indices inside the cache are stale now.
            cache.valid = false;
        }
    });
}

/// `BOND_BREAK_MIN_STEP` environment override for `breakage_start_step`
/// (kept for the July yarn scripts that set it).
fn env_break_min_step() -> u64 {
    static MIN: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("BOND_BREAK_MIN_STEP")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproject_keeps_magnitude_and_direction_class() {
        let n = [0.0, 0.0, 1.0];
        let mut v = [3.0, 0.0, 4.0];
        reproject(&mut v, n, false);
        assert!((norm(v) - 5.0).abs() < 1e-12);
        assert!(v[2].abs() < 1e-12 && v[0] > 0.0);
        let mut t = [0.0, 3.0, -4.0];
        reproject(&mut t, n, true);
        assert!((t[2] + 5.0).abs() < 1e-12, "sign of the axial part is kept");
        let mut z = [0.0, 0.0, 2.0];
        reproject(&mut z, n, false);
        assert_eq!(z, [0.0, 0.0, 2.0], "zero projection leaves v unchanged");
    }
}
