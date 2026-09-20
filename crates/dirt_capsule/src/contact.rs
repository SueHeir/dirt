//! Capsule contact: Hertz normal + Mindlin history-tangential + Tsuji damping.
//!
//! The contact geometry is the only capsule-specific part — the closest points
//! between the two centerline segments (a sphere is the degenerate zero-length
//! segment) give a penetration `delta`, effective radius `r_eff`, and unit normal
//! exactly as sphere–sphere contact does. The force law then reuses the **same
//! precomputed material coefficients** as `dirt_granular` (`e_eff_ij`, `g_eff_ij`,
//! `beta_ij`, `friction_ij`, `limit_damping`), so a capsule gas and a rod-of-
//! spheres gas of the same material see an identical contact law.
//!
//! Unlike sphere contact, the contact point of a capsule is generally **off the
//! line of centers**, so the normal force exerts a torque about the COM; this
//! module applies the *total* contact force at the contact point and takes
//! `arm × F` for the torque, which reduces to the sphere result when the contact
//! is central.
//!
//! Tangential history is kept in a [`CapsuleContactHistory`] resource keyed by the
//! sorted tag pair — correct for the single-rank periodic gas used in the
//! rod-vs-capsule benchmark. (Multi-rank runs would need the history forwarded
//! like `dirt_granular::ContactHistoryStore`; that is future work.)

use std::collections::{HashMap, HashSet};

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use soil_core::{
    Accum, Atom, BondStore, Neighbor, Optional, ParticleSimScheduleSet, ParticlesWith, Read, Write,
};

use dirt_atom::{DemAtom, MaterialTable, SQRT_5_6};
use dirt_schedule::CONTACT_FORCE;

use crate::geometry::{capsule_overlap, closest_segment_segment};
use crate::{capsule_endpoints, is_capsule, CapsuleAtom};

const TANGENTIAL_EPS: f64 = 1e-12;

/// Persistent tangential-spring history keyed by the sorted `(tag_lo, tag_hi)`
/// pair. The stored spring is in the `lo → hi` orientation (canonical), matching
/// the sign convention `dirt_granular` uses.
#[derive(Default)]
pub struct CapsuleContactHistory {
    springs: HashMap<(u32, u32), [f64; 3]>,
}

/// Plugin adding capsule–capsule / capsule–sphere contact.
///
/// Intended for capsule-only (or capsule-dominated) simulations: it processes
/// every neighbor pair in which at least one atom is a capsule. Do **not** also
/// add it for pure-sphere pairs already covered by `HertzMindlinContactPlugin`.
pub struct CapsuleContactPlugin;

impl Plugin for CapsuleContactPlugin {
    fn build(&self, app: &mut App) {
        app.add_resource(CapsuleContactHistory::default());
        app.add_update_system(
            capsule_contact_force.label(CONTACT_FORCE),
            ParticleSimScheduleSet::Force,
        );
    }
}

/// Cross product `a × b`.
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Result of the pure force law for one contact.
struct PairForce {
    /// Total contact force on atom `i` (space frame).
    force_i: [f64; 3],
    /// Updated canonical (lo→hi) tangential spring.
    spring_canonical: [f64; 3],
}

/// Hertz normal + Mindlin history tangential + Tsuji damping for one contact,
/// given geometry, material coefficients, contact-point kinematics, and the
/// stored canonical spring. Mirrors the `dirt_granular` Hertz path (no adhesion /
/// cohesion / MDR), which is the model the rod-vs-capsule benchmark uses.
#[allow(clippy::too_many_arguments)]
fn hertz_mindlin_pair(
    delta: f64,
    r_eff: f64,
    normal: [f64; 3], // unit, points i → j
    v_rel: [f64; 3],  // v_contact_j − v_contact_i
    e_eff: f64,
    g_eff: f64,
    beta: f64,
    mu: f64,
    limit_damping: bool,
    m_r: f64,
    sign: f64, // +1 if tag_i < tag_j else −1 (maps canonical spring to i→j frame)
    spring_canonical: [f64; 3],
    dt: f64,
) -> PairForce {
    let n = normal;
    let sdr = (delta * r_eff).sqrt();
    let s_n = 2.0 * e_eff * sdr;
    let k_n = 4.0 / 3.0 * e_eff * sdr;
    let k_t = 8.0 * g_eff * sdr;

    let v_n = v_rel[0] * n[0] + v_rel[1] * n[1] + v_rel[2] * n[2];

    // Normal force (repulsive positive along n): Hertz spring − Tsuji damping.
    let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
    let mut f_n = k_n * delta - f_diss_n;
    if limit_damping {
        f_n = f_n.max(0.0);
    }
    let f_t_max = mu * f_n.abs();

    // Tangential relative velocity (remove normal component).
    let vt = [
        v_rel[0] - v_n * n[0],
        v_rel[1] - v_n * n[1],
        v_rel[2] - v_n * n[2],
    ];

    // Bring stored spring into the i→j frame.
    let mut s = [
        sign * spring_canonical[0],
        sign * spring_canonical[1],
        sign * spring_canonical[2],
    ];
    // Rotate spring into the current tangent plane (drop normal component).
    let s_dot_n = s[0] * n[0] + s[1] * n[1] + s[2] * n[2];
    for d in 0..3 {
        s[d] -= s_dot_n * n[d];
    }
    // Integrate tangential velocity into the spring displacement.
    for d in 0..3 {
        s[d] += vt[d] * dt;
    }
    // Coulomb cap on the spring: |k_t s| ≤ μ|F_n|.
    let s_mag = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    let f_spring = k_t * s_mag;
    if f_spring > f_t_max && f_spring > TANGENTIAL_EPS {
        let scale = f_t_max / f_spring;
        for d in 0..3 {
            s[d] *= scale;
        }
    }

    // Tangential force: elastic spring + Tsuji tangential damping.
    let gamma_t = 2.0 * SQRT_5_6 * beta * (k_t * m_r).sqrt();
    let mut ft = [
        k_t * s[0] + gamma_t * vt[0],
        k_t * s[1] + gamma_t * vt[1],
        k_t * s[2] + gamma_t * vt[2],
    ];
    // Coulomb cap on the total tangential force.
    let ft_mag = (ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2]).sqrt();
    if ft_mag > f_t_max && ft_mag > TANGENTIAL_EPS {
        let scale = f_t_max / ft_mag;
        for d in 0..3 {
            ft[d] *= scale;
        }
    }

    // Total force on i: normal pushes i away from j (−n), plus tangential.
    let force_i = [
        -f_n * n[0] + ft[0],
        -f_n * n[1] + ft[1],
        -f_n * n[2] + ft[2],
    ];
    // Store the spring back in canonical (lo→hi) orientation.
    let spring_canonical = [sign * s[0], sign * s[1], sign * s[2]];
    PairForce {
        force_i,
        spring_canonical,
    }
}

/// Compute capsule contact forces and torques for all neighbor pairs that
/// involve at least one capsule.
fn capsule_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    mut history: ResMut<CapsuleContactHistory>,
    material_table: Res<MaterialTable>,
    particles: ParticlesWith<'_, (Write<DemAtom>, Read<CapsuleAtom>, Optional<Read<BondStore>>)>,
) {
    let dt = atoms.dt;
    let newton = neighbor.newton;
    let nlocal = atoms.nlocal() as usize;

    particles.with(|(mut dem, caps, bonds)| {
        let mut seen: HashSet<(u32, u32)> = HashSet::new();

        for (i, j) in neighbor.pairs(nlocal) {
            let i_cap = is_capsule(&caps, i);
            let j_cap = is_capsule(&caps, j);
            if !i_cap && !j_cap {
                continue; // pure sphere pair — not ours
            }
            // Only the directly/second-neighbor bonded capsules are excluded.
            // Distant segments of the same fiber must retain self-contact.
            if dirt_atom::same_body(&dem, i, j)
                || bonds
                    .as_ref()
                    .map(|store| store.are_excluded(i, j, &atoms.tag))
                    .unwrap_or(false)
            {
                continue;
            }

            let pos_i = [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ];
            let pos_j = [
                atoms.pos[j][0] as f64,
                atoms.pos[j][1] as f64,
                atoms.pos[j][2] as f64,
            ];
            let r_i = dem.radius[i];
            let r_j = dem.radius[j];

            // Centerline segments (a sphere degenerates to its center point).
            let (ai0, ai1) = if i_cap {
                capsule_endpoints(pos_i, caps.axis[i], caps.half_length[i])
            } else {
                (pos_i, pos_i)
            };
            let (aj0, aj1) = if j_cap {
                capsule_endpoints(pos_j, caps.axis[j], caps.half_length[j])
            } else {
                (pos_j, pos_j)
            };

            let (_s, _t, ca, cb) = closest_segment_segment(ai0, ai1, aj0, aj1);
            let (delta, r_eff, n) = match capsule_overlap(ca, cb, r_i, r_j) {
                Some(v) => v,
                None => continue,
            };

            // Contact point: midway between the two capsule surfaces along n.
            let cp = [
                ca[0] + (r_i - 0.5 * delta) * n[0],
                ca[1] + (r_i - 0.5 * delta) * n[1],
                ca[2] + (r_i - 0.5 * delta) * n[2],
            ];
            let arm_i = [cp[0] - pos_i[0], cp[1] - pos_i[1], cp[2] - pos_i[2]];
            let arm_j = [cp[0] - pos_j[0], cp[1] - pos_j[1], cp[2] - pos_j[2]];

            // Contact-point velocities (rigid body: v + ω × arm).
            let vi = [
                atoms.vel[i][0] as f64,
                atoms.vel[i][1] as f64,
                atoms.vel[i][2] as f64,
            ];
            let vj = [
                atoms.vel[j][0] as f64,
                atoms.vel[j][1] as f64,
                atoms.vel[j][2] as f64,
            ];
            let wi = cross(dem.omega[i], arm_i);
            let wj = cross(dem.omega[j], arm_j);
            let v_rel = [
                (vj[0] + wj[0]) - (vi[0] + wi[0]),
                (vj[1] + wj[1]) - (vi[1] + wi[1]),
                (vj[2] + wj[2]) - (vi[2] + wi[2]),
            ];

            // Material coefficients (shared with dirt_granular).
            let mi = atoms.atom_type[i] as usize;
            let mj = atoms.atom_type[j] as usize;
            let e_eff = material_table.e_eff_ij[mi][mj];
            let g_eff = material_table.g_eff_ij[mi][mj];
            let beta = material_table.beta_ij[mi][mj];
            let mu = material_table.friction_ij[mi][mj];

            // Reduced mass (fall back to real mass if inv_mass is zeroed).
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

            let tag_i = atoms.tag[i];
            let tag_j = atoms.tag[j];
            let (key, sign) = if tag_i < tag_j {
                ((tag_i, tag_j), 1.0)
            } else {
                ((tag_j, tag_i), -1.0)
            };
            let spring0 = history.springs.get(&key).copied().unwrap_or([0.0; 3]);

            let out = hertz_mindlin_pair(
                delta,
                r_eff,
                n,
                v_rel,
                e_eff,
                g_eff,
                beta,
                mu,
                material_table.limit_damping,
                m_r,
                sign,
                spring0,
                dt,
            );
            history.springs.insert(key, out.spring_canonical);
            seen.insert(key);

            let f_i = out.force_i;
            let f_j = [-f_i[0], -f_i[1], -f_i[2]];
            let t_i = cross(arm_i, f_i);
            let t_j = cross(arm_j, f_j);

            for d in 0..3 {
                atoms.force[i][d] += f_i[d] as Accum;
                dem.torque[i][d] += t_i[d];
            }
            if newton {
                for d in 0..3 {
                    atoms.force[j][d] += f_j[d] as Accum;
                    dem.torque[j][d] += t_j[d];
                }
            }
        }

        // Drop springs for contacts that ended this step (no longer neighbors).
        history.springs.retain(|k, _| seen.contains(k));
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A central capsule–capsule collision (contact on the line of centers) must
    /// reduce to the plain Hertz normal force with no torque and no tangential.
    #[test]
    fn central_normal_matches_hertz() {
        let (delta, r_eff) = (0.001, 0.25);
        let e_eff = 3.0e7;
        // Head-on approach: v_rel purely along +n, magnitude 1.0.
        let n = [1.0, 0.0, 0.0];
        let out = hertz_mindlin_pair(
            delta,
            r_eff,
            n,
            [1.0, 0.0, 0.0],
            e_eff,
            1.0e7,
            0.0,
            0.3,
            true,
            1e-6,
            1.0,
            [0.0; 3],
            1e-7,
        );
        // Expected Hertz spring minus a positive damping term (β=0 ⇒ no damping).
        let sdr = (delta * r_eff).sqrt();
        let k_n = 4.0 / 3.0 * e_eff * sdr;
        let expected = -k_n * delta; // force on i is −f_n n
        assert!((out.force_i[0] - expected).abs() / expected.abs() < 1e-12);
        assert!(out.force_i[1].abs() < 1e-12 && out.force_i[2].abs() < 1e-12);
    }

    /// With β = 0 and a purely normal approach, the tangential spring stays zero.
    #[test]
    fn no_tangential_on_pure_normal() {
        let out = hertz_mindlin_pair(
            0.001,
            0.25,
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -0.5],
            3.0e7,
            1.0e7,
            0.2,
            0.3,
            true,
            1e-6,
            1.0,
            [0.0; 3],
            1e-7,
        );
        assert!(out.spring_canonical.iter().all(|x| x.abs() < 1e-15));
    }
}
