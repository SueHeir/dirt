use super::*;
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

pub(super) fn mdr_normal_force(
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
