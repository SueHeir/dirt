pub(super) fn mdr_nonadhesive_force(delta: f64, a_inv: f64, e_eff: f64, a: f64, b: f64) -> f64 {
    if delta <= 0.0 || a <= 0.0 || b <= 0.0 {
        return 0.0;
    }
    let x = (delta * a_inv).clamp(0.0, 1.0);
    let root = (x - x * x).max(0.0).sqrt();
    0.25 * e_eff * a * b * ((1.0 - 2.0 * x).acos() - (2.0 - 4.0 * x) * root)
}

pub(super) fn mdr_round_up_negative_epsilon(value: f64) -> f64 {
    if value < 0.0 && value > -1.0e-20 {
        0.0
    } else {
        value
    }
}
