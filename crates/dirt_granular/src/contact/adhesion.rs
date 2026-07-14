/// Willett et al. (2000) dimensionless closed-form pendular bridge force for
/// two spheres, using the common equal-sphere effective-radius extension:
///
/// `F = 2*pi*R*gamma*cos(theta) / (1 + 1.05*s_hat + 2.5*s_hat^2)`,
/// `s_hat = s*sqrt(R*/V)`, active for `0 <= s <= s_rupture`.
///
/// Returns a positive tensile magnitude; the contact force path subtracts it
/// from the normal force so it acts attractively.
pub fn willett2000_liquid_bridge_force(
    separation: f64,
    r_eff: f64,
    volume: f64,
    surface_tension: f64,
    contact_angle: f64,
    rupture_distance: f64,
) -> f64 {
    if separation < 0.0 || r_eff <= 0.0 || volume <= 0.0 || surface_tension <= 0.0 {
        return 0.0;
    }
    let rupture = if rupture_distance > 0.0 {
        rupture_distance
    } else {
        // Lian, Thornton & Adams' volume-scaled rupture estimate, also used
        // with the Willett force in DEM implementations.
        (1.0 + 0.5 * contact_angle) * volume.cbrt()
    };
    if separation > rupture {
        return 0.0;
    }
    let s_hat = separation * (r_eff / volume).sqrt();
    let denom = 1.0 + 1.05 * s_hat + 2.5 * s_hat * s_hat;
    let wetting = contact_angle.cos().max(0.0);
    2.0 * std::f64::consts::PI * r_eff * surface_tension * wetting / denom
}
