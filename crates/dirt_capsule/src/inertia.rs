//! Analytic mass and principal moments of inertia of a homogeneous
//! sphero-cylinder (capsule): a cylinder of radius `r` and length `L` capped by
//! two hemispheres of radius `r`, taken about its center of mass.
//!
//! With the symmetry (body) axis along **z**, the capsule is axisymmetric so the
//! two transverse moments are equal (`I_x = I_y`) and the body frame is already a
//! principal frame. Splitting the mass into the cylinder (`m_c = ρ π r² L`) and
//! the two hemispheres, which together form one sphere (`m_s = ρ (4/3) π r³`):
//!
//! ```text
//! I_z = m_c (r²/2)             + m_s (2/5) r²
//! I_x = m_c (L²/12 + r²/4)     + m_s [ (2/5) r² + (1/4) L² + (3/8) r L ]
//! ```
//!
//! The cap transverse term is the two hemispheres' moment about the capsule
//! center: each hemisphere's COM sits `3r/8` from its flat face, giving
//! `m_s[(2/5)r² + L²/4 + (3/8)rL]` after the parallel-axis shift. Limits check
//! out: `L → 0` recovers a solid sphere `(2/5) m_s r²` on every axis; `r → 0`
//! recovers a thin rod `m_c L²/12` transverse.

use std::f64::consts::PI;

/// Mass and principal moments of a capsule of radius `r`, cylindrical length `L`,
/// and density `rho`, about its center of mass, with the symmetry axis along z.
///
/// Returns `(mass, [I_x, I_y, I_z])` with `I_x == I_y`.
pub fn capsule_mass_and_moments(radius: f64, length: f64, density: f64) -> (f64, [f64; 3]) {
    let r = radius;
    let l = length;
    let m_c = density * PI * r * r * l; // cylinder
    let m_s = density * (4.0 / 3.0) * PI * r * r * r; // both hemispheres = one sphere
    let mass = m_c + m_s;

    let i_z = m_c * (0.5 * r * r) + m_s * (0.4 * r * r);
    let i_x =
        m_c * (l * l / 12.0 + 0.25 * r * r) + m_s * (0.4 * r * r + 0.25 * l * l + 0.375 * r * l);

    (mass, [i_x, i_x, i_z])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sphere_limit() {
        // L = 0 → solid sphere: mass = 4/3 π r³ ρ, all moments = 2/5 m r².
        let (m, i) = capsule_mass_and_moments(0.5, 0.0, 2500.0);
        let expected_m = 2500.0 * (4.0 / 3.0) * PI * 0.5f64.powi(3);
        let expected_i = 0.4 * expected_m * 0.25;
        assert!((m - expected_m).abs() / expected_m < 1e-12);
        for a in 0..3 {
            assert!((i[a] - expected_i).abs() / expected_i < 1e-12, "axis {a}");
        }
    }

    #[test]
    fn thin_rod_limit() {
        // r → 0 (tiny) → transverse moment approaches thin-rod m_c L²/12,
        // and axial moment → 0.
        let r = 1e-6;
        let l = 1.0;
        let rho = 2500.0;
        let (_m, i) = capsule_mass_and_moments(r, l, rho);
        let m_c = rho * PI * r * r * l;
        let rod_transverse = m_c * l * l / 12.0;
        assert!((i[0] - rod_transverse).abs() / rod_transverse < 1e-3);
        assert!(i[2] < i[0] * 1e-6);
    }

    #[test]
    fn axisymmetric_transverse_equal() {
        let (_m, i) = capsule_mass_and_moments(0.5, 2.0, 1000.0);
        assert_eq!(i[0], i[1]);
        // Elongated capsule: transverse > axial.
        assert!(i[0] > i[2]);
    }

    #[test]
    fn matches_montecarlo() {
        // Cross-check the analytic transverse/axial moments against a crude
        // Monte-Carlo integration over the capsule volume.
        let r = 0.4;
        let l = 1.2;
        let rho = 1500.0;
        let (mass, i) = capsule_mass_and_moments(r, l, rho);

        // Bounding box: x,y in [-r,r], z in [-(l/2+r), l/2+r].
        let hz = 0.5 * l + r;
        let (bx, bz) = (r, hz);
        let n = 4_000_000u64;
        // deterministic LCG sampler
        let mut state = 0x2545F4914F6CDD1Du64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / ((1u64 << 53) as f64)
        };
        let mut hits = 0u64;
        let (mut sx, mut sz) = (0.0f64, 0.0f64); // Σ(y²+z²) and Σ(x²+y²) accumulators
        for _ in 0..n {
            let x = (next() * 2.0 - 1.0) * bx;
            let y = (next() * 2.0 - 1.0) * bx;
            let z = (next() * 2.0 - 1.0) * bz;
            let inside = if z.abs() <= 0.5 * l {
                x * x + y * y <= r * r // cylinder body
            } else {
                let dz = z.abs() - 0.5 * l;
                x * x + y * y + dz * dz <= r * r // hemispherical cap
            };
            if inside {
                hits += 1;
                sx += y * y + z * z; // transverse (about x)
                sz += x * x + y * y; // axial (about z)
            }
        }
        let box_vol = (2.0 * bx) * (2.0 * bx) * (2.0 * bz);
        let dv = box_vol / n as f64;
        let mc_mass = rho * dv * hits as f64;
        let mc_ix = rho * dv * sx;
        let mc_iz = rho * dv * sz;

        assert!(
            (mc_mass - mass).abs() / mass < 0.02,
            "mass {mc_mass} vs {mass}"
        );
        assert!(
            (mc_ix - i[0]).abs() / i[0] < 0.02,
            "I_x {mc_ix} vs {}",
            i[0]
        );
        assert!(
            (mc_iz - i[2]).abs() / i[2] < 0.02,
            "I_z {mc_iz} vs {}",
            i[2]
        );
    }
}
