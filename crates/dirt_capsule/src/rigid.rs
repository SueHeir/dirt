//! Free-atom rigid-body rotational integration for anisotropic particles.
//!
//! A capsule has a non-spherical (though axisymmetric) inertia tensor, which the
//! scalar per-atom integrator in `dirt_granular` cannot advance. This module
//! provides the angular-momentum / Richardson-extrapolated quaternion update used
//! for a single free particle, mirroring the validated `dirt_clump::body` scheme:
//!
//! - **Angular momentum `L` is the integrated state** — half-kicked by the torque
//!   each half-step (`L += ½ dt τ`); angular velocity `ω` is *derived* from `L`
//!   via [`angmom_to_omega`], never integrated directly. Stable for anisotropic
//!   inertia where integrating Euler's equation for `ω` is fragile.
//! - **Richardson-extrapolated quaternion update** — one full step and two half
//!   steps, `q ← 2 q_half − q_full`, renormalized (second order).
//!
//! The math is intentionally duplicated from `dirt_clump::body` rather than shared
//! through a dependency, so this crate stays self-contained.

/// Hamilton product of two quaternions `[w, x, y, z]`.
#[inline]
pub fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

/// Conjugate (inverse for unit quaternions): `[w, -x, -y, -z]`.
#[inline]
pub fn quat_conj(q: [f64; 4]) -> [f64; 4] {
    [q[0], -q[1], -q[2], -q[3]]
}

/// Normalize a quaternion to unit length (identity on a (near-)zero input).
#[inline]
pub fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if n > 1e-30 {
        let inv = 1.0 / n;
        [q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv]
    } else {
        [1.0, 0.0, 0.0, 0.0]
    }
}

/// Rotate a 3-vector by quaternion `q = [w, x, y, z]`.
#[inline]
pub fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let (w, qx, qy, qz) = (q[0], q[1], q[2], q[3]);
    let cx = qy * v[2] - qz * v[1];
    let cy = qz * v[0] - qx * v[2];
    let cz = qx * v[1] - qy * v[0];
    [
        v[0] + 2.0 * (w * cx + qy * cz - qz * cy),
        v[1] + 2.0 * (w * cy + qz * cx - qx * cz),
        v[2] + 2.0 * (w * cz + qx * cy - qy * cx),
    ]
}

/// Columns `(ex, ey, ez)` of the rotation matrix for quaternion `q`.
#[inline]
fn quat_to_rotation_columns(q: [f64; 4]) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let (w, x, y, z) = (q[0], q[1], q[2], q[3]);
    let ex = [
        w * w + x * x - y * y - z * z,
        2.0 * (x * y + w * z),
        2.0 * (x * z - w * y),
    ];
    let ey = [
        2.0 * (x * y - w * z),
        w * w - x * x + y * y - z * z,
        2.0 * (y * z + w * x),
    ];
    let ez = [
        2.0 * (x * z + w * y),
        2.0 * (y * z - w * x),
        w * w - x * x - y * y + z * z,
    ];
    (ex, ey, ez)
}

/// Convert space-frame angular momentum to space-frame angular velocity.
///
/// `quat_principal_to_space` must be `quaternion * principal_axes` (composed by
/// the caller). Projects `L` into the principal frame, divides by `moments`, and
/// rotates the result back to space. Mirrors `MathExtra::angmom_to_omega`.
#[inline]
pub fn angmom_to_omega(
    angmom: [f64; 3],
    quat_principal_to_space: [f64; 4],
    moments: [f64; 3],
) -> [f64; 3] {
    let (ex, ey, ez) = quat_to_rotation_columns(quat_principal_to_space);
    let lpx = ex[0] * angmom[0] + ex[1] * angmom[1] + ex[2] * angmom[2];
    let lpy = ey[0] * angmom[0] + ey[1] * angmom[1] + ey[2] * angmom[2];
    let lpz = ez[0] * angmom[0] + ez[1] * angmom[1] + ez[2] * angmom[2];
    let wpx = if moments[0] > 1e-30 {
        lpx / moments[0]
    } else {
        0.0
    };
    let wpy = if moments[1] > 1e-30 {
        lpy / moments[1]
    } else {
        0.0
    };
    let wpz = if moments[2] > 1e-30 {
        lpz / moments[2]
    } else {
        0.0
    };
    [
        ex[0] * wpx + ey[0] * wpy + ez[0] * wpz,
        ex[1] * wpx + ey[1] * wpy + ez[1] * wpz,
        ex[2] * wpx + ey[2] * wpy + ez[2] * wpz,
    ]
}

/// Quaternion kinematic derivative `dq/dt` components for `dq/dt = ½ ω⊗q`.
#[inline]
fn vecquat(w: [f64; 3], q: [f64; 4]) -> [f64; 4] {
    [
        -w[0] * q[1] - w[1] * q[2] - w[2] * q[3],
        q[0] * w[0] + w[1] * q[3] - w[2] * q[2],
        q[0] * w[1] + w[2] * q[1] - w[0] * q[3],
        q[0] * w[2] + w[0] * q[2] - w[1] * q[1],
    ]
}

/// Richardson-extrapolated quaternion advance over `dtq` (second order).
///
/// `q` is body→space; `principal_axes` is body→principal (composed internally for
/// the half-step omega recompute). Mirrors `FixMultisphere::richardson`.
pub fn richardson(
    q: &mut [f64; 4],
    angmom: [f64; 3],
    omega: [f64; 3],
    moments: [f64; 3],
    dtq: f64,
    principal_axes: [f64; 4],
) {
    let half = 0.5 * dtq;
    let dq = vecquat(omega, *q);
    let q_full = quat_normalize([
        q[0] + dtq * dq[0],
        q[1] + dtq * dq[1],
        q[2] + dtq * dq[2],
        q[3] + dtq * dq[3],
    ]);
    let mut q_half = quat_normalize([
        q[0] + half * dq[0],
        q[1] + half * dq[1],
        q[2] + half * dq[2],
        q[3] + half * dq[3],
    ]);
    let omega_half = angmom_to_omega(angmom, quat_mul(q_half, principal_axes), moments);
    let dq2 = vecquat(omega_half, q_half);
    q_half = quat_normalize([
        q_half[0] + half * dq2[0],
        q_half[1] + half * dq2[1],
        q_half[2] + half * dq2[2],
        q_half[3] + half * dq2[3],
    ]);
    *q = quat_normalize([
        2.0 * q_half[0] - q_full[0],
        2.0 * q_half[1] - q_full[1],
        2.0 * q_half[2] - q_full[2],
        2.0 * q_half[3] - q_full[3],
    ]);
}

/// Rigid-body rotational state carried per capsule atom (a view over DemAtom
/// columns plus the fixed principal data from `CapsuleAtom`).
///
/// `principal_axes` is body→principal; for an axisymmetric capsule with its
/// symmetry axis along the body z-axis this is the identity `[1,0,0,0]`, but the
/// integrator supports any principal frame.
pub struct RigidRotState<'a> {
    /// Orientation quaternion (body→space), advanced in place.
    pub quaternion: &'a mut [f64; 4],
    /// Angular velocity (space frame), derived from `angmom`.
    pub omega: &'a mut [f64; 3],
    /// Angular momentum (space frame), the integrated state.
    pub angmom: &'a mut [f64; 3],
    /// Principal moments `[I_x, I_y, I_z]`.
    pub moments: [f64; 3],
    /// Body→principal rotation quaternion.
    pub principal_axes: [f64; 4],
}

impl RigidRotState<'_> {
    /// Initial half-step: half-kick `L` from torque, derive `ω`, Richardson-advance
    /// the quaternion, then recompute `ω` at the new orientation.
    pub fn integrate_initial(&mut self, torque: [f64; 3], dt: f64) {
        let half = 0.5 * dt;
        for d in 0..3 {
            self.angmom[d] += half * torque[d];
        }
        let q_ps = quat_mul(*self.quaternion, self.principal_axes);
        *self.omega = angmom_to_omega(*self.angmom, q_ps, self.moments);
        richardson(
            self.quaternion,
            *self.angmom,
            *self.omega,
            self.moments,
            half,
            self.principal_axes,
        );
        let q_ps = quat_mul(*self.quaternion, self.principal_axes);
        *self.omega = angmom_to_omega(*self.angmom, q_ps, self.moments);
    }

    /// Final half-step: half-kick `L` from the new torque and recompute `ω`.
    pub fn integrate_final(&mut self, torque: [f64; 3], dt: f64) {
        let half = 0.5 * dt;
        for d in 0..3 {
            self.angmom[d] += half * torque[d];
        }
        let q_ps = quat_mul(*self.quaternion, self.principal_axes);
        *self.omega = angmom_to_omega(*self.angmom, q_ps, self.moments);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Free asymmetric top: torque-free rotation must conserve rotational KE.
    #[test]
    fn torque_free_conserves_energy() {
        let moments = [1.0, 1.5, 2.0]; // anisotropic (like an elongated capsule)
        let mut q = [1.0, 0.0, 0.0, 0.0];
        let mut omega = [0.7, 1.3, 0.9]; // tumbling about all axes
                                         // L = I ω in the principal (= body, identity axes) frame, rotated to space.
        let mut angmom = [
            moments[0] * omega[0],
            moments[1] * omega[1],
            moments[2] * omega[2],
        ];
        let pa = [1.0, 0.0, 0.0, 0.0];

        let ke = |om: [f64; 3], qq: [f64; 4]| {
            let ob = quat_rotate(quat_conj(qq), om);
            0.5 * (moments[0] * ob[0] * ob[0]
                + moments[1] * ob[1] * ob[1]
                + moments[2] * ob[2] * ob[2])
        };
        let ke0 = ke(omega, q);

        let dt = 1e-4;
        for _ in 0..20000 {
            let mut st = RigidRotState {
                quaternion: &mut q,
                omega: &mut omega,
                angmom: &mut angmom,
                moments,
                principal_axes: pa,
            };
            st.integrate_initial([0.0; 3], dt);
            st.integrate_final([0.0; 3], dt);
        }
        let ke1 = ke(omega, q);
        let rel = (ke1 - ke0).abs() / ke0;
        assert!(rel < 1e-2, "rot KE drift {rel} (ke0={ke0}, ke1={ke1})");

        // Angular momentum magnitude (space frame) is exactly conserved (torque-free).
        let l0 =
            (1.0f64 * 0.7 * 1.0f64 * 0.7 + 1.5 * 1.3 * 1.5 * 1.3 + 2.0 * 0.9 * 2.0 * 0.9).sqrt();
        let l1 = (angmom[0] * angmom[0] + angmom[1] * angmom[1] + angmom[2] * angmom[2]).sqrt();
        assert!((l1 - l0).abs() / l0 < 1e-10, "|L| drift: {l0} -> {l1}");
    }

    #[test]
    fn pure_torque_spins_up_about_axis() {
        let moments = [2.0, 3.0, 4.0];
        let mut q = [1.0, 0.0, 0.0, 0.0];
        let mut omega = [0.0; 3];
        let mut angmom = [0.0; 3];
        let dt = 1e-5;
        let tau = [10.0, 0.0, 0.0];
        let mut st = RigidRotState {
            quaternion: &mut q,
            omega: &mut omega,
            angmom: &mut angmom,
            moments,
            principal_axes: [1.0, 0.0, 0.0, 0.0],
        };
        st.integrate_initial(tau, dt);
        st.integrate_final(tau, dt);
        // After one full step ω_x ≈ dt τ_x / I_x.
        assert!((omega[0] - dt * 10.0 / 2.0).abs() < 1e-10);
    }
}
