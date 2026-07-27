//! Pure geometry for capsule (sphero-cylinder) contact — no engine types, so it
//! is trivially unit-testable in isolation. All vectors are `[f64; 3]`.
//!
//! A capsule is a centerline segment `[p0, p1]` swept with radius `r`. Two
//! capsules are in contact when the closest points between their centerlines are
//! nearer than the sum of the radii; a capsule and a sphere are the degenerate
//! case where one segment collapses to a point.

/// A 3-vector in world space.
pub type V3 = [f64; 3];

const EPS: f64 = 1e-12;

/// Component-wise difference `a - b`.
#[inline]
pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Component-wise sum `a + b`.
#[inline]
pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scalar multiple `s · a`.
#[inline]
pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product `a · b`.
#[inline]
pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean length `|a|`.
#[inline]
pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Closest points between segments `[pa0, pa1]` and `[pb0, pb1]`.
///
/// Returns `(s, t, closest_on_A, closest_on_B)` with `s, t ∈ [0, 1]` the
/// parametric positions along A and B. Follows the clamped formulation from
/// Ericson, *Real-Time Collision Detection*. The near-parallel branch
/// (`denom ≈ 0`) resolves deterministically to `s = 0` so an exactly-parallel
/// pair does not jitter.
pub fn closest_segment_segment(pa0: V3, pa1: V3, pb0: V3, pb1: V3) -> (f64, f64, V3, V3) {
    let d1 = sub(pa1, pa0); // direction of A
    let d2 = sub(pb1, pb0); // direction of B
    let r = sub(pa0, pb0);
    let a = dot(d1, d1); // |A|^2
    let e = dot(d2, d2); // |B|^2
    let f = dot(d2, r);

    // Both segments degenerate to points.
    if a <= EPS && e <= EPS {
        return (0.0, 0.0, pa0, pb0);
    }
    // A degenerate to a point.
    if a <= EPS {
        let t = (f / e).clamp(0.0, 1.0);
        return (0.0, t, pa0, add(pb0, scale(d2, t)));
    }

    let c = dot(d1, r);
    let (s, t);
    if e <= EPS {
        // B degenerate to a point.
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else {
        let b = dot(d1, d2);
        let denom = a * e - b * b; // always >= 0
        let s0 = if denom > EPS {
            ((b * f - c * e) / denom).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let mut t0 = (b * s0 + f) / e;
        if t0 < 0.0 {
            t0 = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else if t0 > 1.0 {
            t0 = 1.0;
            s = ((b - c) / a).clamp(0.0, 1.0);
        } else {
            s = s0;
        }
        t = t0;
    }
    (s, t, add(pa0, scale(d1, s)), add(pb0, scale(d2, t)))
}

/// Closest point on segment `[p0, p1]` to point `q`.
///
/// Returns `(s, closest_on_segment)` with `s ∈ [0, 1]`. This is the
/// capsule–sphere degenerate case of [`closest_segment_segment`].
pub fn closest_segment_point(p0: V3, p1: V3, q: V3) -> (f64, V3) {
    let d = sub(p1, p0);
    let len2 = dot(d, d);
    if len2 <= EPS {
        return (0.0, p0);
    }
    let s = (dot(sub(q, p0), d) / len2).clamp(0.0, 1.0);
    (s, add(p0, scale(d, s)))
}

/// Overlap of two capsules from their centerline closest points and cross-
/// section radii.
///
/// Returns `None` when the capsules are separated. On contact returns
/// `(delta, r_eff, normal)`:
/// - `delta > 0` is the penetration depth (feeds the Hertz law unchanged),
/// - `r_eff = ra·rb / (ra + rb)` — the same effective radius as sphere–sphere,
///   because at the closest point each capsule presents its circular section,
/// - `normal` is the unit contact normal pointing from A toward B.
pub fn capsule_overlap(ca: V3, cb: V3, ra: f64, rb: f64) -> Option<(f64, f64, V3)> {
    let d = sub(cb, ca);
    let dist = norm(d);
    if dist <= EPS {
        return None;
    }
    let delta = (ra + rb) - dist;
    if delta <= 0.0 {
        return None;
    }
    let r_eff = ra * rb / (ra + rb);
    Some((delta, r_eff, scale(d, 1.0 / dist)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perpendicular_crossing_meets_at_midpoints() {
        // A along x at z=0, B along y at z=1, crossing over the origin.
        let (s, t, ca, cb) = closest_segment_segment(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 1.0],
            [0.0, 1.0, 1.0],
        );
        assert!((s - 0.5).abs() < 1e-9);
        assert!((t - 0.5).abs() < 1e-9);
        assert!(norm(sub(ca, [0.0, 0.0, 0.0])) < 1e-9);
        assert!(norm(sub(cb, [0.0, 0.0, 1.0])) < 1e-9);
    }

    #[test]
    fn segment_point_projects_onto_interior() {
        let (s, cp) = closest_segment_point([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 3.0, 0.0]);
        assert!((s - 0.5).abs() < 1e-9);
        assert!(norm(sub(cp, [1.0, 0.0, 0.0])) < 1e-9);
    }

    #[test]
    fn overlap_detects_touch() {
        // Closest points 1.5 apart, radii sum 2.0 -> penetration 0.5.
        let hit = capsule_overlap([0.0, 0.0, 0.0], [0.0, 0.0, 1.5], 1.0, 1.0);
        let (delta, r_eff, n) = hit.expect("should be in contact");
        assert!((delta - 0.5).abs() < 1e-9);
        assert!((r_eff - 0.5).abs() < 1e-9);
        assert!((n[2] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn overlap_none_when_separated() {
        assert!(capsule_overlap([0.0, 0.0, 0.0], [0.0, 0.0, 3.0], 1.0, 1.0).is_none());
    }
}
