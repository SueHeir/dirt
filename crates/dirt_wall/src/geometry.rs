//! Runtime wall geometry and state.

use crate::motion::WallMotion;
use soil_core::region::Region;

/// Runtime representation of an infinite plane wall.
///
/// Contact geometry for a plane wall:
/// ```text
///          ·  particle center (px, py, pz)
///         /|
///        / |  distance = dot(P - W, n̂)
///       /  |
///   n̂ ↑   |  delta (overlap) = radius - distance
///      |   |
///  ────W───┴──── wall surface
///      (point_x, point_y, point_z)
/// ```
///
/// The signed distance from the particle center to the wall plane is
/// `dot(pos - point, normal)`. Overlap is `delta = radius - distance`.
/// Force is applied whenever the particle overlaps the plane.  If an
/// integration step has carried a centre through the plane, the same inward
/// normal still restores it; a containment wall must not silently lose contact
/// on its forbidden side.
pub struct WallPlane {
    /// Current X-coordinate of a point on the plane (updated by motion).
    pub point_x: f64,
    /// Current Y-coordinate of a point on the plane (updated by motion).
    pub point_y: f64,
    /// Current Z-coordinate of a point on the plane (updated by motion).
    pub point_z: f64,
    /// X-component of the unit outward normal.
    pub normal_x: f64,
    /// Y-component of the unit outward normal.
    pub normal_y: f64,
    /// Z-component of the unit outward normal.
    pub normal_z: f64,
    /// Index into [`MaterialTable`] for this wall's material properties.
    pub material_index: usize,
    /// Optional name for runtime enable/disable.
    pub name: Option<String>,
    /// Lower X bound of the wall's active region.
    pub bound_x_low: f64,
    /// Upper X bound of the wall's active region.
    pub bound_x_high: f64,
    /// Lower Y bound of the wall's active region.
    pub bound_y_low: f64,
    /// Upper Y bound of the wall's active region.
    pub bound_y_high: f64,
    /// Lower Z bound of the wall's active region.
    pub bound_z_low: f64,
    /// Upper Z bound of the wall's active region.
    pub bound_z_high: f64,
    /// Current wall velocity `[vx, vy, vz]` (m/s).
    pub velocity: [f64; 3],
    /// Motion mode (static, constant velocity, oscillating, or servo).
    pub motion: WallMotion,
    /// Initial position of the wall point (reference for oscillation).
    pub origin: [f64; 3],
    /// Accumulated scalar contact force this timestep (along normal),
    /// used by the servo controller to compute the feedback error.
    pub force_accumulator: f64,
    /// Wall temperature in K (None = no wall heat transfer).
    pub temperature: Option<f64>,
}

impl WallPlane {
    /// Check if atom position is within the wall's bounding region.
    #[inline]
    pub fn in_bounds(&self, x: f64, y: f64, z: f64) -> bool {
        x >= self.bound_x_low
            && x <= self.bound_x_high
            && y >= self.bound_y_low
            && y <= self.bound_y_high
            && z >= self.bound_z_low
            && z <= self.bound_z_high
    }
}

/// Runtime representation of a cylindrical wall aligned to one axis.
///
/// Contact geometry (cross-section, `inside = true`):
/// ```text
///        ╭──────────╮  cylinder surface (radius R)
///       ╱            ╲
///      │   ← n̂  ·    │  radial_dist from axis to particle
///      │        ↑     │  gap = R - radial_dist
///       ╲      particle╱  delta = particle_radius - gap
///        ╰──────────╯
/// ```
///
/// For `inside = true`, the normal points inward (toward the axis).
/// For `inside = false`, the normal points outward (away from the axis).
/// Particles outside the `[lo, hi]` axial range are ignored.
pub struct WallCylinder {
    /// Axis index: 0 = X, 1 = Y, 2 = Z.
    pub axis: usize,
    /// Center of the cylinder in the 2D plane perpendicular to the axis.
    /// For a Z-cylinder, this is `[center_x, center_y]`.
    pub center: [f64; 2],
    /// Cylinder radius (meters).
    pub radius: f64,
    /// Lower axial bound — particles below this are ignored.
    pub lo: f64,
    /// Upper axial bound — particles above this are ignored.
    pub hi: f64,
    /// If `true`, particles live inside the cylinder and the contact normal
    /// points inward (toward the axis).
    pub inside: bool,
    /// Index into [`MaterialTable`] for this wall's material properties.
    pub material_index: usize,
    /// Optional name for runtime enable/disable.
    pub name: Option<String>,
    /// Accumulated scalar contact force this timestep.
    pub force_accumulator: f64,
    /// Wall temperature in K (None = no wall heat transfer).
    pub temperature: Option<f64>,
}

/// Runtime representation of a spherical wall.
///
/// Contact geometry (`inside = true`):
/// ```text
///         ╭─────╮  sphere surface (radius R)
///        ╱       ╲
///       │  ← n̂  · │  dist = |pos - center|
///       │       ↑  │  gap = R - dist
///        ╲  particle╱  delta = particle_radius - gap
///         ╰─────╯
/// ```
///
/// For `inside = true`, the normal points inward (toward the center).
/// For `inside = false`, the normal points outward (away from the center).
pub struct WallSphere {
    /// Sphere center `[x, y, z]` (meters).
    pub center: [f64; 3],
    /// Sphere radius (meters).
    pub radius: f64,
    /// If `true`, particles live inside the sphere and the contact normal
    /// points inward (toward the center).
    pub inside: bool,
    /// Index into [`MaterialTable`] for this wall's material properties.
    pub material_index: usize,
    /// Optional name for runtime enable/disable.
    pub name: Option<String>,
    /// Accumulated scalar contact force this timestep.
    pub force_accumulator: f64,
    /// Wall temperature in K (None = no wall heat transfer).
    pub temperature: Option<f64>,
}

/// Runtime representation of a region-surface wall.
///
/// Uses any [`Region`] shape (block, sphere, cylinder, cone, union,
/// intersect, etc.) as a wall surface. Contact detection delegates to
/// [`Region::closest_point_on_surface`] to get the signed distance and
/// outward normal for each particle.
pub struct WallRegion {
    /// The region whose surface acts as the wall.
    pub region: Region,
    /// If `true`, particles live inside the region and the contact normal
    /// points inward (away from the surface, toward the interior).
    pub inside: bool,
    /// Index into [`MaterialTable`] for this wall's material properties.
    pub material_index: usize,
    /// Optional name for runtime enable/disable.
    pub name: Option<String>,
    /// Accumulated scalar contact force this timestep.
    pub force_accumulator: f64,
    /// Wall temperature in K (None = no wall heat transfer).
    pub temperature: Option<f64>,
}

/// Collection of all wall types with per-wall active/inactive flags.
///
/// Stored as a resource in the [`App`]. Individual walls can be disabled or
/// enabled at runtime via [`deactivate_by_name`](Self::deactivate_by_name) and
/// [`activate_by_name`](Self::activate_by_name).
pub struct Walls {
    /// Plane walls.
    pub planes: Vec<WallPlane>,
    /// Per-plane active flags (parallel to `planes`).
    pub active: Vec<bool>,
    /// Cylinder walls.
    pub cylinders: Vec<WallCylinder>,
    /// Per-cylinder active flags (parallel to `cylinders`).
    pub cylinder_active: Vec<bool>,
    /// Sphere walls.
    pub spheres: Vec<WallSphere>,
    /// Per-sphere active flags (parallel to `spheres`).
    pub sphere_active: Vec<bool>,
    /// Region-surface walls.
    pub regions: Vec<WallRegion>,
    /// Per-region active flags (parallel to `regions`).
    pub region_active: Vec<bool>,
    /// Elapsed simulation time (seconds), used for oscillation phase tracking.
    pub time: f64,
    /// Per-contact Mindlin tangential-spring history for wall friction, keyed by
    /// `(wall_kind, wall_index, particle_tag)` where wall_kind is
    /// 0=plane, 1=cylinder, 2=sphere, 3=region. Rebuilt each step so contacts
    /// that end are pruned automatically.
    pub tangential_springs: std::collections::HashMap<(u8, usize, u32), [f64; 3]>,
    /// Per-contact rolling-displacement history for the SDS rolling-resistance
    /// model (same key scheme as `tangential_springs`). Empty under the default
    /// `constant` rolling model, which is stateless.
    pub rolling_springs: std::collections::HashMap<(u8, usize, u32), [f64; 3]>,
}

impl Walls {
    /// Deactivate all walls (of any type) whose `name` matches the given string.
    ///
    /// Deactivated walls are skipped during contact force computation and
    /// motion updates. This is useful for removing walls mid-simulation
    /// (e.g., removing a lid after compaction).
    pub fn deactivate_by_name(&mut self, name: &str) {
        for (i, wall) in self.planes.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.active[i] = false;
            }
        }
        for (i, wall) in self.cylinders.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.cylinder_active[i] = false;
            }
        }
        for (i, wall) in self.spheres.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.sphere_active[i] = false;
            }
        }
        for (i, wall) in self.regions.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.region_active[i] = false;
            }
        }
    }

    /// Activate all walls (of any type) whose `name` matches the given string.
    ///
    /// Reactivated walls participate in contact force computation and, for
    /// moving plane walls, motion updates again. Unknown names are ignored,
    /// matching [`deactivate_by_name`](Self::deactivate_by_name).
    pub fn activate_by_name(&mut self, name: &str) {
        for (i, wall) in self.planes.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.active[i] = true;
            }
        }
        for (i, wall) in self.cylinders.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.cylinder_active[i] = true;
            }
        }
        for (i, wall) in self.spheres.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.sphere_active[i] = true;
            }
        }
        for (i, wall) in self.regions.iter().enumerate() {
            if wall.name.as_deref() == Some(name) {
                self.region_active[i] = true;
            }
        }
    }
}
