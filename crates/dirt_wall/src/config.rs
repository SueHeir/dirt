//! Declarative wall configuration and validation helpers.

use dirt_atom::MaterialTable;
use serde::Deserialize;
use soil_core::region::Region;

fn default_neg_inf() -> f64 {
    f64::NEG_INFINITY
}
fn default_pos_inf() -> f64 {
    f64::INFINITY
}
fn default_wall_type() -> String {
    "plane".to_string()
}

pub(crate) fn curved_or_region_wall_surface_energy_warning(
    wall: &WallDef,
    material_table: &MaterialTable,
    material_index: usize,
) -> Option<String> {
    let geometry = match wall.wall_type.as_str() {
        "cylinder" => "cylinder",
        "sphere" => "sphere",
        "region" => "region",
        _ => return None,
    };
    let surface_energy = *material_table.surface_energy.get(material_index)?;
    if surface_energy <= 0.0 {
        return None;
    }

    let wall_name = wall
        .name
        .as_ref()
        .map(|name| format!(" named '{name}'"))
        .unwrap_or_default();
    Some(format!(
        "WARNING: {geometry} wall{wall_name} uses material '{}' with \
         surface_energy = {surface_energy}, but JKR/DMT `surface_energy` is \
         plane-wall-only and is ignored by cylinder, sphere, and region walls. \
         Use a plane wall for JKR/DMT wall adhesion, or use `cohesion_energy` \
         for unchanged SJKR cohesion on curved/region walls.",
        wall.material
    ))
}

// ── Config structs ──────────────────────────────────────────────────────────

/// Sinusoidal oscillation parameters for a wall.
///
/// The wall displaces along its normal as `amplitude * sin(2π * frequency * t)`.
/// Velocity is computed analytically as the time derivative.
///
/// # TOML
/// ```toml
/// oscillate = { amplitude = 0.001, frequency = 50.0 }
/// ```
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct OscillateDef {
    /// Peak displacement from the origin position (meters).
    pub amplitude: f64,
    /// Oscillation frequency (Hz).
    pub frequency: f64,
}

/// Proportional servo controller parameters for a wall.
///
/// Each timestep the servo computes `error = target_force - measured_force`,
/// then sets `velocity = clamp(gain * error, -max_velocity, max_velocity)`
/// along the wall normal. This drives the wall toward a steady-state contact
/// force equal to `target_force`.
///
/// # TOML
/// ```toml
/// servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }
/// ```
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct ServoDef {
    /// Desired total contact force on this wall (N).
    pub target_force: f64,
    /// Maximum wall velocity magnitude (m/s), prevents overshooting.
    pub max_velocity: f64,
    /// Proportional gain (m/s per N of force error).
    pub gain: f64,
}

/// TOML definition of a single wall entry (`[[wall]]`).
///
/// This is a union struct: different fields are relevant depending on `type`.
/// See the [module-level docs](self) for TOML examples of each wall type.
#[derive(Deserialize, Clone)]
pub struct WallDef {
    /// Wall type: `"plane"` (default), `"cylinder"`, `"sphere"`, or `"region"`.
    #[serde(default = "default_wall_type", rename = "type")]
    pub wall_type: String,

    // ── Plane fields ────────────────────────────────────────────────────
    /// X-coordinate of a point on the plane (default: 0.0).
    #[serde(default)]
    pub point_x: f64,
    /// Y-coordinate of a point on the plane (default: 0.0).
    #[serde(default)]
    pub point_y: f64,
    /// Z-coordinate of a point on the plane (default: 0.0).
    #[serde(default)]
    pub point_z: f64,
    /// X-component of the outward normal vector (will be normalized).
    #[serde(default)]
    pub normal_x: f64,
    /// Y-component of the outward normal vector (will be normalized).
    #[serde(default)]
    pub normal_y: f64,
    /// Z-component of the outward normal vector (will be normalized).
    #[serde(default)]
    pub normal_z: f64,

    // ── Cylinder/sphere fields ──────────────────────────────────────────
    /// Cylinder axis: `"x"`, `"y"`, or `"z"` (default: `"z"`).
    #[serde(default)]
    pub axis: Option<String>,
    /// Center coordinates: `[c0, c1]` for cylinder (in the plane ⊥ to axis),
    /// or `[x, y, z]` for sphere.
    #[serde(default)]
    pub center: Option<Vec<f64>>,
    /// Radius of the cylinder or sphere wall surface (meters).
    #[serde(default)]
    pub radius: Option<f64>,
    /// Lower axial bound for cylinder (default: −∞).
    #[serde(default)]
    pub lo: Option<f64>,
    /// Upper axial bound for cylinder (default: +∞).
    #[serde(default)]
    pub hi: Option<f64>,
    /// If `true`, particles live inside the wall surface and the contact
    /// normal points inward. If `false` (default), particles are outside.
    #[serde(default)]
    pub inside: Option<bool>,

    // ── Common fields ───────────────────────────────────────────────────
    /// Material name — must match a `[[dem.materials]]` entry.
    pub material: String,
    /// Optional name for runtime enable/disable via [`Walls::deactivate_by_name`]
    /// and [`Walls::activate_by_name`].
    #[serde(default)]
    pub name: Option<String>,
    /// Lower X bound for the plane wall's active region (default: −∞).
    #[serde(default = "default_neg_inf")]
    pub bound_x_low: f64,
    /// Upper X bound for the plane wall's active region (default: +∞).
    #[serde(default = "default_pos_inf")]
    pub bound_x_high: f64,
    /// Lower Y bound for the plane wall's active region (default: −∞).
    #[serde(default = "default_neg_inf")]
    pub bound_y_low: f64,
    /// Upper Y bound for the plane wall's active region (default: +∞).
    #[serde(default = "default_pos_inf")]
    pub bound_y_high: f64,
    /// Lower Z bound for the plane wall's active region (default: −∞).
    #[serde(default = "default_neg_inf")]
    pub bound_z_low: f64,
    /// Upper Z bound for the plane wall's active region (default: +∞).
    #[serde(default = "default_pos_inf")]
    pub bound_z_high: f64,
    /// Constant wall velocity `[vx, vy, vz]` (m/s). Plane walls only.
    #[serde(default)]
    pub velocity: Option<[f64; 3]>,
    /// Sinusoidal oscillation parameters. Plane walls only.
    #[serde(default)]
    pub oscillate: Option<OscillateDef>,
    /// Servo controller parameters. Plane walls only.
    #[serde(default)]
    pub servo: Option<ServoDef>,
    /// Region definition for `type = "region"` walls.
    #[serde(default)]
    pub region: Option<Region>,
    /// Wall temperature in K (None = no wall heat transfer).
    #[serde(default)]
    pub temperature: Option<f64>,
}
