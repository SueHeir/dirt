use super::*;

pub(super) fn default_source() -> String {
    "random".to_string()
}

/// Column index mapping for CSV file-based insertion.
#[derive(Deserialize, Clone, Debug)]
#[serde(deny_unknown_fields)]
pub struct ColumnMapping {
    /// Zero-based column index of the x coordinate.
    #[serde(default)]
    pub x: Option<usize>,
    /// Zero-based column index of the y coordinate.
    #[serde(default)]
    pub y: Option<usize>,
    /// Zero-based column index of the z coordinate.
    #[serde(default)]
    pub z: Option<usize>,
    /// Zero-based column index of the particle radius.
    #[serde(default)]
    pub radius: Option<usize>,
    /// Zero-based column index of the x velocity component.
    #[serde(default)]
    pub vx: Option<usize>,
    /// Zero-based column index of the y velocity component.
    #[serde(default)]
    pub vy: Option<usize>,
    /// Zero-based column index of the z velocity component.
    #[serde(default)]
    pub vz: Option<usize>,
    /// Zero-based column index of the integer atom type.
    #[serde(default)]
    pub atom_type: Option<usize>,
}

impl Default for ColumnMapping {
    fn default() -> Self {
        ColumnMapping {
            x: Some(0),
            y: Some(1),
            z: Some(2),
            radius: Some(3),
            vx: None,
            vy: None,
            vz: None,
            atom_type: None,
        }
    }
}

#[derive(Deserialize, Clone)]
/// A single particle insertion block from `[[particles.insert]]`.
///
/// Three modes determined by config fields:
/// - **Random** (default): `source = "random"` (or omitted), requires `material`, `count`, `radius`, `density`
/// - **Rate-based**: random insertion with `rate` field present — registers for periodic insertion
/// - **File-based**: `source = "file"`, requires `file` and `format`
pub struct InsertConfig {
    /// Insertion source: `"random"` (default) or `"file"`.
    #[serde(default = "default_source")]
    pub source: String,
    /// Material name (must match a `[[dem.materials]]` entry). Required for random/rate modes.
    pub material: Option<String>,
    /// Number of particles to insert at setup time. Required for random mode without `rate`.
    pub count: Option<u32>,
    /// Particle radius: fixed value or distribution. Required for random/rate modes.
    pub radius: Option<RadiusSpec>,
    /// Particle density (kg/m³). Required for random/rate modes.
    pub density: Option<f64>,
    /// Random velocity magnitude (Gaussian distribution).
    pub velocity: Option<f64>,
    /// Directional velocity components (additive with random velocity).
    pub velocity_x: Option<f64>,
    /// Directional y velocity component (m/s), additive with random velocity.
    pub velocity_y: Option<f64>,
    /// Directional z velocity component (m/s), additive with random velocity.
    pub velocity_z: Option<f64>,
    /// Insertion region. Defaults to domain bounds (inset by particle radius).
    #[serde(default)]
    pub region: Option<Region>,
    // ── Rate-based insertion fields ──
    /// Particles to insert per interval. Presence of this field triggers rate mode.
    pub rate: Option<u32>,
    /// Insert every N timesteps (default: 1).
    pub rate_interval: Option<usize>,
    /// First timestep to begin insertion (default: 0).
    pub rate_start: Option<usize>,
    /// Last timestep for insertion (optional, no default = run forever).
    pub rate_end: Option<usize>,
    /// Maximum total particles to insert (optional).
    pub rate_limit: Option<u32>,
    // ── File-based insertion fields ──
    /// Path to particle data file.
    pub file: Option<String>,
    /// File format: `"csv"` or `"lammps_dump"`.
    pub format: Option<String>,
    /// Column index mapping for CSV files.
    pub columns: Option<ColumnMapping>,
    /// Explicit mapping from integer atom types (in file) to named materials.
    /// e.g. `type_map = { 1 = "glass", 2 = "steel" }`.
    /// When present, overrides implicit type-to-material mapping.
    /// Keys are strings in TOML but parsed as u32 integers.
    #[serde(default)]
    pub type_map: Option<HashMap<String, String>>,
    /// LAMMPS data file atom style: `"atomic"`, `"sphere"`, `"bpm/sphere"`.
    /// Auto-detected from `Atoms # style` header if not specified.
    pub atom_style: Option<String>,
    /// Seed for the deterministic insertion RNG. Defaults to 0.
    ///
    /// Insertion runs on every MPI rank with the SAME seed so each rank
    /// generates the identical candidate stream (positions, radii, velocities,
    /// tags) and therefore the identical global packing; a rank only stores the
    /// atoms whose position falls inside its own subdomain. A fixed default
    /// makes packings reproducible across runs and across rank counts.
    #[serde(default)]
    pub seed: Option<u64>,
}

/// TOML `[particles]` — contains a list of insertion blocks.
#[derive(Deserialize, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct ParticlesConfig {
    /// The ordered list of `[[particles.insert]]` blocks to process.
    pub insert: Option<Vec<InsertConfig>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum InsertFileError {
    MissingField {
        source: &'static str,
        field: &'static str,
    },
    UnknownFormat {
        format: String,
    },
    InvalidTypeMapKey {
        key: String,
    },
    FileOpen {
        path: String,
        source: String,
    },
    FileRead {
        path: String,
        line: usize,
        source: String,
    },
    MissingSection {
        path: String,
        section: &'static str,
    },
    MissingColumn {
        path: String,
        line: usize,
        field: &'static str,
    },
    MissingDefault {
        path: String,
        field: &'static str,
        context: &'static str,
    },
    ParseField {
        path: String,
        line: usize,
        field: String,
        value: String,
        source: String,
    },
    RowTooShort {
        path: String,
        line: usize,
        style: String,
        expected: usize,
        found: usize,
    },
    UnsupportedAtomStyle {
        path: String,
        style: String,
    },
}

impl fmt::Display for InsertFileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InsertFileError::MissingField { source, field } => write!(
                f,
                "{} source = \"file\" requires '{}' in [[particles.insert]]",
                source, field
            ),
            InsertFileError::UnknownFormat { format } => write!(
                f,
                "Unknown file format '{}' in [[particles.insert]]. Supported: csv, lammps_dump, lammps_data",
                format
            ),
            InsertFileError::InvalidTypeMapKey { key } => {
                write!(f, "type_map key '{}' is not a valid integer atom type", key)
            }
            InsertFileError::FileOpen { path, source } => {
                write!(f, "failed to open particle file '{}': {}", path, source)
            }
            InsertFileError::FileRead { path, line, source } => {
                write!(f, "failed to read line {} of '{}': {}", line, path, source)
            }
            InsertFileError::MissingSection { path, section } => {
                write!(f, "no '{}' section found in LAMMPS data file '{}'", section, path)
            }
            InsertFileError::MissingColumn { path, line, field } => write!(
                f,
                "missing or invalid '{}' column at line {} of '{}'",
                field, line, path
            ),
            InsertFileError::MissingDefault {
                path,
                field,
                context,
            } => write!(
                f,
                "'{}' required in config for {} while reading '{}'",
                field, context, path
            ),
            InsertFileError::ParseField {
                path,
                line,
                field,
                value,
                source,
            } => write!(
                f,
                "failed to parse {} {:?} at line {} of '{}': {}",
                field, value, line, path, source
            ),
            InsertFileError::RowTooShort {
                path,
                line,
                style,
                expected,
                found,
            } => write!(
                f,
                "expected at least {} columns for {} style at line {} of '{}', found {}",
                expected, style, line, path, found
            ),
            InsertFileError::UnsupportedAtomStyle { path, style } => write!(
                f,
                "unsupported atom_style '{}' in LAMMPS data file '{}'. Supported: atomic, sphere, bpm/sphere",
                style, path
            ),
        }
    }
}

impl std::error::Error for InsertFileError {}

// ── Rate-based insertion state ──────────────────────────────────────────────

/// Tracks a single rate-based insertion configuration and its progress.
pub struct RateInsertEntry {
    /// The originating insertion configuration for this rate stream.
    pub config: InsertConfig,
    /// Validated, domain-aware insertion inputs shared with immediate insertion.
    pub(super) prepared: PreparedInsert,
    /// Running count of particles inserted so far by this entry.
    pub total_inserted: u32,
}

/// Resource holding all active rate-based insertion entries.
pub struct RateInsertState {
    /// One entry per rate-based `[[particles.insert]]` block still active.
    pub entries: Vec<RateInsertEntry>,
}

impl Default for RateInsertState {
    fn default() -> Self {
        RateInsertState {
            entries: Vec::new(),
        }
    }
}

/// Inputs common to random immediate and rate insertion after validation.
///
/// Keeping this private prevents a rate entry from retaining an unvalidated
/// configuration and makes both insertion policies use exactly the same
/// material lookup, region construction, radius bounds, and velocity setup.
#[derive(Clone, Debug)]
pub(super) struct PreparedInsert {
    pub(super) mat_idx: u32,
    pub(super) density: f64,
    pub(super) radius: RadiusSpec,
    pub(super) max_radius: f64,
    pub(super) region: Region,
    pub(super) cutoff_padding: f64,
    pub(super) velocity_base: [f64; 3],
    pub(super) velocity_normal: Option<Normal<f64>>,
    pub(super) seed: u64,
}

pub(super) fn prepare_random_insert(
    insert: &InsertConfig,
    materials: &MaterialTable,
    domain: &Domain,
    context: &str,
) -> Result<PreparedInsert, String> {
    let material = insert
        .material
        .as_deref()
        .ok_or_else(|| format!("{context} requires 'material'"))?;
    let mat_idx = resolve_material(materials, material)?;
    let radius = insert
        .radius
        .as_ref()
        .ok_or_else(|| format!("{context} requires 'radius'"))?
        .clone();
    let density = insert
        .density
        .ok_or_else(|| format!("{context} requires 'density'"))?;
    let max_radius = radius
        .try_max_radius()
        .map_err(|error| format!("invalid radius in {context}: {error}"))?;
    let random_velocity = insert.velocity.unwrap_or(0.0);
    validate_insert_velocity(random_velocity, context)?;
    let region = insert
        .region
        .clone()
        .unwrap_or(default_insert_region(domain, max_radius, context)?);
    validate_insert_region(&region, context)?;
    Ok(PreparedInsert {
        mat_idx,
        density,
        radius,
        max_radius,
        region,
        cutoff_padding: materials.liquid_bridge_cutoff_padding(mat_idx),
        velocity_base: [
            insert.velocity_x.unwrap_or(0.0),
            insert.velocity_y.unwrap_or(0.0),
            insert.velocity_z.unwrap_or(0.0),
        ],
        velocity_normal: (random_velocity > 0.0).then(|| {
            Normal::new(0.0, random_velocity)
                .expect("validated insertion velocity must create a Normal distribution")
        }),
        seed: insert.seed.unwrap_or(0),
    })
}

pub(super) fn validate_insert_velocity(rand_vel: f64, context: &str) -> Result<(), String> {
    if !rand_vel.is_finite() || rand_vel < 0.0 {
        return Err(format!(
            "velocity in {context} must be finite and non-negative, got {rand_vel}"
        ));
    }
    Ok(())
}

pub(super) fn is_rate_insert_config(insert: &InsertConfig) -> bool {
    insert.rate.is_some()
        || insert.rate_interval.is_some()
        || insert.rate_start.is_some()
        || insert.rate_end.is_some()
        || insert.rate_limit.is_some()
}

pub(super) fn validate_rate_insert_config<'a>(
    insert: &'a InsertConfig,
    context: &str,
) -> Result<(u32, &'a RadiusSpec, f64), String> {
    let rate = insert
        .rate
        .ok_or_else(|| format!("{context} requires 'rate' for rate-based insertion"))?;
    let radius = insert
        .radius
        .as_ref()
        .ok_or_else(|| format!("{context} requires 'radius' for rate-based insertion"))?;
    let density = insert
        .density
        .ok_or_else(|| format!("{context} requires 'density' for rate-based insertion"))?;
    Ok((rate, radius, density))
}

pub(super) fn default_insert_region(
    domain: &Domain,
    max_r: f64,
    context: &str,
) -> Result<Region, String> {
    let min = [
        domain.boundaries_low[0] + max_r,
        domain.boundaries_low[1] + max_r,
        domain.boundaries_low[2] + max_r,
    ];
    let max = [
        domain.boundaries_high[0] - max_r,
        domain.boundaries_high[1] - max_r,
        domain.boundaries_high[2] - max_r,
    ];

    for axis in 0..3 {
        if min[axis] > max[axis] {
            let extent = domain.boundaries_high[axis] - domain.boundaries_low[axis];
            return Err(format!(
                "{context} default insertion region is smaller than particle: radius {max_r} exceeds domain extent {extent} on axis {axis}"
            ));
        }
    }

    Ok(Region::Block { min, max })
}

pub(super) fn validate_insert_region(region: &Region, context: &str) -> Result<(), String> {
    match region {
        Region::Block { min, max } => {
            for axis in 0..3 {
                if min[axis] >= max[axis] {
                    return Err(format!(
                        "{context} insertion region is empty or degenerate on axis {axis}: min {} must be less than max {}",
                        min[axis], max[axis]
                    ));
                }
            }
            Ok(())
        }
        Region::Sphere { radius, .. } => {
            if !radius.is_finite() || *radius <= 0.0 {
                return Err(format!(
                    "{context} sphere insertion region radius must be finite and > 0, got {radius}"
                ));
            }
            Ok(())
        }
        Region::Cylinder { radius, lo, hi, .. } => {
            if !radius.is_finite() || *radius <= 0.0 {
                return Err(format!(
                    "{context} cylinder insertion region radius must be finite and > 0, got {radius}"
                ));
            }
            if lo >= hi {
                return Err(format!(
                    "{context} cylinder insertion region is empty or degenerate: lo {lo} must be less than hi {hi}"
                ));
            }
            Ok(())
        }
        Region::Cone {
            rad_lo,
            rad_hi,
            lo,
            hi,
            ..
        } => {
            if !rad_lo.is_finite() || *rad_lo < 0.0 {
                return Err(format!(
                    "{context} cone insertion region rad_lo must be finite and >= 0, got {rad_lo}"
                ));
            }
            if !rad_hi.is_finite() || *rad_hi < 0.0 {
                return Err(format!(
                    "{context} cone insertion region rad_hi must be finite and >= 0, got {rad_hi}"
                ));
            }
            if *rad_lo == 0.0 && *rad_hi == 0.0 {
                return Err(format!(
                    "{context} cone insertion region is empty or degenerate: at least one end radius must be > 0"
                ));
            }
            if lo >= hi {
                return Err(format!(
                    "{context} cone insertion region is empty or degenerate: lo {lo} must be less than hi {hi}"
                ));
            }
            Ok(())
        }
        Region::Plane { .. } => Err(format!(
            "{context} plane insertion region is unbounded and cannot be sampled"
        )),
        Region::Union { regions } => {
            if regions.is_empty() {
                return Err(format!("{context} union insertion region is empty"));
            }
            for child in regions {
                validate_insert_region(child, context)?;
            }
            Ok(())
        }
        Region::Intersect { regions } => {
            if regions.is_empty() {
                return Err(format!("{context} intersect insertion region is empty"));
            }
            for child in regions {
                validate_insert_region(child, context)?;
            }
            Ok(())
        }
    }
}
