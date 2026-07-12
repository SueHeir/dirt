//! Particle insertion: random, rate-based, and file-based from `[[particles.insert]]` config.
//!
//! # Parallel insertion model (born-in-owner)
//!
//! Insertion is **not** done on one rank and scattered. Instead it runs on
//! *every* MPI rank, and every rank seeds its RNG with the **same** `seed`, so
//! all ranks generate the bit-identical candidate stream — the same positions,
//! radii, velocities, and tags, in the same order. Each rank then keeps **only**
//! the candidates whose position falls inside its own subdomain, tested with a
//! half-open interval (`low ≤ pos < high`, via the `owns_position` helper) that matches
//! `Domain::exchange()` ownership exactly.
//!
//! Three properties follow:
//!
//! - **Born in its owner.** Each atom is materialized only on the rank that owns
//!   its position, so the first post-insertion `exchange()` never has to migrate
//!   it — no insertion-time communication is needed.
//! - **Exactly-once.** The half-open convention guarantees every position is
//!   claimed by one rank, never two and never zero (no duplicated or dropped
//!   atoms at subdomain boundaries).
//! - **Deterministic across rank counts.** Because the candidate stream depends
//!   only on `seed` (not on the decomposition), the *global* packing is
//!   identical whether you run on 1 rank or 64 — only the partitioning of that
//!   packing changes. Overlap rejection uses a global spatial hash replicated on
//!   every rank, so packings are reproducible run-to-run and rank-count to
//!   rank-count.
//!
//! File-based insertion follows the same rule: the file is parsed on every rank
//! and filtered to the local subdomain, with tags advanced identically so they
//! stay globally consistent.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use serde::Deserialize;

use grass_scheduler::prelude::CurrentState;
use soil_core::{
    deep_merge, Atom, AtomDataRegistry, CommResource, CommState, Config, Domain, DomainConfig,
    ParticleSimScheduleSet, ParticleStore, Real, Region, RunConfig, RunState, ScheduleSetupSet,
    StageOverrides,
};

use crate::{DemAtom, MaterialTable, RadiusSpec};

// ── Particle insertion ─────────────────────────────────────────────────────

fn default_source() -> String {
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
enum InsertFileError {
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
    /// Resolved material-table index for the inserted particles.
    pub mat_idx: u32,
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

fn validate_insert_velocity(rand_vel: f64, context: &str) -> Result<(), String> {
    if !rand_vel.is_finite() || rand_vel < 0.0 {
        return Err(format!(
            "velocity in {context} must be finite and non-negative, got {rand_vel}"
        ));
    }
    Ok(())
}

fn is_rate_insert_config(insert: &InsertConfig) -> bool {
    insert.rate.is_some()
        || insert.rate_interval.is_some()
        || insert.rate_start.is_some()
        || insert.rate_end.is_some()
        || insert.rate_limit.is_some()
}

fn validate_rate_insert_config<'a>(
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

fn default_insert_region(domain: &Domain, max_r: f64, context: &str) -> Result<Region, String> {
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

fn validate_insert_region(region: &Region, context: &str) -> Result<(), String> {
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

/// Samples one candidate point, treating SOIL's bounded rejection-sampling
/// exhaustion as a rejected insertion attempt.
///
/// A composite region can be structurally valid while its overlap has such a
/// small measure that SOIL exhausts its finite inner rejection budget.  No
/// particle can be constructed in that case.  The outer insertion loop already
/// has a deterministic attempt limit, so retrying there is the same physical
/// outcome as an overlap rejection: successful draws keep the pre-existing RNG
/// stream and MPI ownership rules, while an exhausted draw cannot panic or
/// terminate just one simulation rank.
fn try_sample_insertion_point(region: &Region, rng: &mut impl rand::Rng) -> Option<[f64; 3]> {
    region.random_point_inside(rng).ok()
}

// ── SpatialHash for O(1) overlap checking ───────────────────────────────────

/// Grid-based spatial hash for fast overlap detection during particle insertion.
///
/// Divides space into cubic cells of size `cell_size` (typically ~2× max particle diameter).
/// Overlap queries check the 3×3×3 neighborhood of the candidate cell, ensuring all
/// potential overlaps are found without a full O(N²) scan.
/// Periodic boundary info for minimum-image overlap checks during insertion.
struct PeriodicBox {
    is_periodic: [bool; 3],
    box_size: [f64; 3],
}

impl PeriodicBox {
    fn from_domain(domain: &Domain) -> Self {
        PeriodicBox {
            is_periodic: domain.periodic_flags(),
            box_size: domain.size,
        }
    }

    /// Compute minimum-image squared distance between two positions.
    fn min_image_dist_sq(&self, a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let mut dist_sq = 0.0;
        for d in 0..3 {
            let mut delta = a[d] - b[d];
            if self.is_periodic[d] {
                let half = 0.5 * self.box_size[d];
                if delta > half {
                    delta -= self.box_size[d];
                } else if delta < -half {
                    delta += self.box_size[d];
                }
            }
            dist_sq += delta * delta;
        }
        dist_sq
    }
}

struct SpatialHash {
    cell_size: f64,
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl SpatialHash {
    fn new(cell_size: f64) -> Self {
        SpatialHash {
            cell_size,
            cells: HashMap::new(),
        }
    }

    fn cell_key(&self, pos: &[f64; 3]) -> (i64, i64, i64) {
        (
            (pos[0] / self.cell_size).floor() as i64,
            (pos[1] / self.cell_size).floor() as i64,
            (pos[2] / self.cell_size).floor() as i64,
        )
    }

    fn insert(&mut self, idx: usize, pos: &[f64; 3]) {
        let key = self.cell_key(pos);
        self.cells.entry(key).or_default().push(idx);
    }

    fn has_overlap(
        &self,
        pos: &[f64; 3],
        radius: f64,
        positions: &[[f64; 3]],
        radii: &[f64],
        pbc: &PeriodicBox,
    ) -> bool {
        // Collect all cell keys to check: 3x3x3 neighborhood plus periodic images
        let key = self.cell_key(pos);
        let min_dist_check = radius * 2.2; // conservative check radius

        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let neighbor_key = (key.0 + di, key.1 + dj, key.2 + dk);
                    if let Some(indices) = self.cells.get(&neighbor_key) {
                        for &idx in indices {
                            let dist_sq = pbc.min_image_dist_sq(pos, &positions[idx]);
                            let min_dist = (radius + radii[idx]) * 1.1;
                            if dist_sq <= min_dist * min_dist {
                                return true;
                            }
                        }
                    }
                }
            }
        }

        // For periodic axes where the box is small (< 3 cell sizes), the standard
        // 3x3x3 neighborhood may miss periodic images. Do a brute-force check
        // against all atoms using minimum-image distances.
        let needs_pbc_check = (0..3)
            .any(|d| pbc.is_periodic[d] && pbc.box_size[d] < 3.0 * self.cell_size + min_dist_check);
        if needs_pbc_check {
            for idx in 0..positions.len() {
                let dist_sq = pbc.min_image_dist_sq(pos, &positions[idx]);
                let min_dist = (radius + radii[idx]) * 1.1;
                if dist_sq <= min_dist * min_dist {
                    return true;
                }
            }
        }

        false
    }
}

// ── DemAtomInsertPlugin ─────────────────────────────────────────────────────

/// Inserts DEM particles at setup time and registers rate-based insertion for runtime.
pub struct DemAtomInsertPlugin;

impl Plugin for DemAtomInsertPlugin {
    fn default_config(&self) -> Option<&str> {
        Some(
            r#"# Particle insertion blocks (one per material/group)
[[particles.insert]]
material = "glass"          # must match a [[dem.materials]] name
count = 100
radius = 0.001
density = 2500.0
# velocity = 0.1            # random velocity magnitude (Gaussian)
# velocity_x = 0.0          # directional velocity (additive with random)
# velocity_y = 0.0
# velocity_z = 0.0
# region = { type = "block", min = [0.0, 0.0, 0.0], max = [1.0, 1.0, 1.0] }  # defaults to domain bounds
#
# Size distributions (instead of fixed radius):
# radius = { distribution = "uniform", min = 0.0008, max = 0.0012 }
# radius = { distribution = "gaussian", mean = 0.001, std = 0.0001 }
# radius = { distribution = "lognormal", mean = 0.001, std = 0.0001 }
# radius = { distribution = "discrete", values = [0.001, 0.0015], weights = [0.7, 0.3] }
#
# Rate-based insertion (insert particles over time):
# rate = 10              # particles per interval
# rate_interval = 100    # insert every N timesteps
# rate_start = 0         # first timestep (default 0)
# rate_end = 500000      # last timestep (optional)
# rate_limit = 5000      # total max particles (optional)
#
# File-based insertion:
# source = "file"
# file = "particles.csv"
# format = "csv"
# material = "glass"
# density = 2500.0
# columns = { x = 0, y = 1, z = 2, radius = 3 }"#,
        )
    }

    fn build(&self, app: &mut App) {
        app.add_resource(RateInsertState::default());
        // Insertion must run AFTER domain decomposition so each rank's subdomain
        // bounds (sub_domain_low/high) are populated before atoms are placed —
        // parallel insertion filters candidates against those bounds.
        app.add_setup_system(
            dem_insert_atoms.after("domain_read_input"),
            ScheduleSetupSet::Setup,
        )
        .add_setup_system(calculate_delta_time, ScheduleSetupSet::PostSetup)
        .add_update_system(
            dem_rate_insert,
            ParticleSimScheduleSet::PreInitialIntegration,
        );
    }

    fn try_build(&self, app: &mut App) -> Result<(), AppError> {
        let config = Config::try_load::<ParticlesConfig>(app, "particles")
            .map_err(|error| AppError::message(error.to_string()))?;
        let domain_config = Config::try_load::<DomainConfig>(app, "domain")
            .map_err(|error| AppError::message(error.to_string()))?;
        let materials = app
            .get_resource_ref::<MaterialTable>()
            .ok_or_else(|| AppError::message("DemAtomInsertPlugin requires DemAtomPlugin"))?;
        validate_particles_config(&config, &materials, &domain_config)
            .map_err(AppError::message)?;
        validate_stage_particles_configs(app, &materials, &domain_config)
            .map_err(AppError::message)?;
        drop(materials);
        self.build(app);
        Ok(())
    }
}

/// Validates insertion inputs before scheduling, so malformed input is returned
/// by `App::try_add_plugins` rather than reaching a simulation system.
fn validate_particles_config(
    config: &ParticlesConfig,
    materials: &MaterialTable,
    domain_config: &DomainConfig,
) -> Result<(), String> {
    for insert in config.insert.as_deref().unwrap_or(&[]) {
        match insert.source.as_str() {
            "file" => {
                let name = insert
                    .material
                    .as_deref()
                    .ok_or("file [[particles.insert]] requires 'material'")?;
                resolve_material(materials, name)?;
                insert
                    .file
                    .as_deref()
                    .ok_or("file [[particles.insert]] requires 'file'")?;
                insert
                    .format
                    .as_deref()
                    .ok_or("file [[particles.insert]] requires 'format'")?;
                if let Some(map) = &insert.type_map {
                    resolve_type_map(map, materials).map_err(|e| e.to_string())?;
                }
                // Parse the complete file during fallible plugin assembly.  Merely
                // checking its path leaves malformed rows to the legacy setup
                // scheduler, where they cannot be returned to the runner.
                validate_file_insert(insert, materials).map_err(|e| e.to_string())?;
            }
            "random" => {
                let name = insert
                    .material
                    .as_deref()
                    .ok_or("[[particles.insert]] requires 'material'")?;
                resolve_material(materials, name)?;
                if is_rate_insert_config(insert) {
                    let (_, radius, _) =
                        validate_rate_insert_config(insert, "rate-based [[particles.insert]]")?;
                    let max_radius = radius.try_max_radius().map_err(|error| {
                        format!("invalid radius in rate-based [[particles.insert]]: {error}")
                    })?;
                    validate_insert_velocity(
                        insert.velocity.unwrap_or(0.0),
                        "rate-based [[particles.insert]]",
                    )?;
                    validate_preflight_region(
                        insert,
                        domain_config,
                        max_radius,
                        "rate-based [[particles.insert]]",
                    )?;
                } else {
                    insert
                        .count
                        .ok_or("[[particles.insert]] requires 'count' for random insertion")?;
                    let radius = insert
                        .radius
                        .as_ref()
                        .ok_or("[[particles.insert]] requires 'radius' for random insertion")?;
                    insert
                        .density
                        .ok_or("[[particles.insert]] requires 'density' for random insertion")?;
                    let max_radius = radius.try_max_radius().map_err(|error| {
                        format!("invalid radius in [[particles.insert]]: {error}")
                    })?;
                    validate_insert_velocity(
                        insert.velocity.unwrap_or(0.0),
                        "[[particles.insert]]",
                    )?;
                    validate_preflight_region(
                        insert,
                        domain_config,
                        max_radius,
                        "[[particles.insert]]",
                    )?;
                }
            }
            other => {
                return Err(format!(
                    "unknown particles.insert source '{other}'; supported: random, file"
                ))
            }
        }
    }
    Ok(())
}

/// Preflights every stage-level `[particles]` override before insertion systems
/// are registered.  `dem_insert_atoms` reads a deep-merged stage table at setup
/// time, so only checking the top-level section would leave later-stage errors
/// to the infallible scheduler path.
fn validate_stage_particles_configs(
    app: &App,
    materials: &MaterialTable,
    domain_config: &DomainConfig,
) -> Result<(), String> {
    let config_table = app
        .get_resource_ref::<Config>()
        .map(|config| config.table.clone())
        .unwrap_or_default();
    let run_config = app
        .get_resource_ref::<RunConfig>()
        .ok_or("DemAtomInsertPlugin requires RunPlugin")?;

    for (stage_index, stage) in run_config.stages.iter().enumerate() {
        if !stage.overrides.contains_key("particles") {
            continue;
        }

        let mut merged = config_table.clone();
        deep_merge(&mut merged, &stage.overrides);
        let particles = merged
            .get("particles")
            .cloned()
            .ok_or_else(|| format!("stage {} has an invalid [particles] override", stage_index))?
            .try_into::<ParticlesConfig>()
            .map_err(|error| {
                format!(
                    "failed to parse [particles] override for [[run]] stage {}: {}",
                    stage_index, error
                )
            })?;
        validate_particles_config(&particles, materials, domain_config).map_err(|error| {
            format!(
                "invalid [particles] override for [[run]] stage {}: {}",
                stage_index, error
            )
        })?;
    }
    Ok(())
}

fn validate_preflight_region(
    insert: &InsertConfig,
    domain_config: &DomainConfig,
    max_radius: f64,
    context: &str,
) -> Result<(), String> {
    if let Some(region) = &insert.region {
        return validate_insert_region(region, context);
    }
    let mut domain = Domain::default();
    domain.boundaries_low = [
        domain_config.x_low,
        domain_config.y_low,
        domain_config.z_low,
    ];
    domain.boundaries_high = [
        domain_config.x_high,
        domain_config.y_high,
        domain_config.z_high,
    ];
    default_insert_region(&domain, max_radius, context).map(|_| ())
}

/// Parses a file insertion into scratch stores during plugin preflight.
///
/// The real insertion still happens after domain decomposition, so ownership
/// filtering and tags retain their normal setup-time semantics.  The scratch
/// parse deliberately exercises the same readers, making file open/read/row
/// errors typed `AppError`s from `try_add_plugins`/`try_start`.
fn validate_file_insert(
    insert: &InsertConfig,
    materials: &MaterialTable,
) -> Result<(), InsertFileError> {
    let mut atom = Atom::default();
    let mut registry = AtomDataRegistry::new();
    registry
        .try_register(DemAtom::default(), 0)
        .expect("fresh registry accepts DemAtom");
    let mut max_tag = 0;
    insert_from_file(
        insert,
        &mut atom,
        &registry,
        materials,
        &Domain::default(),
        &mut max_tag,
    )
}

// ── Helper: insert a single particle ────────────────────────────────────────

/// Appends a single DEM particle to both the shared `Atom` arrays and the `DemAtom` extension.
///
/// Computes mass from density and radius (solid sphere: m = ρ·4/3·π·r³), and inverse
/// moment of inertia (I = 2/5·m·r² for a solid sphere). Initializes quaternion to identity
/// and angular velocity/momentum/torque to zero.
/// Typed initialization record for a newly-created DEM SoA row.  The store owns
/// structural synchronization; this record only supplies DIRT-specific values.
#[derive(Clone, Copy, Debug)]
struct DemParticle {
    pos: [f64; 3],
    vel: [f64; 3],
    radius: f64,
    cutoff_padding: f64,
    density: f64,
    mat_idx: u32,
    tag: u32,
}

impl DemParticle {
    fn mass(self) -> f64 {
        self.density * 4.0 / 3.0 * PI * self.radius.powi(3)
    }

    /// Populate the already-reserved core row.  Structural changes belong to
    /// `ParticleStore`; this method deliberately only writes an existing row.
    fn write_core(self, atom: &mut Atom, i: usize, mass: f64) {
        atom.tag[i] = self.tag;
        atom.origin_index[i] = 0;
        atom.cutoff_radius[i] = (self.radius + self.cutoff_padding.max(0.0)) as Real;
        atom.image[i] = [0, 0, 0];
        atom.is_ghost[i] = false;
        atom.pos[i] = [
            self.pos[0] as Real,
            self.pos[1] as Real,
            self.pos[2] as Real,
        ];
        atom.vel[i] = [
            self.vel[0] as Real,
            self.vel[1] as Real,
            self.vel[2] as Real,
        ];
        atom.force[i] = [0.0; 3];
        atom.mass[i] = mass as Real;
        atom.inv_mass[i] = (1.0 / mass) as Real;
        atom.atom_type[i] = self.mat_idx;
    }

    /// Populate the matching, already-reserved DEM extension row.
    fn write_dem(self, dem: &mut DemAtom, i: usize, mass: f64) {
        dem.radius[i] = self.radius;
        dem.density[i] = self.density;
        dem.inv_inertia[i] = 1.0 / (0.4 * mass * self.radius * self.radius);
        dem.quaternion[i] = [1.0, 0.0, 0.0, 0.0];
        dem.omega[i] = [0.0; 3];
        dem.ang_mom[i] = [0.0; 3];
        dem.torque[i] = [0.0; 3];
        dem.body_id[i] = 0.0;
    }
}

fn insert_single_particle(atom: &mut Atom, registry: &AtomDataRegistry, row: DemParticle) {
    let global_natoms = atom
        .natoms
        .checked_add(1)
        .expect("global particle count overflow during DEM insertion");
    ParticleStore::new(atom, registry)
        .push_default_local(global_natoms)
        .expect("registered DEM rows must accept transactional insertion");
    let i = atom.len() - 1;
    let mass = row.mass();
    row.write_core(atom, i, mass);
    let mut dem_data = registry.expect_mut::<DemAtom>("insert_single_particle");
    row.write_dem(&mut dem_data, i, mass);
}

// ── Helper: subdomain ownership ─────────────────────────────────────────────

/// Whether `pos` lies inside this rank's subdomain, using a half-open interval
/// `[sub_domain_low, sub_domain_high)` per axis — consistent with how
/// `exchange()` defines ownership (an atom is sent low if `pos < low`, sent high
/// if `pos >= high`). The half-open convention guarantees every position is
/// claimed by exactly one rank, never two and never none.
fn owns_position(domain: &Domain, pos: &[f64; 3]) -> bool {
    (0..3).all(|d| pos[d] >= domain.sub_domain_low[d] && pos[d] < domain.sub_domain_high[d])
}

// ── Helper: resolve material index ──────────────────────────────────────────

fn resolve_material(material_table: &MaterialTable, name: &str) -> Result<u32, String> {
    material_table.find_material(name).ok_or_else(|| {
        format!(
            "unknown material '{}' in [[particles.insert]]. Available: {:?}",
            name, material_table.names
        )
    })
}

fn resolve_file_material(
    material_table: &MaterialTable,
    name: &str,
) -> Result<u32, InsertFileError> {
    resolve_material(material_table, name).map_err(|message| InsertFileError::ParseField {
        path: "config".to_string(),
        line: 0,
        field: "material".to_string(),
        value: name.to_string(),
        source: message,
    })
}

// ── Helper: resolve type_map to index map ────────────────────────────────────

/// Validates material names in `type_map` and builds a `HashMap<u32, u32>` mapping
/// file atom types to material indices. Called once per file load.
fn resolve_type_map(
    type_map: &HashMap<String, String>,
    material_table: &MaterialTable,
) -> Result<HashMap<u32, u32>, InsertFileError> {
    let mut index_map = HashMap::new();
    for (key_str, mat_name) in type_map {
        let file_type: u32 = key_str
            .parse()
            .map_err(|_| InsertFileError::InvalidTypeMapKey {
                key: key_str.clone(),
            })?;
        let mat_idx = resolve_material(material_table, mat_name).map_err(|message| {
            InsertFileError::ParseField {
                path: "config".to_string(),
                line: 0,
                field: "type_map material".to_string(),
                value: mat_name.clone(),
                source: message,
            }
        })?;
        index_map.insert(file_type, mat_idx);
    }
    Ok(index_map)
}

/// Look up material index for a given file atom type.
/// Checks type_map first, then falls back to the default material index.
fn lookup_material_for_type(
    file_type: u32,
    type_index_map: Option<&HashMap<u32, u32>>,
    default_mat_idx: u32,
) -> u32 {
    if let Some(map) = type_index_map {
        if let Some(&idx) = map.get(&file_type) {
            return idx;
        }
    }
    default_mat_idx
}

// ── Setup system: dem_insert_atoms ──────────────────────────────────────────

/// Setup system that processes all `[[particles.insert]]` blocks at simulation start.
///
/// For each block: immediate random insertion places particles with overlap checking,
/// file-based insertion loads from CSV/LAMMPS files, and rate-based insertion registers
/// entries in [`RateInsertState`] for periodic insertion during the run.
pub fn dem_insert_atoms(
    comm: Res<CommResource>,
    domain: Res<Domain>,
    mut atom: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    stage_overrides: Res<StageOverrides>,
    run_config: Res<RunConfig>,
    scheduler_manager: Res<SchedulerManager>,
    mut rate_state: ResMut<RateInsertState>,
) {
    let index = scheduler_manager.index;

    // Determine if this stage should insert particles:
    // - First stage: use top-level [particles] (backward compat) or stage overrides
    // - Later stages: only if the stage's [[run]] block explicitly has particles
    let has_stage_particles = index < run_config.num_stages()
        && run_config
            .current_stage(index)
            .overrides
            .contains_key("particles");

    let particles_config: ParticlesConfig = if has_stage_particles || index == 0 {
        stage_overrides.section("particles")
    } else {
        ParticlesConfig::default()
    };

    // Insert particles per insert block.
    //
    // Insertion runs on EVERY rank. Random insertion is fully deterministic
    // (seeded RNG): every rank generates the identical global packing, but only
    // stores the atoms whose position lies inside its own subdomain. This means
    // each atom is born inside its owner's subdomain, so the per-step `exchange()`
    // only ever needs a single hop. File-based insertion is parsed on every rank
    // and likewise filtered to the local subdomain.
    if let Some(ref inserts) = particles_config.insert {
        {
            // Tags must be globally unique and identical across ranks, so seed the
            // running tag from the global max (reduced) rather than the local max.
            // (No all_reduce_max in the backend, so reduce -max with min and negate.)
            let local_max_tag = atom.get_max_tag() as f64;
            let mut max_tag = (-comm.all_reduce_min_f64(-local_max_tag)) as u32;

            for insert in inserts {
                if insert.source == "file" {
                    // ── File-based insertion ──
                    insert_from_file(
                        insert,
                        &mut atom,
                        &registry,
                        &material_table,
                        &domain,
                        &mut max_tag,
                    )
                    .expect("file insertion was fully parsed during fallible plugin preflight");
                } else if is_rate_insert_config(insert) {
                    // ── Rate-based: register for runtime insertion ──
                    let mat_name = insert.material.as_deref().expect(
                        "rate insertion material was validated during fallible plugin preflight",
                    );
                    let mat_idx = resolve_material(&material_table, mat_name)
                        .expect("rate insertion was validated before setup");
                    let (rate, _, _) =
                        validate_rate_insert_config(insert, "Rate-based [[particles.insert]]")
                            .expect(
                                "rate insertion was validated during fallible plugin preflight",
                            );
                    println!(
                        "DemAtomInsert: registering rate-based insertion for material '{}' (rate={}/every {})",
                        mat_name,
                        rate,
                        insert.rate_interval.unwrap_or(1),
                    );
                    rate_state.entries.push(RateInsertEntry {
                        config: insert.clone(),
                        mat_idx,
                        total_inserted: 0,
                    });
                } else {
                    // ── Immediate random insertion (deterministic, born-in-owner) ──
                    let mat_name = insert.material.as_deref().expect(
                        "random insertion material was validated during fallible plugin preflight",
                    );
                    let mat_idx = resolve_material(&material_table, mat_name)
                        .expect("random insertion was validated before setup");
                    let cutoff_padding = material_table.liquid_bridge_cutoff_padding(mat_idx);
                    let count = insert.count.expect(
                        "random insertion count was validated during fallible plugin preflight",
                    );
                    let radius_spec = insert.radius.as_ref().expect(
                        "random insertion radius was validated during fallible plugin preflight",
                    );
                    let density = insert.density.expect(
                        "random insertion density was validated during fallible plugin preflight",
                    );

                    let max_r = radius_spec.try_max_radius().expect(
                        "random insertion radius was validated during fallible plugin preflight",
                    );
                    if comm.rank() == 0 {
                        println!(
                            "DemAtomInsert: inserting {} particles of material '{}' (r={}, rho={}, E={}, nu={})",
                            count,
                            mat_name,
                            max_r,
                            density,
                            material_table.youngs_mod[mat_idx as usize],
                            material_table.poisson_ratio[mat_idx as usize]
                        );
                    }

                    // Use explicit region or default to domain bounds inset by max radius.
                    let region = insert
                        .region
                        .clone()
                        .map(Ok)
                        .unwrap_or_else(|| {
                            default_insert_region(&domain, max_r, "[[particles.insert]]")
                        })
                        .expect("random insertion region was validated during fallible plugin preflight");

                    // Velocity setup (drawn deterministically per accepted atom).
                    let rand_vel = insert.velocity.unwrap_or(0.0);
                    let normal = (rand_vel > 0.0).then(|| {
                        Normal::new(0.0, rand_vel)
                            .expect("insert velocity was validated before Normal construction")
                    });
                    let vx = insert.velocity_x.unwrap_or(0.0);
                    let vy = insert.velocity_y.unwrap_or(0.0);
                    let vz = insert.velocity_z.unwrap_or(0.0);

                    // Seeded RNG: identical candidate stream on every rank.
                    let seed = insert.seed.unwrap_or(0);
                    let mut rng = StdRng::seed_from_u64(seed);

                    // Replicated scratch of ALL accepted (position, radius) — every
                    // rank maintains the full packing for globally-correct overlap
                    // checks. A spatial hash keeps the cost O(N)/rank instead of
                    // P×O(N²). Cell size ~ 2× max diameter so the 3×3×3 neighborhood
                    // covers every possible overlap.
                    let pbc = PeriodicBox::from_domain(&domain);
                    let cell_size = (2.0 * max_r * 1.1).max(1e-10);
                    let mut spatial_hash = SpatialHash::new(cell_size);
                    let mut all_pos: Vec<[f64; 3]> = Vec::with_capacity(count as usize);
                    let mut all_rad: Vec<f64> = Vec::with_capacity(count as usize);

                    let mut inserted = 0u32;
                    let mut attempts = 0u64;
                    let max_attempts = count as u64 * 1_000_000;
                    while inserted < count && attempts < max_attempts {
                        attempts += 1;
                        // Advance the shared RNG identically on every rank.
                        let Some([x, y, z]) = try_sample_insertion_point(&region, &mut rng) else {
                            continue;
                        };
                        let radius = radius_spec.try_sample(&mut rng)
                            .expect("random insertion radius was validated during fallible plugin preflight");
                        let candidate = [x, y, z];

                        if spatial_hash.has_overlap(&candidate, radius, &all_pos, &all_rad, &pbc) {
                            continue;
                        }

                        // Accepted: draw velocity + assign tag in GLOBAL order so the
                        // RNG advances identically on every rank.
                        let mut vel = [vx, vy, vz];
                        if let Some(ref n) = normal {
                            vel[0] += n.sample(&mut rng);
                            vel[1] += n.sample(&mut rng);
                            vel[2] += n.sample(&mut rng);
                        }
                        let tag = max_tag;
                        max_tag += 1;

                        // Replicate into the global packing scratch (all ranks).
                        let new_idx = all_pos.len();
                        spatial_hash.insert(new_idx, &candidate);
                        all_pos.push(candidate);
                        all_rad.push(radius);

                        // Materialize into the Atom arrays only if this rank owns the
                        // position (half-open interval, matching exchange() ownership).
                        if owns_position(&domain, &candidate) {
                            insert_single_particle(
                                &mut atom,
                                &registry,
                                DemParticle {
                                    pos: candidate,
                                    vel,
                                    radius,
                                    cutoff_padding,
                                    density,
                                    mat_idx,
                                    tag,
                                },
                            );
                        }
                        inserted += 1;
                    }
                    if inserted < count && comm.rank() == 0 {
                        eprintln!(
                            "WARNING: Could only insert {}/{} particles after {} attempts. \
                             Increase domain size or reduce particle count.",
                            inserted, count, max_attempts
                        );
                    }
                }
            }
        }
    }
}

// ── File-based insertion ────────────────────────────────────────────────────

fn insert_from_file(
    insert: &InsertConfig,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let file_path = insert
        .file
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "particle",
            field: "file",
        })?;
    let format = insert
        .format
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "particle",
            field: "format",
        })?;

    match format {
        "csv" => read_csv_particles(
            insert,
            file_path,
            atom,
            registry,
            material_table,
            domain,
            max_tag,
        ),
        "lammps_dump" => read_lammps_dump_particles(
            insert,
            file_path,
            atom,
            registry,
            material_table,
            domain,
            max_tag,
        ),
        "lammps_data" => read_lammps_data_particles(
            insert,
            file_path,
            atom,
            registry,
            material_table,
            domain,
            max_tag,
        ),
        other => Err(InsertFileError::UnknownFormat {
            format: other.to_string(),
        }),
    }
}

fn read_csv_particles(
    insert: &InsertConfig,
    file_path: &str,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let mat_name = insert
        .material
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "CSV",
            field: "material",
        })?;
    let mat_idx = resolve_file_material(material_table, mat_name)?;

    let type_index_map = insert
        .type_map
        .as_ref()
        .map(|tm| resolve_type_map(tm, material_table))
        .transpose()?;

    // Open before checking fields that are needed only to decode rows. This
    // reports a missing input file at the fallible boundary instead of hiding
    // it behind a later configuration omission.
    let file = File::open(file_path).map_err(|e| InsertFileError::FileOpen {
        path: file_path.to_string(),
        source: e.to_string(),
    })?;

    let density = insert.density.ok_or(InsertFileError::MissingField {
        source: "CSV",
        field: "density",
    })?;

    let cols = insert.columns.clone().unwrap_or_default();
    let col_x = cols.x.unwrap_or(0);
    let col_y = cols.y.unwrap_or(1);
    let col_z = cols.z.unwrap_or(2);
    let col_radius = cols.radius;
    let col_vx = cols.vx;
    let col_vy = cols.vy;
    let col_vz = cols.vz;
    let col_atom_type = cols.atom_type;

    let default_radius = match &insert.radius {
        Some(RadiusSpec::Fixed(r)) => Some(*r),
        _ => None,
    };

    let reader = BufReader::new(file);
    let mut count = 0u32;

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| InsertFileError::FileRead {
            path: file_path.to_string(),
            line: line_num + 1,
            source: e.to_string(),
        })?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // Skip header line if it starts with a letter
        if line_num == 0 && trimmed.chars().next().map_or(false, |c| c.is_alphabetic()) {
            continue;
        }

        let fields: Vec<&str> = trimmed.split(',').map(|s| s.trim()).collect();
        let parse = |idx: usize, name: &'static str| -> Result<f64, InsertFileError> {
            fields
                .get(idx)
                .ok_or_else(|| InsertFileError::MissingColumn {
                    path: file_path.to_string(),
                    line: line_num + 1,
                    field: name,
                })
                .and_then(|s| {
                    s.parse()
                        .map_err(|e: std::num::ParseFloatError| InsertFileError::ParseField {
                            path: file_path.to_string(),
                            line: line_num + 1,
                            field: format!("{} (column {})", name, idx),
                            value: (*s).to_string(),
                            source: e.to_string(),
                        })
                })
        };

        let x = parse(col_x, "x")?;
        let y = parse(col_y, "y")?;
        let z = parse(col_z, "z")?;
        let radius = col_radius
            .map(|c| parse(c, "radius"))
            .transpose()?
            .or(default_radius)
            .ok_or(InsertFileError::MissingDefault {
                path: file_path.to_string(),
                field: "radius",
                context: "CSV file insertion with no radius column",
            })?;
        let vx = col_vx.map(|c| parse(c, "vx")).transpose()?.unwrap_or(0.0);
        let vy = col_vy.map(|c| parse(c, "vy")).transpose()?.unwrap_or(0.0);
        let vz = col_vz.map(|c| parse(c, "vz")).transpose()?.unwrap_or(0.0);

        // Determine material: type_map lookup (if atom_type column present) → default material
        let row_mat_idx = match col_atom_type {
            Some(col) => {
                let file_type = parse(col, "atom_type")? as u32;
                lookup_material_for_type(file_type, type_index_map.as_ref(), mat_idx)
            }
            None => mat_idx,
        };
        let cutoff_padding = material_table.liquid_bridge_cutoff_padding(row_mat_idx);

        // Tag advances for every file particle (keeps tags globally consistent
        // across ranks); the atom is only stored if it lies in this subdomain.
        if owns_position(domain, &[x, y, z]) {
            insert_single_particle(
                atom,
                registry,
                DemParticle {
                    pos: [x, y, z],
                    vel: [vx, vy, vz],
                    radius,
                    cutoff_padding,
                    density,
                    mat_idx: row_mat_idx,
                    tag: *max_tag,
                },
            );
            count += 1;
        }
        *max_tag += 1;
    }

    println!(
        "DemAtomInsert: loaded {} local particles from CSV '{}'",
        count, file_path
    );
    Ok(())
}

fn read_lammps_dump_particles(
    insert: &InsertConfig,
    file_path: &str,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let mat_name = insert
        .material
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "lammps_dump",
            field: "material",
        })?;
    let mat_idx = resolve_file_material(material_table, mat_name)?;

    let type_index_map = insert
        .type_map
        .as_ref()
        .map(|tm| resolve_type_map(tm, material_table))
        .transpose()?;

    let density = insert.density.ok_or(InsertFileError::MissingField {
        source: "lammps_dump",
        field: "density",
    })?;

    let default_radius = match &insert.radius {
        Some(RadiusSpec::Fixed(r)) => Some(*r),
        _ => None,
    };

    let file = File::open(file_path).map_err(|e| InsertFileError::FileOpen {
        path: file_path.to_string(),
        source: e.to_string(),
    })?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Parse LAMMPS dump format
    let mut n_atoms: usize = 0;
    let mut column_names: Vec<String> = Vec::new();
    let mut reading_atoms = false;
    let mut count = 0u32;

    // Helper to find column index by name
    let find_col =
        |names: &[String], name: &str| -> Option<usize> { names.iter().position(|n| n == name) };

    let mut line_num = 0usize;
    while let Some(line) = lines.next() {
        line_num += 1;
        let line = line.map_err(|e| InsertFileError::FileRead {
            path: file_path.to_string(),
            line: line_num,
            source: e.to_string(),
        })?;
        let trimmed = line.trim();

        if trimmed == "ITEM: NUMBER OF ATOMS" {
            if let Some(next) = lines.next() {
                line_num += 1;
                let next = next.map_err(|e| InsertFileError::FileRead {
                    path: file_path.to_string(),
                    line: line_num,
                    source: e.to_string(),
                })?;
                n_atoms = next.trim().parse().map_err(|e: std::num::ParseIntError| {
                    InsertFileError::ParseField {
                        path: file_path.to_string(),
                        line: line_num,
                        field: "number of atoms".to_string(),
                        value: next.trim().to_string(),
                        source: e.to_string(),
                    }
                })?;
            }
            continue;
        }

        if trimmed.starts_with("ITEM: ATOMS") {
            // Parse column names from header: "ITEM: ATOMS id type x y z ..."
            column_names = trimmed
                .strip_prefix("ITEM: ATOMS")
                .unwrap_or("")
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            reading_atoms = true;
            continue;
        }

        if trimmed.starts_with("ITEM:") {
            reading_atoms = false;
            continue;
        }

        if reading_atoms && !trimmed.is_empty() {
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() < column_names.len() {
                continue;
            }

            let parse_col = |name: &'static str| -> Result<Option<f64>, InsertFileError> {
                let Some(i) = find_col(&column_names, name) else {
                    return Ok(None);
                };
                let Some(value) = fields.get(i) else {
                    return Err(InsertFileError::MissingColumn {
                        path: file_path.to_string(),
                        line: line_num,
                        field: name,
                    });
                };
                value
                    .parse()
                    .map(Some)
                    .map_err(|e: std::num::ParseFloatError| InsertFileError::ParseField {
                        path: file_path.to_string(),
                        line: line_num,
                        field: name.to_string(),
                        value: (*value).to_string(),
                        source: e.to_string(),
                    })
            };

            let x = parse_col("x")?.ok_or(InsertFileError::MissingColumn {
                path: file_path.to_string(),
                line: line_num,
                field: "x",
            })?;
            let y = parse_col("y")?.ok_or(InsertFileError::MissingColumn {
                path: file_path.to_string(),
                line: line_num,
                field: "y",
            })?;
            let z = parse_col("z")?.ok_or(InsertFileError::MissingColumn {
                path: file_path.to_string(),
                line: line_num,
                field: "z",
            })?;
            let vx = parse_col("vx")?.unwrap_or(0.0);
            let vy = parse_col("vy")?.unwrap_or(0.0);
            let vz = parse_col("vz")?.unwrap_or(0.0);
            let radius =
                parse_col("radius")?
                    .or(default_radius)
                    .ok_or(InsertFileError::MissingDefault {
                        path: file_path.to_string(),
                        field: "radius",
                        context: "LAMMPS dump file insertion with no radius column",
                    })?;

            // Determine material: type_map override → default material
            let row_mat_idx = match parse_col("type")? {
                Some(t) => lookup_material_for_type(t as u32, type_index_map.as_ref(), mat_idx),
                None => mat_idx,
            };
            let cutoff_padding = material_table.liquid_bridge_cutoff_padding(row_mat_idx);

            if owns_position(domain, &[x, y, z]) {
                insert_single_particle(
                    atom,
                    registry,
                    DemParticle {
                        pos: [x, y, z],
                        vel: [vx, vy, vz],
                        radius,
                        cutoff_padding,
                        density,
                        mat_idx: row_mat_idx,
                        tag: *max_tag,
                    },
                );
                count += 1;
            }
            *max_tag += 1;
        }
    }

    let _ = n_atoms; // used for format validation if needed
    println!(
        "DemAtomInsert: loaded {} local particles from LAMMPS dump '{}'",
        count, file_path
    );
    Ok(())
}

/// Parse a field from a LAMMPS data file, with a user-friendly error on failure.
fn parse_field<T: std::str::FromStr>(
    value: &str,
    field_name: &str,
    line_num: usize,
    file_path: &str,
) -> Result<T, InsertFileError>
where
    T::Err: std::fmt::Display,
{
    value.parse::<T>().map_err(|e| InsertFileError::ParseField {
        path: file_path.to_string(),
        line: line_num,
        field: field_name.to_string(),
        value: value.to_string(),
        source: e.to_string(),
    })
}

fn read_lammps_data_particles(
    insert: &InsertConfig,
    file_path: &str,
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    domain: &Domain,
    max_tag: &mut u32,
) -> Result<(), InsertFileError> {
    let mat_name = insert
        .material
        .as_deref()
        .ok_or(InsertFileError::MissingField {
            source: "lammps_data",
            field: "material",
        })?;
    let mat_idx = resolve_file_material(material_table, mat_name)?;

    let type_index_map = insert
        .type_map
        .as_ref()
        .map(|tm| resolve_type_map(tm, material_table))
        .transpose()?;

    let default_density = insert.density;
    let default_radius = match &insert.radius {
        Some(RadiusSpec::Fixed(r)) => Some(*r),
        _ => None,
    };

    let file = File::open(file_path).map_err(|e| InsertFileError::FileOpen {
        path: file_path.to_string(),
        source: e.to_string(),
    })?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader
        .lines()
        .enumerate()
        .map(|(i, l)| {
            l.map_err(|e| InsertFileError::FileRead {
                path: file_path.to_string(),
                line: i + 1,
                source: e.to_string(),
            })
        })
        .collect::<Result<_, _>>()?;

    // Detect atom style from config or from "Atoms # style" header
    let config_style = insert.atom_style.as_deref();

    // Find section start indices
    let mut atoms_start = None;
    let mut atoms_style = None;
    let mut velocities_start = None;

    for (i, line) in lines.iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.starts_with("Atoms") {
            atoms_start = Some(i + 1);
            // Try to detect style from "Atoms # style" comment
            if let Some(comment) = trimmed.strip_prefix("Atoms") {
                let comment = comment.trim();
                if let Some(style) = comment.strip_prefix('#') {
                    let style = style.trim();
                    if !style.is_empty() {
                        atoms_style = Some(style.to_string());
                    }
                }
            }
        } else if trimmed == "Velocities" {
            velocities_start = Some(i + 1);
        }
    }

    let atom_style = config_style
        .map(|s| s.to_string())
        .or(atoms_style)
        .unwrap_or_else(|| "atomic".to_string());

    let atoms_start = atoms_start.ok_or(InsertFileError::MissingSection {
        path: file_path.to_string(),
        section: "Atoms",
    })?;

    // Parse Atoms section
    struct ParsedAtom {
        id: u32,
        atom_type: u32,
        pos: [f64; 3],
        radius: f64,
        density: f64,
    }

    let section_headers = [
        "Atoms",
        "Velocities",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "Masses",
        "Pair Coeffs",
    ];
    let is_section_header = |line: &str| -> bool {
        let trimmed = line.trim();
        section_headers.iter().any(|h| trimmed.starts_with(h))
    };

    let mut parsed_atoms: Vec<ParsedAtom> = Vec::new();

    for i in atoms_start..lines.len() {
        let trimmed = lines[i].trim();
        if trimmed.is_empty() {
            continue;
        }
        if is_section_header(trimmed) {
            break;
        }
        // Skip comment lines
        if trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();

        match atom_style.as_str() {
            "atomic" => {
                // id type x y z
                if fields.len() < 5 {
                    return Err(InsertFileError::RowTooShort {
                        path: file_path.to_string(),
                        line: i + 1,
                        style: "atomic".to_string(),
                        expected: 5,
                        found: fields.len(),
                    });
                }
                let id: u32 = parse_field(fields[0], "atom id", i + 1, file_path)?;
                let atype: u32 = parse_field(fields[1], "atom type", i + 1, file_path)?;
                let x: f64 = parse_field(fields[2], "x coordinate", i + 1, file_path)?;
                let y: f64 = parse_field(fields[3], "y coordinate", i + 1, file_path)?;
                let z: f64 = parse_field(fields[4], "z coordinate", i + 1, file_path)?;
                let radius = default_radius.ok_or(InsertFileError::MissingDefault {
                    path: file_path.to_string(),
                    field: "radius",
                    context: "atomic style LAMMPS data",
                })?;
                let density = default_density.ok_or(InsertFileError::MissingDefault {
                    path: file_path.to_string(),
                    field: "density",
                    context: "atomic style LAMMPS data",
                })?;
                parsed_atoms.push(ParsedAtom {
                    id,
                    atom_type: atype,
                    pos: [x, y, z],
                    radius,
                    density,
                });
            }
            "sphere" | "bpm/sphere" => {
                // id type diameter density x y z
                if fields.len() < 7 {
                    return Err(InsertFileError::RowTooShort {
                        path: file_path.to_string(),
                        line: i + 1,
                        style: atom_style.clone(),
                        expected: 7,
                        found: fields.len(),
                    });
                }
                let id: u32 = parse_field(fields[0], "atom id", i + 1, file_path)?;
                let atype: u32 = parse_field(fields[1], "atom type", i + 1, file_path)?;
                let diameter: f64 = parse_field(fields[2], "diameter", i + 1, file_path)?;
                let density: f64 = parse_field(fields[3], "density", i + 1, file_path)?;
                let x: f64 = parse_field(fields[4], "x coordinate", i + 1, file_path)?;
                let y: f64 = parse_field(fields[5], "y coordinate", i + 1, file_path)?;
                let z: f64 = parse_field(fields[6], "z coordinate", i + 1, file_path)?;
                parsed_atoms.push(ParsedAtom {
                    id,
                    atom_type: atype,
                    pos: [x, y, z],
                    radius: diameter / 2.0,
                    density,
                });
            }
            other => {
                return Err(InsertFileError::UnsupportedAtomStyle {
                    path: file_path.to_string(),
                    style: other.to_string(),
                });
            }
        }
    }

    // Parse Velocities section (optional) — build id → [vx, vy, vz] map
    let mut velocity_map: HashMap<u32, [f64; 3]> = HashMap::new();
    if let Some(vel_start) = velocities_start {
        for i in vel_start..lines.len() {
            let trimmed = lines[i].trim();
            if trimmed.is_empty() {
                continue;
            }
            if is_section_header(trimmed) {
                break;
            }
            if trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() >= 4 {
                let id: u32 = parse_field(fields[0], "atom id (Velocities)", i + 1, file_path)?;
                let vx: f64 = parse_field(fields[1], "vx", i + 1, file_path)?;
                let vy: f64 = parse_field(fields[2], "vy", i + 1, file_path)?;
                let vz: f64 = parse_field(fields[3], "vz", i + 1, file_path)?;
                velocity_map.insert(id, [vx, vy, vz]);
            }
        }
    }

    // Insert all parsed atoms (only those owned by this subdomain).
    let mut count = 0usize;
    for pa in parsed_atoms {
        let vel = velocity_map.get(&pa.id).copied().unwrap_or([0.0; 3]);
        let row_mat_idx = lookup_material_for_type(pa.atom_type, type_index_map.as_ref(), mat_idx);
        let cutoff_padding = material_table.liquid_bridge_cutoff_padding(row_mat_idx);
        if owns_position(domain, &pa.pos) {
            insert_single_particle(
                atom,
                registry,
                DemParticle {
                    pos: pa.pos,
                    vel,
                    radius: pa.radius,
                    cutoff_padding,
                    density: pa.density,
                    mat_idx: row_mat_idx,
                    tag: *max_tag,
                },
            );
            count += 1;
        }
        *max_tag += 1;
    }

    println!(
        "DemAtomInsert: loaded {} local particles from LAMMPS data file '{}' (style: {})",
        count, file_path, atom_style
    );
    Ok(())
}

// ── Update system: rate-based insertion ─────────────────────────────────────

/// Update system for rate-based particle insertion during the simulation run.
///
/// Checks each registered [`RateInsertEntry`] against the current timestep, interval,
/// start/end bounds, and total limit. Uses a [`SpatialHash`] for O(1) overlap detection
/// when placing new particles. Runs in `ParticleSimScheduleSet::PreInitialIntegration`.
#[allow(clippy::too_many_arguments)]
pub fn dem_rate_insert(
    comm: Res<CommResource>,
    domain: Res<Domain>,
    mut atom: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    material_table: Res<MaterialTable>,
    mut rate_state: ResMut<RateInsertState>,
    mut comm_state: ResMut<CurrentState<CommState>>,
) {
    // Rate insertion runs on EVERY rank (born-in-owner): each rank generates the
    // identical candidate stream from a step-derived seed and stores only the
    // candidates that fall inside its own subdomain, so a new atom is born inside
    // its owner and never needs a multi-hop exchange.
    if rate_state.entries.is_empty() {
        return;
    }

    let step = run_state.total_cycle;
    let mut any_to_insert = false;

    // Quick check if any entry needs insertion this step (before stripping ghosts)

    // Quick check if any entry needs insertion this step
    for entry in rate_state.entries.iter() {
        let interval = entry.config.rate_interval.unwrap_or(1);
        let start = entry.config.rate_start.unwrap_or(0);
        if step < start {
            continue;
        }
        if let Some(end) = entry.config.rate_end {
            if step > end {
                continue;
            }
        }
        if let Some(limit) = entry.config.rate_limit {
            if entry.total_inserted >= limit {
                continue;
            }
        }
        let steps_since_start = step - start;
        if interval == 0 || steps_since_start % interval == 0 {
            any_to_insert = true;
            break;
        }
    }

    if !any_to_insert {
        return;
    }

    // Strip ghost atoms before inserting new local atoms.
    // New atoms are appended at atom.len(), which must equal nlocal so that
    // the subsequent borders() truncate_to_nlocal() doesn't discard them.
    if atom.nghost > 0 {
        // Keep the core ghost suffix and every DIRT extension in one
        // transaction.  In particular, a plugin registered after setup has a
        // row here too; truncating Atom and AtomDataRegistry separately would
        // leave a panic/error window between the two mutations.
        ParticleStore::new(&mut atom, &registry)
            .discard_ghosts()
            .expect("rate insertion requires an aligned local/ghost particle layout");
    }

    // Base tag must be globally consistent across ranks. (No all_reduce_max in
    // the backend, so reduce -max with min and negate.) Each attempt this step
    // consumes one tag slot regardless of acceptance so tags stay unique across
    // ranks without any extra collective.
    let local_max_tag = if atom.tag.is_empty() {
        -1.0
    } else {
        atom.get_max_tag() as f64
    };
    let base_tag = (-comm.all_reduce_min_f64(-local_max_tag)) as i64 + 1;
    let mut tag_cursor: u32 = base_tag.max(0) as u32;

    for entry_idx in 0..rate_state.entries.len() {
        let interval = rate_state.entries[entry_idx]
            .config
            .rate_interval
            .unwrap_or(1);
        let start = rate_state.entries[entry_idx].config.rate_start.unwrap_or(0);
        let (rate, radius_spec, density) = validate_rate_insert_config(
            &rate_state.entries[entry_idx].config,
            "rate-based [[particles.insert]]",
        )
        .expect("rate insertion was validated during fallible plugin preflight");

        if step < start {
            continue;
        }
        if let Some(end) = rate_state.entries[entry_idx].config.rate_end {
            if step > end {
                continue;
            }
        }
        if let Some(limit) = rate_state.entries[entry_idx].config.rate_limit {
            if rate_state.entries[entry_idx].total_inserted >= limit {
                continue;
            }
        }
        let steps_since_start = step - start;
        if interval > 0 && steps_since_start % interval != 0 {
            continue;
        }

        // How many to insert this step
        let mut to_insert = rate;
        if let Some(limit) = rate_state.entries[entry_idx].config.rate_limit {
            let remaining = limit - rate_state.entries[entry_idx].total_inserted;
            to_insert = to_insert.min(remaining);
        }

        let mat_idx = rate_state.entries[entry_idx].mat_idx;
        let cutoff_padding = material_table.liquid_bridge_cutoff_padding(mat_idx);

        let max_r = radius_spec
            .try_max_radius()
            .expect("rate insertion radius was validated during fallible plugin preflight");
        let region = rate_state.entries[entry_idx]
            .config
            .region
            .clone()
            .map(Ok)
            .unwrap_or_else(|| {
                default_insert_region(&domain, max_r, "rate-based [[particles.insert]]")
            })
            .expect("rate insertion region was validated during fallible plugin preflight");

        // Velocity parameters (drawn deterministically per accepted candidate).
        let config_seed = rate_state.entries[entry_idx].config.seed.unwrap_or(0);
        let rand_vel = rate_state.entries[entry_idx].config.velocity.unwrap_or(0.0);
        let vel_normal = (rand_vel > 0.0).then(|| {
            Normal::new(0.0, rand_vel)
                .expect("rate insertion velocity was validated before Normal construction")
        });
        let vx = rate_state.entries[entry_idx]
            .config
            .velocity_x
            .unwrap_or(0.0);
        let vy = rate_state.entries[entry_idx]
            .config
            .velocity_y
            .unwrap_or(0.0);
        let vz = rate_state.entries[entry_idx]
            .config
            .velocity_z
            .unwrap_or(0.0);

        // Seed the candidate stream from (config seed, step, entry) so it is
        // identical on every rank yet varies between insertion events.
        let mut rng = StdRng::seed_from_u64(
            config_seed
                ^ (step as u64).wrapping_mul(0x9E3779B97F4A7C15)
                ^ (entry_idx as u64).wrapping_mul(0xD1B54A32D192ED03),
        );

        // Replicated overlap scratch: ONLY the positions/radii accepted THIS step
        // (the global set of new atoms). It is identical on every rank because the
        // candidate stream is seeded identically and the scratch is the same on
        // all ranks. This is what keeps accept/reject — and therefore the RNG
        // advancement and the global accept count `inserted` — in lock-step across
        // ranks, so the collective borders()/exchange() triggered below stay
        // synchronized.
        //
        // NOTE: existing local atoms are deliberately NOT added to the scratch.
        // They differ per rank, so including them would make accept/reject (and
        // the RNG stream) diverge across ranks and desync the collectives. New
        // atoms may therefore be born overlapping already-present particles; the
        // contact model resolves that initial overlap via repulsion. Rate-insert
        // regions are normally placed in free space (e.g. above a settled bed),
        // so this is rare in practice.
        let cell_size = (2.0 * max_r * 1.1).max(1e-10);
        let pbc = PeriodicBox::from_domain(&domain);
        let mut spatial_hash = SpatialHash::new(cell_size);
        let mut all_pos: Vec<[f64; 3]> = Vec::new();
        let mut all_rad: Vec<f64> = Vec::new();

        let mut inserted = 0u32; // accepted globally
        let mut local_inserted = 0u32; // stored on this rank
        let mut attempts = 0u32;
        let max_attempts = to_insert * 100;

        while inserted < to_insert && attempts < max_attempts {
            // Every attempt consumes one tag slot so tags stay globally unique.
            let tag = tag_cursor;
            tag_cursor = tag_cursor.wrapping_add(1);
            attempts += 1;

            let Some([x, y, z]) = try_sample_insertion_point(&region, &mut rng) else {
                continue;
            };
            let radius = radius_spec
                .try_sample(&mut rng)
                .expect("rate insertion radius was validated during fallible plugin preflight");
            let candidate = [x, y, z];

            if spatial_hash.has_overlap(&candidate, radius, &all_pos, &all_rad, &pbc) {
                continue;
            }

            // Accepted: draw velocity (advances RNG identically on every rank).
            let mut vel = [vx, vy, vz];
            if let Some(ref n) = vel_normal {
                vel[0] += n.sample(&mut rng);
                vel[1] += n.sample(&mut rng);
                vel[2] += n.sample(&mut rng);
            }

            // Replicate into the global new-atom scratch on every rank.
            let scratch_idx = all_pos.len();
            spatial_hash.insert(scratch_idx, &candidate);
            all_pos.push(candidate);
            all_rad.push(radius);

            // Store only if this rank owns the position.
            if owns_position(&domain, &candidate) {
                insert_single_particle(
                    &mut atom,
                    &registry,
                    DemParticle {
                        pos: candidate,
                        vel,
                        radius,
                        cutoff_padding,
                        density,
                        mat_idx,
                        tag,
                    },
                );
                local_inserted += 1;
            }
            inserted += 1;
        }

        rate_state.entries[entry_idx].total_inserted += inserted;
        let _ = local_inserted;

        if inserted > 0 {
            // Force full ghost rebuild on EVERY rank if any rank inserted, so the
            // collective borders()/exchange() stay in lock-step. `inserted` is the
            // global accept count and is identical on all ranks (replicated stream).
            comm_state.0 = CommState::FullRebuild;
        }
        if inserted > 0 && attempts >= max_attempts && comm.rank() == 0 {
            eprintln!(
                "WARNING: Rate insertion at step {} only placed {}/{} particles (max attempts reached)",
                step, inserted, to_insert
            );
        }
    }
}

// ── Delta time calculation ──────────────────────────────────────────────────

/// Computes a stable timestep from the Rayleigh wave speed criterion.
///
/// For each particle, estimates the Rayleigh wave transit time across the particle
/// diameter using `dt_R = π·r / α · √(ρ/G)`, where α ≈ 0.1631·ν + 0.8766 and
/// G = E / (2·(1+ν)). The final timestep is 15% of the minimum across all particles.
fn calculate_delta_time(
    comm: Res<CommResource>,
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    run_config: Res<RunConfig>,
    scheduler_manager: Res<SchedulerManager>,
) {
    // If the current stage specifies an explicit dt, use it directly.
    let index = scheduler_manager.index;
    let config_dt = run_config.current_stage(index).dt;
    if config_dt > 0.0 {
        atoms.dt = config_dt;
        if comm.rank() == 0 {
            println!("Using {} for delta time (from config)", config_dt);
        }
        return;
    }

    // Auto-compute from Rayleigh wave speed criterion.
    let dem = registry.expect::<DemAtom>("calculate_delta_time");
    let mut dt: f64 = 0.001;

    for i in 0..dem.radius.len() {
        let mat_idx = atoms.atom_type[i] as usize;
        let youngs_mod = material_table.youngs_mod[mat_idx];
        let poisson_ratio = material_table.poisson_ratio[mat_idx];
        let g = youngs_mod / (2.0 * (1.0 + poisson_ratio));
        let alpha = 0.1631 * poisson_ratio + 0.876605;
        let delta = PI * dem.radius[i] / alpha * (dem.density[i] / g).sqrt();
        dt = delta.min(dt);
    }

    dt = comm.all_reduce_min_f64(dt);

    if comm.rank() == 0 {
        println!("Using {} for delta time", dt * 0.15);
    }
    atoms.dt = dt * 0.15;
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RadiusDistribution;
    use soil_core::{toml, AtomData, AtomDataRegistry, ParticleStoreError};
    use soil_derive::AtomData;

    /// A second extension registered after construction.  This represents an
    /// optional DIRT plugin arriving after particles were inserted.
    #[derive(Default, AtomData)]
    struct LateProbe {
        rows: Vec<f64>,
    }

    /// Deliberately refuses to create a default row, exercising the facade's
    /// rollback boundary from the DIRT insertion caller's side.
    #[derive(Default)]
    struct BrokenDefaults;

    impl AtomData for BrokenDefaults {
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }
        fn as_any_mut(&mut self) -> &mut dyn std::any::Any {
            self
        }
        fn snapshot(&self) -> Box<dyn AtomData> {
            Box::new(Self)
        }
        fn len(&self) -> usize {
            0
        }
        fn push_default(&mut self) {}
        fn truncate(&mut self, _: usize) {}
        fn swap_remove(&mut self, _: usize) {}
        fn pack(&self, _: usize, _: &mut Vec<f64>) {}
        fn unpack(&mut self, _: &[f64]) -> usize {
            0
        }
        fn apply_permutation(&mut self, _: &[usize], _: usize) {}
    }

    fn test_dem_registry() -> AtomDataRegistry {
        let mut registry = AtomDataRegistry::new();
        registry.try_register(DemAtom::new(), 0).unwrap();
        registry
    }

    fn rate_config(seed: u64) -> InsertConfig {
        InsertConfig {
            source: "random".to_string(),
            material: Some("glass".to_string()),
            count: None,
            radius: Some(RadiusSpec::Fixed(0.001)),
            density: Some(2500.0),
            velocity: None,
            velocity_x: Some(0.0),
            velocity_y: Some(0.0),
            velocity_z: Some(0.0),
            region: Some(Region::Block {
                min: [0.0; 3],
                max: [1.0; 3],
            }),
            rate: Some(8),
            rate_interval: Some(1),
            rate_start: Some(0),
            rate_end: Some(0),
            rate_limit: Some(8),
            file: None,
            format: None,
            columns: None,
            type_map: None,
            atom_style: None,
            seed: Some(seed),
        }
    }

    fn run_rate_once(
        mut atom: Atom,
        registry: AtomDataRegistry,
        domain: Domain,
        seed: u64,
    ) -> (Atom, usize, u32) {
        let mut materials = MaterialTable::new();
        materials.add_material("glass", 8.7e9, 0.3, 0.9, 0.5, 0.0, 0.0);
        materials.build_pair_tables();
        let mut app = App::new();
        app.add_resource(CommResource(Box::new(soil_core::SingleProcessComm::new())));
        app.add_resource(domain);
        app.add_resource(std::mem::take(&mut atom));
        app.add_resource(registry);
        app.add_resource(RunState::new());
        app.add_resource(materials);
        app.add_resource(RateInsertState {
            entries: vec![RateInsertEntry {
                config: rate_config(seed),
                mat_idx: 0,
                total_inserted: 0,
            }],
        });
        app.add_resource(CurrentState(CommState::CommunicateOnly));
        app.add_update_system(
            dem_rate_insert,
            ParticleSimScheduleSet::PreInitialIntegration,
        );
        app.organize_systems();
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap().clone();
        let late_rows = app
            .get_resource_ref::<AtomDataRegistry>()
            .unwrap()
            .get::<LateProbe>()
            .map_or(0, |probe| probe.rows.len());
        let inserted = app.get_resource_ref::<RateInsertState>().unwrap().entries[0].total_inserted;
        (atom, late_rows, inserted)
    }

    fn unit_domain(low_x: f64, high_x: f64) -> Domain {
        let mut domain = Domain::new();
        domain.sub_domain_low = [low_x, 0.0, 0.0];
        domain.sub_domain_high = [high_x, 1.0, 1.0];
        domain.boundaries_low = [0.0; 3];
        domain.boundaries_high = [1.0; 3];
        domain.size = [1.0; 3];
        domain.sub_length = [high_x - low_x, 1.0, 1.0];
        domain.volume = high_x - low_x;
        domain
    }

    #[test]
    fn production_rate_insert_partitions_exact_deterministic_tag_rows_and_cleans_late_ghosts() {
        // This invokes the public production system, not its materialization helper.
        // First create the exact local+ghost layout that a rate step receives.
        let mut atom = Atom::new();
        let mut registry = test_dem_registry();
        insert_single_particle(
            &mut atom,
            &registry,
            DemParticle {
                pos: [0.1, 0.5, 0.5],
                vel: [0.0; 3],
                radius: 0.001,
                cutoff_padding: 0.0,
                density: 2500.0,
                mat_idx: 0,
                tag: 41,
            },
        );
        let mut ghost = Atom::new();
        let ghost_registry = test_dem_registry();
        insert_single_particle(
            &mut ghost,
            &ghost_registry,
            DemParticle {
                pos: [0.9, 0.5, 0.5],
                vel: [0.0; 3],
                radius: 0.003,
                cutoff_padding: 0.0,
                density: 2500.0,
                mat_idx: 0,
                tag: 42,
            },
        );
        let mut packed = Vec::new();
        ParticleStore::new(&mut ghost, &ghost_registry)
            .pack_migrant(0, &mut packed)
            .unwrap();
        ParticleStore::new(&mut atom, &registry)
            .append_ghost_records(&packed, 1)
            .unwrap();
        registry
            .try_register(LateProbe::default(), atom.len())
            .unwrap();

        let (atom, late_rows, inserted) =
            run_rate_once(atom, registry, unit_domain(0.0, 1.0), 20260712);
        assert_eq!(inserted, 8);
        assert_eq!((atom.nlocal, atom.nghost, atom.len()), (9, 0, 9));
        assert_eq!(atom.tag[0], 41, "the local prefix survives ghost removal");
        assert_eq!(late_rows, atom.len());

        let (_, _, full_inserted) =
            run_rate_once(Atom::new(), test_dem_registry(), unit_domain(0.0, 1.0), 99);
        assert_eq!(full_inserted, 8);
        let (low, _, _) =
            run_rate_once(Atom::new(), test_dem_registry(), unit_domain(0.0, 0.5), 99);
        let (high, _, _) =
            run_rate_once(Atom::new(), test_dem_registry(), unit_domain(0.5, 1.0), 99);
        let (full, _, _) =
            run_rate_once(Atom::new(), test_dem_registry(), unit_domain(0.0, 1.0), 99);
        let mut partitioned: Vec<_> = low
            .tag
            .iter()
            .zip(&low.pos)
            .chain(high.tag.iter().zip(&high.pos))
            .map(|(tag, pos)| (*tag, *pos))
            .collect();
        let mut serial: Vec<_> = full
            .tag
            .iter()
            .zip(&full.pos)
            .map(|(tag, pos)| (*tag, *pos))
            .collect();
        partitioned.sort_by_key(|row| row.0);
        serial.sort_by_key(|row| row.0);
        assert_eq!(
            partitioned, serial,
            "fixed-seed production rows/tags must be rank-count invariant"
        );
        for (_, pos) in &partitioned {
            assert_eq!(
                owns_position(&unit_domain(0.0, 0.5), pos) as u8
                    + owns_position(&unit_domain(0.5, 1.0), pos) as u8,
                1
            );
        }
    }

    #[test]
    fn particle_store_construction_covers_immediate_and_rate_rows() {
        let mut atoms = Atom::new();
        let registry = test_dem_registry();

        // This is the shared materialization endpoint reached by both the
        // immediate and periodic rate candidate loops.  Use distinct defaults
        // so a future path-specific field regression is observable here.
        for (tag, pos, velocity, radius, material) in [
            (7, [0.1, 0.2, 0.3], [1.0, 0.0, -1.0], 0.002, 3),
            (8, [0.4, 0.5, 0.6], [0.0, 2.0, -2.0], 0.003, 4),
        ] {
            insert_single_particle(
                &mut atoms,
                &registry,
                DemParticle {
                    pos,
                    vel: velocity,
                    radius,
                    cutoff_padding: 0.0004,
                    density: 2500.0,
                    mat_idx: material,
                    tag,
                },
            );
        }

        let dem = registry.expect::<DemAtom>("particle-store construction test");
        assert_eq!((atoms.nlocal, atoms.nghost, atoms.natoms), (2, 0, 2));
        assert_eq!(atoms.tag, vec![7, 8]);
        assert_eq!(atoms.atom_type, vec![3, 4]);
        assert_eq!(atoms.cutoff_radius, vec![0.002 + 0.0004, 0.003 + 0.0004]);
        assert_eq!(dem.radius, vec![0.002, 0.003]);
        assert_eq!(dem.body_id, vec![0.0, 0.0]);
        assert_eq!(dem.quaternion, vec![[1.0, 0.0, 0.0, 0.0]; 2]);
        assert!(registry.validate_rows(atoms.len()));
    }

    #[test]
    fn particle_store_construction_backfills_late_extensions_and_rolls_back() {
        let mut atoms = Atom::new();
        let mut registry = test_dem_registry();
        insert_single_particle(
            &mut atoms,
            &registry,
            DemParticle {
                pos: [0.0; 3],
                vel: [0.0; 3],
                radius: 0.001,
                cutoff_padding: 0.0,
                density: 2500.0,
                mat_idx: 0,
                tag: 1,
            },
        );
        registry
            .try_register(LateProbe::default(), atoms.len())
            .unwrap();
        assert_eq!(
            registry.expect::<LateProbe>("late extension").rows,
            vec![0.0]
        );

        insert_single_particle(
            &mut atoms,
            &registry,
            DemParticle {
                pos: [1.0; 3],
                vel: [0.0; 3],
                radius: 0.001,
                cutoff_padding: 0.0,
                density: 2500.0,
                mat_idx: 0,
                tag: 2,
            },
        );
        assert_eq!(
            registry.expect::<LateProbe>("late extension").rows,
            vec![0.0, 0.0]
        );
        assert!(registry.validate_rows(atoms.len()));

        let mut rollback_atoms = Atom::new();
        let mut rollback_registry = AtomDataRegistry::new();
        rollback_registry.try_register(BrokenDefaults, 0).unwrap();
        assert_eq!(
            ParticleStore::new(&mut rollback_atoms, &rollback_registry).push_default_local(1),
            Err(ParticleStoreError::MalformedExtensionRecord)
        );
        assert!(rollback_atoms.is_empty());
        assert_eq!((rollback_atoms.nlocal, rollback_atoms.natoms), (0, 0));
        assert!(rollback_registry.validate_rows(0));
    }

    #[test]
    fn particle_store_restart_rejection_preserves_dem_construction() {
        let mut atoms = Atom::new();
        let registry = test_dem_registry();
        insert_single_particle(
            &mut atoms,
            &registry,
            DemParticle {
                pos: [0.25; 3],
                vel: [0.0; 3],
                radius: 0.001,
                cutoff_padding: 0.0,
                density: 2500.0,
                mat_idx: 0,
                tag: 17,
            },
        );
        let before_tags = atoms.tag.clone();
        let before_radius = registry
            .expect::<DemAtom>("restart snapshot")
            .radius
            .clone();
        let mut malformed = atoms.clone();
        malformed.mass.clear();
        assert_eq!(
            ParticleStore::new(&mut atoms, &registry).replace_from_restart(malformed, &[]),
            Err(ParticleStoreError::InvalidStructuralOperation)
        );
        assert_eq!(atoms.tag, before_tags);
        assert_eq!(
            registry.expect::<DemAtom>("restart rollback").radius,
            before_radius
        );
        assert!(registry.validate_rows(atoms.len()));
    }

    #[test]
    fn rate_insertion_ghost_cleanup_keeps_dem_rows_synchronized() {
        // Reproduce the layout rate insertion sees after a communication pass:
        // one local row followed by a received ghost carrying real DemAtom
        // fields.  This uses the SOIL framing path rather than manufacturing a
        // matching extension vector by hand.
        let source_registry = test_dem_registry();
        let mut source = Atom::new();
        insert_single_particle(
            &mut source,
            &source_registry,
            DemParticle {
                pos: [0.75, 0.0, 0.0],
                vel: [0.0; 3],
                radius: 0.003,
                cutoff_padding: 0.0,
                density: 2500.0,
                mat_idx: 2,
                tag: 22,
            },
        );
        let mut packed = Vec::new();
        ParticleStore::new(&mut source, &source_registry)
            .pack_migrant(0, &mut packed)
            .unwrap();

        let registry = test_dem_registry();
        let mut atoms = Atom::new();
        insert_single_particle(
            &mut atoms,
            &registry,
            DemParticle {
                pos: [0.25, 0.0, 0.0],
                vel: [0.0; 3],
                radius: 0.001,
                cutoff_padding: 0.0,
                density: 1000.0,
                mat_idx: 1,
                tag: 11,
            },
        );
        ParticleStore::new(&mut atoms, &registry)
            .append_ghost_records(&packed, 1)
            .unwrap();
        assert_eq!((atoms.nlocal, atoms.nghost), (1, 1));
        assert_eq!(
            registry.expect::<DemAtom>("ghost setup").radius,
            vec![0.001, 0.003]
        );

        ParticleStore::new(&mut atoms, &registry)
            .discard_ghosts()
            .unwrap();
        assert_eq!((atoms.nlocal, atoms.nghost, atoms.len()), (1, 0, 1));
        assert_eq!(
            registry.expect::<DemAtom>("ghost cleanup").radius,
            vec![0.001]
        );
        assert!(registry.validate_rows(atoms.len()));
    }

    #[test]
    fn ownership_partition_is_exact_at_multirank_boundaries() {
        let mut low = Domain::new();
        low.sub_domain_low = [0.0, 0.0, 0.0];
        low.sub_domain_high = [0.5, 1.0, 1.0];
        let mut high = Domain::new();
        high.sub_domain_low = [0.5, 0.0, 0.0];
        high.sub_domain_high = [1.0, 1.0, 1.0];
        for point in [
            [0.0, 0.5, 0.5],
            [0.499999, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [0.999999, 0.5, 0.5],
        ] {
            assert_eq!(
                owns_position(&low, &point) as u8 + owns_position(&high, &point) as u8,
                1
            );
        }
    }

    #[test]
    fn bounded_sampling_draws_a_candidate_and_exhaustion_rejects_one() {
        let block = Region::Block {
            min: [-1.0, -2.0, -3.0],
            max: [1.0, 2.0, 3.0],
        };
        let mut immediate_rng = StdRng::seed_from_u64(17);
        let point = try_sample_insertion_point(&block, &mut immediate_rng)
            .expect("a bounded region must produce an immediate-insertion candidate");
        assert!(block.contains(&point));

        let disjoint = Region::Intersect {
            regions: vec![
                Region::Sphere {
                    center: [0.0, 0.0, 0.0],
                    radius: 1.0,
                },
                Region::Sphere {
                    center: [4.0, 0.0, 0.0],
                    radius: 1.0,
                },
            ],
        };
        let mut rate_rng = StdRng::seed_from_u64(23);
        assert!(
            try_sample_insertion_point(&disjoint, &mut rate_rng).is_none(),
            "SOIL rejection-budget exhaustion must become a rejected rate-insertion candidate"
        );
    }

    // ── SpatialHash tests ───────────────────────────────────────────────────

    fn no_pbc() -> PeriodicBox {
        PeriodicBox {
            is_periodic: [false; 3],
            box_size: [1.0; 3],
        }
    }

    #[test]
    fn spatial_hash_no_overlap() {
        let mut hash = SpatialHash::new(0.1);
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let radii = vec![0.01, 0.01];
        for (i, pos) in positions.iter().enumerate() {
            hash.insert(i, pos);
        }
        // Far away — no overlap
        assert!(!hash.has_overlap(&[1.0, 1.0, 1.0], 0.01, &positions, &radii, &no_pbc()));
    }

    #[test]
    fn spatial_hash_detects_overlap() {
        let mut hash = SpatialHash::new(0.1);
        let positions = vec![[0.0, 0.0, 0.0]];
        let radii = vec![0.05];
        hash.insert(0, &positions[0]);
        // Close enough to overlap
        assert!(hash.has_overlap(&[0.05, 0.0, 0.0], 0.05, &positions, &radii, &no_pbc()));
    }

    #[test]
    fn spatial_hash_near_boundary() {
        let mut hash = SpatialHash::new(0.1);
        let positions = vec![[0.09, 0.0, 0.0]];
        let radii = vec![0.04];
        hash.insert(0, &positions[0]);
        // Just across cell boundary — should still detect overlap
        assert!(hash.has_overlap(&[0.11, 0.0, 0.0], 0.04, &positions, &radii, &no_pbc()));
    }

    #[test]
    fn spatial_hash_periodic_overlap() {
        // Particles near opposite edges of a periodic box should overlap
        let pbc = PeriodicBox {
            is_periodic: [false, true, false],
            box_size: [1.0, 0.1, 1.0],
        };
        let mut hash = SpatialHash::new(0.05);
        let positions = vec![[0.0, 0.005, 0.0]]; // near y=0 edge
        let radii = vec![0.02];
        hash.insert(0, &positions[0]);
        // Near y=0.095 edge — through PBC, distance is 0.01 < 2*0.02
        assert!(hash.has_overlap(&[0.0, 0.095, 0.0], 0.02, &positions, &radii, &pbc));
    }

    // ── InsertConfig deserialization tests ───────────────────────────────────

    #[test]
    fn insert_config_backward_compat() {
        let toml_str = r#"
material = "glass"
count = 100
radius = 0.001
density = 2500.0
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.material.as_deref(), Some("glass"));
        assert_eq!(config.count, Some(100));
        assert_eq!(config.density, Some(2500.0));
        assert!(matches!(config.radius, Some(RadiusSpec::Fixed(r)) if (r - 0.001).abs() < 1e-15));
        assert!(config.rate.is_none());
        assert_eq!(config.source, "random");
    }

    #[test]
    fn insert_config_with_distribution() {
        let toml_str = r#"
material = "glass"
count = 500
density = 2500.0
radius = { distribution = "uniform", min = 0.0008, max = 0.0012 }
velocity_z = -1.0
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert!(matches!(
            config.radius,
            Some(RadiusSpec::Distribution(RadiusDistribution::Uniform { .. }))
        ));
        assert_eq!(config.velocity_z, Some(-1.0));
    }

    #[test]
    fn insert_config_rate_based() {
        let toml_str = r#"
material = "glass"
density = 2500.0
radius = { distribution = "uniform", min = 0.0008, max = 0.0012 }
velocity_z = -1.0
rate = 10
rate_interval = 100
rate_start = 0
rate_end = 500000
rate_limit = 5000
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.rate, Some(10));
        assert_eq!(config.rate_interval, Some(100));
        assert_eq!(config.rate_start, Some(0));
        assert_eq!(config.rate_end, Some(500000));
        assert_eq!(config.rate_limit, Some(5000));
    }

    #[test]
    fn rate_insert_missing_rate_reports_validation_error() {
        let toml_str = r#"
material = "glass"
density = 2500.0
radius = 0.001
rate_interval = 100
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert!(is_rate_insert_config(&config));
        let err = validate_rate_insert_config(&config, "Rate-based [[particles.insert]]")
            .expect_err("missing rate should be reported before runtime insertion");
        assert!(err.contains("requires 'rate'"));
    }

    #[test]
    fn rate_insert_missing_radius_reports_validation_error() {
        let toml_str = r#"
material = "glass"
density = 2500.0
rate = 10
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        let err = validate_rate_insert_config(&config, "Rate-based [[particles.insert]]")
            .expect_err("missing radius should be reported before runtime insertion");
        assert!(err.contains("requires 'radius'"));
    }

    #[test]
    fn rate_insert_missing_density_reports_validation_error() {
        let toml_str = r#"
material = "glass"
radius = 0.001
rate = 10
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        let err = validate_rate_insert_config(&config, "Rate-based [[particles.insert]]")
            .expect_err("missing density should be reported before runtime insertion");
        assert!(err.contains("requires 'density'"));
    }

    #[test]
    fn negative_random_velocity_reports_validation_error() {
        let err = validate_insert_velocity(-0.1, "[[particles.insert]]")
            .expect_err("negative random velocity should be rejected");
        assert!(err.contains("must be finite and non-negative"));
    }

    #[test]
    fn insert_config_file_based_csv() {
        let toml_str = r#"
source = "file"
file = "particles.csv"
format = "csv"
material = "glass"
density = 2500.0
columns = { x = 0, y = 1, z = 2, radius = 3 }
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.source, "file");
        assert_eq!(config.file.as_deref(), Some("particles.csv"));
        assert_eq!(config.format.as_deref(), Some("csv"));
        let cols = config.columns.unwrap();
        assert_eq!(cols.x, Some(0));
        assert_eq!(cols.radius, Some(3));
    }

    #[test]
    fn insert_config_file_based_lammps() {
        let toml_str = r#"
source = "file"
file = "dump.lammpstrj"
format = "lammps_dump"
material = "glass"
density = 2500.0
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.source, "file");
        assert_eq!(config.format.as_deref(), Some("lammps_dump"));
    }

    #[test]
    fn insert_config_with_type_map() {
        let toml_str = r#"
source = "file"
file = "dump.lammpstrj"
format = "lammps_dump"
material = "glass"
density = 2500.0
type_map = { 1 = "glass", 2 = "steel" }
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        let tm = config.type_map.unwrap();
        assert_eq!(tm.len(), 2);
        assert_eq!(tm["1"], "glass");
        assert_eq!(tm["2"], "steel");
        assert_eq!(config.material.as_deref(), Some("glass"));
    }

    #[test]
    fn insert_config_type_map_with_fallback() {
        let toml_str = r#"
source = "file"
file = "particles.csv"
format = "csv"
density = 2500.0
material = "glass"
type_map = { 2 = "steel" }
columns = { x = 0, y = 1, z = 2, radius = 3, atom_type = 4 }
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.material.as_deref(), Some("glass"));
        let tm = config.type_map.unwrap();
        assert_eq!(tm.len(), 1);
        assert_eq!(tm["2"], "steel");
    }

    #[test]
    fn insert_config_no_type_map_backward_compat() {
        let toml_str = r#"
source = "file"
file = "dump.lammpstrj"
format = "lammps_dump"
material = "glass"
density = 2500.0
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert!(config.type_map.is_none());
    }

    #[test]
    fn insert_config_lammps_data() {
        let toml_str = r#"
source = "file"
file = "data.lammps"
format = "lammps_data"
material = "glass"
density = 2500.0
radius = 0.001
type_map = { 1 = "glass", 2 = "steel" }
atom_style = "atomic"
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.source, "file");
        assert_eq!(config.format.as_deref(), Some("lammps_data"));
        assert_eq!(config.atom_style.as_deref(), Some("atomic"));
        let tm = config.type_map.unwrap();
        assert_eq!(tm.len(), 2);
    }

    #[test]
    fn insert_config_lammps_data_sphere_style() {
        let toml_str = r#"
source = "file"
file = "data.lammps"
format = "lammps_data"
material = "glass"
atom_style = "bpm/sphere"
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(config.atom_style.as_deref(), Some("bpm/sphere"));
        // No density/radius required for sphere style (per-atom in file)
        assert!(config.density.is_none());
        assert!(config.radius.is_none());
    }

    #[test]
    fn insert_config_with_region() {
        let toml_str = r#"
material = "glass"
count = 100
radius = 0.001
density = 2500.0
region = { type = "cylinder", center = [0.01, 0.01], radius = 0.008, axis = "z", lo = 0.04, hi = 0.05 }
"#;
        let config: InsertConfig = toml::from_str(toml_str).unwrap();
        assert!(config.region.is_some());
    }
}
