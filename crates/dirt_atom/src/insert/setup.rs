use super::*;

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
            dem_insert_atoms.after(soil_core::domain_read_input),
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
pub(super) fn validate_particles_config(
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
pub(super) fn validate_stage_particles_configs(
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

pub(super) fn validate_preflight_region(
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
pub(super) fn validate_file_insert(
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

// ── Delta time calculation ──────────────────────────────────────────────────

/// Computes a stable timestep from the Rayleigh wave speed criterion.
///
/// For each particle, estimates the Rayleigh wave transit time across the particle
/// diameter using `dt_R = π·r / α · √(ρ/G)`, where α ≈ 0.1631·ν + 0.8766 and
/// G = E / (2·(1+ν)). The final timestep is 15% of the minimum across all particles.
pub(super) fn calculate_delta_time(
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
