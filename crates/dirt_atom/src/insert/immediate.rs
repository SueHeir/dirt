use super::*;

/// Typed initialization record for a newly-created DEM SoA row. The store owns
/// structural synchronization; this record only supplies DIRT-specific values.
#[derive(Clone, Copy, Debug)]

pub(super) struct DemParticle {
    pub(super) pos: [f64; 3],
    pub(super) vel: [f64; 3],
    pub(super) radius: f64,
    pub(super) cutoff_padding: f64,
    pub(super) density: f64,
    pub(super) mat_idx: u32,
    pub(super) tag: u32,
}

impl DemParticle {
    pub(super) fn mass(self) -> f64 {
        self.density * 4.0 / 3.0 * PI * self.radius.powi(3)
    }

    /// Populate the already-reserved core row.  Structural changes belong to
    /// `ParticleStore`; this method deliberately only writes an existing row.
    pub(super) fn write_core(self, atom: &mut Atom, i: usize, mass: f64) {
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
    pub(super) fn write_dem(self, dem: &mut DemAtom, i: usize, mass: f64) {
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

pub(super) fn insert_single_particle(
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    row: DemParticle,
) {
    let global_natoms = atom
        .natoms()
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
pub(super) fn owns_position(domain: &Domain, pos: &[f64; 3]) -> bool {
    (0..3).all(|d| pos[d] >= domain.sub_domain_low[d] && pos[d] < domain.sub_domain_high[d])
}

// ── Helper: resolve material index ──────────────────────────────────────────

pub(super) fn resolve_material(material_table: &MaterialTable, name: &str) -> Result<u32, String> {
    material_table.find_material(name).ok_or_else(|| {
        format!(
            "unknown material '{}' in [[particles.insert]]. Available: {:?}",
            name, material_table.names
        )
    })
}

pub(super) fn resolve_file_material(
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
pub(super) fn resolve_type_map(
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
pub(super) fn lookup_material_for_type(
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
                    let prepared = prepare_random_insert(
                        insert,
                        &material_table,
                        &domain,
                        "rate-based [[particles.insert]]",
                    )
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
                        prepared,
                        total_inserted: 0,
                    });
                } else {
                    // ── Immediate random insertion (deterministic, born-in-owner) ──
                    let mat_name = insert.material.as_deref().expect(
                        "random insertion material was validated during fallible plugin preflight",
                    );
                    let count = insert.count.expect(
                        "random insertion count was validated during fallible plugin preflight",
                    );
                    let prepared = prepare_random_insert(
                        insert,
                        &material_table,
                        &domain,
                        "[[particles.insert]]",
                    )
                    .expect("random insertion was validated during fallible plugin preflight");
                    if comm.rank() == 0 {
                        println!(
                            "DemAtomInsert: inserting {} particles of material '{}' (r={}, rho={}, E={}, nu={})",
                            count,
                            mat_name,
                            prepared.max_radius,
                            prepared.density,
                            material_table.youngs_mod[prepared.mat_idx as usize],
                            material_table.poisson_ratio[prepared.mat_idx as usize]
                        );
                    }
                    let mut candidates =
                        CandidateGenerator::new(&prepared, &domain, prepared.seed, count as usize);

                    let mut inserted = 0u32;
                    let mut attempts = 0u64;
                    let max_attempts = count as u64 * 1_000_000;
                    while inserted < count && attempts < max_attempts {
                        attempts += 1;
                        let Some(candidate) = candidates.next(&prepared) else {
                            continue;
                        };
                        let tag = max_tag;
                        max_tag += 1;

                        // Materialize into the Atom arrays only if this rank owns the
                        // position (half-open interval, matching exchange() ownership).
                        if owns_position(&domain, &candidate.pos) {
                            insert_single_particle(
                                &mut atom,
                                &registry,
                                candidate.particle(&prepared, tag),
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

pub(super) fn insert_from_file(
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
