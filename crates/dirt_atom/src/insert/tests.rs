use super::*;
use crate::{Elastic, Friction, Material, RadiusDistribution};
use soil_core::{toml, AtomData, AtomDataRegistry, ParticleStoreError};
use soil_derive::AtomData;

fn add_test_glass(materials: &mut MaterialTable) {
    materials
        .add(
            Material::new("glass", Elastic::new(8.7e9, 0.3, 0.9)).with_friction(Friction {
                sliding: 0.5,
                ..Friction::default()
            }),
        )
        .unwrap();
}

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
    unsafe fn push_default(&mut self) {}
    unsafe fn truncate(&mut self, _: usize) {}
    unsafe fn swap_remove(&mut self, _: usize) {}
    fn pack(&self, _: usize, _: &mut Vec<f64>) {}
    unsafe fn unpack(&mut self, _: &[f64]) -> usize {
        0
    }
    unsafe fn apply_permutation(&mut self, _: &[usize], _: usize) {}
}

fn test_dem_registry() -> AtomDataRegistry {
    let mut registry = AtomDataRegistry::new();
    registry.register_for_atoms(DemAtom::new(), &Atom::new()).unwrap();
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
    config: InsertConfig,
) -> (Atom, usize, u32, CommState) {
    let mut materials = MaterialTable::new();
    add_test_glass(&mut materials);
    materials.build_pair_tables();
    let prepared = prepare_random_insert(
        &config,
        &materials,
        &domain,
        "rate-based [[particles.insert]]",
    )
    .unwrap();
    let mut app = App::new();
    app.add_resource(CommResource(Box::new(soil_core::SingleProcessComm::new())));
    app.add_resource(domain);
    app.add_resource(std::mem::take(&mut atom));
    app.add_resource(registry);
    app.add_resource(RunState::new());
    app.add_resource(materials);
    app.add_resource(RateInsertState {
        entries: vec![RateInsertEntry {
            config,
            prepared,
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
    let comm_state = app.get_resource_ref::<CurrentState<CommState>>().unwrap().0;
    (atom, late_rows, inserted, comm_state)
}

fn run_immediate_once(domain: Domain, config: &InsertConfig) -> Atom {
    // Exercise the production setup system through App scheduling, with the
    // same TOML-shaped StageOverrides used by a normal first run stage.
    let mut insert = toml::Table::new();
    insert.insert("source".into(), toml::Value::String("random".into()));
    insert.insert(
        "material".into(),
        toml::Value::String(config.material.clone().unwrap()),
    );
    insert.insert(
        "count".into(),
        toml::Value::Integer(config.count.unwrap().into()),
    );
    insert.insert(
        "radius".into(),
        toml::Value::Float(match config.radius.as_ref().unwrap() {
            RadiusSpec::Fixed(radius) => *radius,
            other => panic!("test requires a fixed radius, got {other:?}"),
        }),
    );
    insert.insert(
        "density".into(),
        toml::Value::Float(config.density.unwrap()),
    );
    insert.insert(
        "seed".into(),
        toml::Value::Integer(config.seed.unwrap().try_into().unwrap()),
    );
    let mut particles = toml::Table::new();
    particles.insert(
        "insert".into(),
        toml::Value::Array(vec![toml::Value::Table(insert)]),
    );
    let mut stage_table = toml::Table::new();
    stage_table.insert("particles".into(), toml::Value::Table(particles));

    let mut materials = MaterialTable::new();
    add_test_glass(&mut materials);
    materials.build_pair_tables();
    let mut app = App::new();
    app.add_resource(CommResource(Box::new(soil_core::SingleProcessComm::new())));
    app.add_resource(domain);
    app.add_resource(Atom::new());
    app.add_resource(test_dem_registry());
    app.add_resource(materials);
    app.add_resource(StageOverrides { table: stage_table });
    app.add_resource(RunConfig::default());
    app.add_resource(SchedulerManager::default());
    app.add_resource(RateInsertState::default());
    app.add_update_system(
        dem_insert_atoms,
        ParticleSimScheduleSet::PreInitialIntegration,
    );
    app.organize_systems();
    app.run();
    let atom = app.get_resource_ref::<Atom>().unwrap().clone();
    atom
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
fn immediate_and_rate_share_fixed_seed_candidate_and_tag_streams() {
    let mut materials = MaterialTable::new();
    add_test_glass(&mut materials);
    materials.build_pair_tables();
    let domain = unit_domain(0.0, 1.0);
    let mut config = rate_config(20260712);
    config.count = Some(6);
    let prepared = prepare_random_insert(&config, &materials, &domain, "test insertion")
        .expect("well-formed insertion should prepare once for both schedules");

    // A rate event at step zero / entry zero uses the configured seed, exactly
    // like immediate insertion. The shared generator must therefore produce
    // bit-identical accepted particles.
    let mut immediate = CandidateGenerator::new(&prepared, &domain, prepared.seed, 6);
    let mut rate = CandidateGenerator::new(&prepared, &domain, prepared.seed, 6);
    let immediate_rows: Vec<_> = (10..16)
        .map(|tag| (tag, immediate.next(&prepared).unwrap()))
        .collect();
    let rate_rows: Vec<_> = (10..16)
        .map(|tag| (tag, rate.next(&prepared).unwrap()))
        .collect();
    assert_eq!(immediate_rows, rate_rows);

    // Ownership filtering cannot affect the replicated stream: combining two
    // half-open partitions recovers each accepted tag exactly once.
    for (tag, candidate) in &immediate_rows {
        assert_eq!(
            owns_position(&unit_domain(0.0, 0.5), &candidate.pos) as u8
                + owns_position(&unit_domain(0.5, 1.0), &candidate.pos) as u8,
            1,
            "tag {tag} must have exactly one owner"
        );
    }
}

#[test]
fn prepared_insert_rejects_malformed_random_config_before_scheduling() {
    let mut materials = MaterialTable::new();
    add_test_glass(&mut materials);
    materials.build_pair_tables();
    let mut config = rate_config(1);
    config.density = None;
    let error = prepare_random_insert(&config, &materials, &unit_domain(0.0, 1.0), "test")
        .expect_err("prepared insertion must reject incomplete config");
    assert!(error.contains("requires 'density'"));
}

#[test]
fn production_immediate_and_rate_match_accepted_rows_tags_after_rejection() {
    // Find a fixed seed with an overlap rejection between two accepted
    // candidates. This specifically guards the former rate-only tag gap:
    // both actual scheduled systems must assign tags after acceptance.
    let domain = unit_domain(0.0, 1.0);
    let mut immediate_config = rate_config(0);
    immediate_config.count = Some(2);
    immediate_config.rate = None;
    immediate_config.rate_interval = None;
    immediate_config.rate_start = None;
    immediate_config.rate_end = None;
    immediate_config.rate_limit = None;
    immediate_config.radius = Some(RadiusSpec::Fixed(0.25));
    immediate_config.region = None;

    let mut materials = MaterialTable::new();
    add_test_glass(&mut materials);
    materials.build_pair_tables();
    let seed = (0..10_000u64)
        .find(|&seed| {
            immediate_config.seed = Some(seed);
            let prepared =
                prepare_random_insert(&immediate_config, &materials, &domain, "test insertion")
                    .unwrap();
            let mut candidates = CandidateGenerator::new(&prepared, &domain, seed, 2);
            candidates.next(&prepared).is_some()
                && candidates.next(&prepared).is_none()
                && candidates.next(&prepared).is_some()
        })
        .expect("a bounded fixed-radius stream should expose an overlap rejection");
    immediate_config.seed = Some(seed);

    let immediate = run_immediate_once(unit_domain(0.0, 1.0), &immediate_config);
    assert_eq!(immediate.tag, vec![0, 1]);

    let mut rate_config = immediate_config.clone();
    rate_config.count = None;
    rate_config.rate = Some(2);
    rate_config.rate_interval = Some(1);
    rate_config.rate_start = Some(0);
    rate_config.rate_end = Some(0);
    rate_config.rate_limit = Some(2);
    let (rate, _, _, serial_comm_state) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.0, 1.0),
        rate_config.clone(),
    );
    assert_eq!(rate.tag, vec![0, 1]);
    assert_eq!(serial_comm_state, CommState::FullRebuild);
    assert_eq!(
        immediate
            .tag
            .iter()
            .zip(immediate.pos.iter())
            .collect::<Vec<_>>(),
        rate.tag.iter().zip(rate.pos.iter()).collect::<Vec<_>>(),
        "the scheduled immediate and rate paths must retain the same accepted tag/row stream"
    );

    // The rate system's replicated acceptance stream must still partition
    // exactly once across two half-open ownership domains.
    let (low, _, _, low_comm_state) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.0, 0.5),
        rate_config.clone(),
    );
    let (high, _, _, high_comm_state) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.5, 1.0),
        rate_config,
    );
    assert_eq!(low_comm_state, serial_comm_state);
    assert_eq!(high_comm_state, serial_comm_state);
    let mut partitioned: Vec<_> = low
        .tag
        .iter()
        .zip(low.pos.iter())
        .chain(high.tag.iter().zip(high.pos.iter()))
        .map(|(tag, pos)| (*tag, *pos))
        .collect();
    let mut serial: Vec<_> = rate
        .tag
        .iter()
        .zip(rate.pos.iter())
        .map(|(tag, pos)| (*tag, *pos))
        .collect();
    partitioned.sort_by_key(|row| row.0);
    serial.sort_by_key(|row| row.0);
    assert_eq!(partitioned, serial);
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
    registry.register_for_atoms(LateProbe::default(), &atom).unwrap();

    let (atom, late_rows, inserted, _) =
        run_rate_once(atom, registry, unit_domain(0.0, 1.0), rate_config(20260712));
    assert_eq!(inserted, 8);
    assert_eq!((atom.nlocal(), atom.nghost(), atom.len()), (9, 0, 9));
    assert_eq!(atom.tag[0], 41, "the local prefix survives ghost removal");
    assert_eq!(late_rows, atom.len());

    let (_, _, full_inserted, _) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.0, 1.0),
        rate_config(99),
    );
    assert_eq!(full_inserted, 8);
    let (low, _, _, _) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.0, 0.5),
        rate_config(99),
    );
    let (high, _, _, _) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.5, 1.0),
        rate_config(99),
    );
    let (full, _, _, _) = run_rate_once(
        Atom::new(),
        test_dem_registry(),
        unit_domain(0.0, 1.0),
        rate_config(99),
    );
    let mut partitioned: Vec<_> = low
        .tag
        .iter()
        .zip(low.pos.iter())
        .chain(high.tag.iter().zip(high.pos.iter()))
        .map(|(tag, pos)| (*tag, *pos))
        .collect();
    let mut serial: Vec<_> = full
        .tag
        .iter()
        .zip(full.pos.iter())
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
    assert_eq!((atoms.nlocal(), atoms.nghost(), atoms.natoms()), (2, 0, 2));
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
    registry.register_for_atoms(LateProbe::default(), &atoms).unwrap();
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
    rollback_registry
        .register_for_atoms(BrokenDefaults, &rollback_atoms)
        .unwrap();
    assert_eq!(
        ParticleStore::new(&mut rollback_atoms, &rollback_registry).push_default_local(1),
        Err(ParticleStoreError::MalformedExtensionRecord)
    );
    assert!(rollback_atoms.is_empty());
    assert_eq!((rollback_atoms.nlocal(), rollback_atoms.natoms()), (0, 0));
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

    // Structural columns are intentionally no longer directly clearable. A
    // mismatched ownership count is the same malformed-restart class and must
    // leave the synchronized core and extension stores untouched.
    // The SOIL facade intentionally does not permit constructing an Atom with
    // inconsistent ownership counters.  An empty restart is the public way to
    // exercise the same rejected structural replacement boundary.
    let malformed = Atom::new();
    assert_eq!(
        ParticleStore::new(&mut atoms, &registry).replace_from_restart(malformed, &[]),
        Err(ParticleStoreError::MalformedExtensionRecord)
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
    assert_eq!((atoms.nlocal(), atoms.nghost()), (1, 1));
    assert_eq!(
        registry.expect::<DemAtom>("ghost setup").radius,
        vec![0.001, 0.003]
    );

    ParticleStore::new(&mut atoms, &registry)
        .discard_ghosts()
        .unwrap();
    assert_eq!((atoms.nlocal(), atoms.nghost(), atoms.len()), (1, 0, 1));
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
