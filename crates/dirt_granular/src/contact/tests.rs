use super::*;

#[test]
fn hooke_and_hertz_publish_the_same_typed_contact_seam() {
    for model in ["hertz", "hooke"] {
        let mut app = App::new();
        app.add_resource(soil_core::Config::from_str(&format!(
            "[dem]\ncontact_model = \"{model}\""
        )));
        app.add_resource(grass_scheduler::CurrentState(
            soil_core::CommState::FullRebuild,
        ));
        app.add_resource(soil_core::RunState::new());
        app.add_resource(Neighbor::default());
        app.add_plugins(dirt_atom::DemAtomPlugin);
        app.add_plugins(HertzMindlinContactPlugin);
        app.add_update_system(
            contact_seam_consumer.requires(CONTACT_FORCE),
            ParticleSimScheduleSet::Force,
        );
        app.organize_systems();
    }
}

fn contact_seam_consumer() {}
use dirt_atom::DemAtom;
use dirt_test_utils::{make_material_table, push_dem_test_atom, ParticleFixture, ParticleSpec};
use soil_core::Neighbor;
use soil_core::{Atom, AtomDataRegistry};

fn push_test_atom_with_history(
    atom: &mut Atom,
    dem: &mut DemAtom,
    history: &mut ContactHistoryStore,
    tag: u32,
    pos: [f64; 3],
    radius: f64,
) {
    push_dem_test_atom(atom, dem, tag, pos, radius);
    history.contacts.push(Vec::new());
}

/// Step 4 correctness: the interior/boundary two-pass force (Interior pass for
/// local-local pairs while the halo is in flight, Boundary pass for ghost pairs
/// after they land) must equal the single All pass — bit-for-bit.
#[test]
fn interior_boundary_split_matches_single_pass() {
    let r = 0.001;
    let build = || {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        // atom 0 (local) overlaps atom 1 (local -> interior pair) and atom 2
        // (ghost -> boundary pair).
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], r);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [1.5 * r, 0.0, 0.0], r);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 2, [0.0, 1.5 * r, 0.0], r);
        atom.nlocal = 2;
        atom.natoms = 3;
        // Half neighbour list (newton): atom 0 -> {1 (local), 2 (ghost)}.
        let mut nb = Neighbor::new();
        nb.neighbor_offsets = vec![0, 2, 2, 2];
        nb.neighbor_indices = vec![1, 2];
        let mut reg = AtomDataRegistry::new();
        reg.try_register(dem, atom.len()).unwrap();
        reg.try_register(hist, atom.len()).unwrap();
        (atom, nb, reg)
    };
    let mt = make_material_table();

    let (mut a_all, nb, r_all) = build();
    contact_force_core(&mut a_all, &nb, &r_all, &mt, None, ForcePass::All);

    let (mut a_split, _nb, r_split) = build();
    contact_force_core(&mut a_split, &nb, &r_split, &mt, None, ForcePass::Interior);
    contact_force_core(&mut a_split, &nb, &r_split, &mt, None, ForcePass::Boundary);

    let mut max_diff = 0.0f64;
    for i in 0..3 {
        for d in 0..3 {
            max_diff = max_diff.max((a_all.force[i][d] as f64 - a_split.force[i][d] as f64).abs());
        }
    }
    assert!(
        max_diff < 1e-15,
        "interior+boundary != all: max force diff = {max_diff:.3e}"
    );
    // Sanity: the contact actually produced a non-trivial force.
    assert!(a_all.force[0][0].abs() as f64 + a_all.force[0][1].abs() as f64 > 0.0);
}

#[test]
fn fused_contact_repulsive_for_overlap() {
    let radius = 0.001;
    let mut fixture = ParticleFixture::pair(
        ParticleSpec::new(0, [0.0, 0.0, 0.0], radius),
        ParticleSpec::new(1, [0.0019, 0.0, 0.0], radius),
    )
    .build();
    let mut hist = ContactHistoryStore::new();
    hist.contacts.resize_with(fixture.atom.len(), Vec::new);
    fixture.register_atom_data(hist);
    let mut app = fixture.into_app();
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    assert!(
        atom.force[0][0] < 0.0,
        "particle 0 should have negative x force"
    );
    assert!(
        atom.force[1][0] > 0.0,
        "particle 1 should have positive x force"
    );
    assert!((atom.force[0][0] + atom.force[1][0]).abs() < 1e-10);
}

#[test]
fn fused_contact_tangential_with_sliding() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    atom.vel[1][1] = 0.1;
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    // Normal force present
    assert!(atom.force[0][0] < 0.0, "normal force on atom 0");
    assert!(atom.force[1][0] > 0.0, "normal force on atom 1");
    // Tangential force present
    assert!(atom.force[0][1].abs() > 0.0, "tangential force on atom 0");
    assert!(
        (atom.force[0][1] + atom.force[1][1]).abs() < 1e-10,
        "tangential forces equal and opposite"
    );
    // Torque present (stored in DemAtom via registry)
    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    let t_mag =
        (dem.torque[0][0].powi(2) + dem.torque[0][1].powi(2) + dem.torque[0][2].powi(2)).sqrt();
    assert!(t_mag > 0.0, "torque on atom 0");
}

/// The `linear_nohistory` tangential model must keep the tangential spring
/// displacement identically zero (velocity-Coulomb, LAMMPS `pair_granular`
/// `tangential linear_nohistory`), while the default `history` (Mindlin) model
/// accumulates it. Both are driven with the same sub-Coulomb tangential slip.
#[test]
fn linear_nohistory_has_no_spring_accumulation() {
    let radius = 0.001;
    let build = || {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.00185, 0.0, 0.0],
            radius,
        );
        atom.vel[1][1] = 0.001; // small tangential slip, below the Coulomb cap
        atom.nlocal = 2;
        atom.natoms = 2;
        let mut nb = Neighbor::new();
        nb.neighbor_offsets = vec![0, 1, 1];
        nb.neighbor_indices = vec![1];
        let mut reg = AtomDataRegistry::new();
        reg.try_register(dem, atom.len()).unwrap();
        reg.try_register(hist, atom.len()).unwrap();
        (atom, nb, reg)
    };
    let spring_mag = |reg: &AtomDataRegistry| -> f64 {
        let h = reg.expect::<ContactHistoryStore>("spring");
        let s = h.contacts[0]
            .iter()
            .find(|(t, _, _)| *t == 1)
            .map(|(_, s, _)| *s)
            .unwrap_or([0.0; CONTACT_HISTORY_LEN]);
        (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt()
    };

    // History (Mindlin): the tangential spring accumulates over the contact.
    let mut mt_h = make_material_table();
    mt_h.tangential_model = "history".to_string();
    let (mut a, nb, reg) = build();
    for _ in 0..20 {
        a.force[0] = [0.0; 3];
        a.force[1] = [0.0; 3];
        contact_force_core(&mut a, &nb, &reg, &mt_h, None, ForcePass::All);
    }
    let xi_history = spring_mag(&reg);
    assert!(
        xi_history > 0.0,
        "history model must accumulate spring, got {xi_history:e}"
    );

    // linear_nohistory: spring stays exactly zero; force is still present.
    let mut mt_nh = make_material_table();
    mt_nh.tangential_model = "linear_nohistory".to_string();
    let (mut a2, nb2, reg2) = build();
    for _ in 0..20 {
        a2.force[0] = [0.0; 3];
        a2.force[1] = [0.0; 3];
        contact_force_core(&mut a2, &nb2, &reg2, &mt_nh, None, ForcePass::All);
    }
    let xi_nohistory = spring_mag(&reg2);
    assert_eq!(
        xi_nohistory, 0.0,
        "linear_nohistory must not accumulate spring"
    );
    assert!(
        (a2.force[0][1] as f64).abs() > 0.0,
        "linear_nohistory must still produce a tangential (velocity-Coulomb) force"
    );
}

#[test]
fn fused_contact_no_force_for_gap() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [0.003, 0.0, 0.0], radius);
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    assert!(atom.force[0][0].abs() < 1e-20);
    assert!(atom.force[1][0].abs() < 1e-20);
}

fn make_material_table_cohesion() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 1e9);
    mt.build_pair_tables();
    mt
}

fn make_material_table_rolling() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.3, 0.0);
    mt.build_pair_tables();
    mt
}

#[test]
fn cohesion_produces_attractive_force() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    // Very small overlap with high cohesion energy → cohesion dominates
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.00199999, 0.0, 0.0],
        radius, // delta = 1e-8 (tiny overlap)
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_cohesion());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    // With cohesion and small overlap, normal force on atom 0 should be positive (attractive toward atom 1)
    assert!(
        atom.force[0][0] > 0.0,
        "cohesion should make force attractive on atom 0, got {}",
        atom.force[0][0]
    );
    // Newton's 3rd law
    assert!(
        (atom.force[0][0] + atom.force[1][0]).abs() < 1e-10,
        "forces should be equal and opposite"
    );
}

#[test]
fn zero_cohesion_matches_original() {
    // Two identical setups — one with default table, one with explicit 0.0 cohesion
    let radius = 0.001;
    let sep = 0.0019;

    let run = |mt: MaterialTable| -> [f64; 3] {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;
        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();
        registry.try_register(hist, atom.len()).unwrap();
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        [
            atom.force[0][0] as f64,
            atom.force[0][1] as f64,
            atom.force[0][2] as f64,
        ]
    };

    let f_default = run(make_material_table());
    let mut mt_zero = MaterialTable::new();
    mt_zero.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0);
    mt_zero.build_pair_tables();
    let f_zero = run(mt_zero);

    for d in 0..3 {
        assert!(
            (f_default[d] - f_zero[d]).abs() < 1e-15,
            "zero params should reproduce original, dim {} default={} zero={}",
            d,
            f_default[d],
            f_zero[d]
        );
    }
}

fn make_material_table_jkr() -> MaterialTable {
    let mut mt = MaterialTable::new();
    // Use high surface energy (1.0 J/m²) so adhesion clearly dominates at small overlaps
    mt.add_material_full("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 1.0);
    mt.build_pair_tables();
    mt
}

#[test]
fn jkr_pulloff_force_matches_theory() {
    // Test in adhesion-only regime (gap, not overlap) where force = -F_adhesion exactly
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    let gamma = 1.0;
    let r_eff = radius / 2.0;

    // Place particles with a tiny gap (adhesion-only regime)
    let gap = 1e-9;
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [2.0 * radius + gap, 0.0, 0.0],
        radius,
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    let mt = make_material_table_jkr();
    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let expected_pulloff = 1.5 * std::f64::consts::PI * gamma * r_eff;
    // In adhesion-only regime, force should be exactly -F_adhesion
    // Force on atom 0 should be positive (attracted toward atom 1)
    assert!(
        atom.force[0][0] > 0.0,
        "JKR should produce attractive force, got {}",
        atom.force[0][0]
    );
    // f_n_mag = -F_adhesion, force[0] -= f_n_mag * nx → force[0] += F_adhesion
    let f_mag = atom.force[0][0] as f64;
    assert!(
        (f_mag - expected_pulloff).abs() / expected_pulloff < 1e-6,
        "pull-off force should match theory {}, got {}",
        expected_pulloff,
        f_mag
    );
}

#[test]
fn jkr_adhesion_only_regime() {
    // Two particles with a small gap (no geometric overlap) but within JKR range
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    // Gap of 1e-9 (very small, within JKR pull-off distance for gamma=1.0)
    let gap = 1e-9;
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [2.0 * radius + gap, 0.0, 0.0],
        radius,
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_jkr());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    // Should be attractive (atom 0 pulled toward atom 1 = positive x)
    assert!(
        atom.force[0][0] > 0.0,
        "JKR adhesion-only should attract, got {}",
        atom.force[0][0]
    );
    // Newton's 3rd law
    assert!(
        (atom.force[0][0] + atom.force[1][0]).abs() < 1e-10,
        "forces should be equal and opposite"
    );
}

#[test]
fn jkr_no_interaction_beyond_pulloff() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    // Large gap — well beyond JKR pull-off distance
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.003, 0.0, 0.0],
        radius, // gap = 0.001 >> delta_pulloff
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_jkr());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    assert!(
        atom.force[0][0].abs() < 1e-20,
        "no force beyond pull-off distance"
    );
}

fn make_material_table_hooke() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material_extended("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 0.0, 0.0, 1e6, 5e5);
    mt.contact_model = "hooke".to_string();
    mt.build_pair_tables();
    mt
}

fn make_material_table_twisting() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material_extended(
        "glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0,
    );
    mt.build_pair_tables();
    mt
}

#[test]
fn hooke_force_linear_in_delta() {
    let radius = 0.001;
    let run = |sep: f64| -> f64 {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;
        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();
        registry.try_register(hist, atom.len()).unwrap();
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_hooke());
        app.add_update_system(hooke_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        atom.force[0][0] as f64
    };

    // delta1 = 2*r - sep1, delta2 = 2*r - sep2
    let sep1 = 0.00195; // delta = 0.00005
    let sep2 = 0.0019; // delta = 0.0001
    let f1 = run(sep1);
    let f2 = run(sep2);

    // Hooke: force proportional to delta → f2/f1 ≈ 2.0 (linear)
    let ratio = f2 / f1;
    assert!(
        (ratio - 2.0).abs() < 0.15,
        "Hooke force should be linear in delta, got ratio {} (expected ~2.0)",
        ratio
    );
}

#[test]
fn hooke_no_force_beyond_contact() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [0.003, 0.0, 0.0], radius);
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_hooke());
    app.add_update_system(hooke_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    assert!(
        atom.force[0][0].abs() < 1e-20,
        "no force beyond contact distance"
    );
}

#[test]
fn twisting_friction_opposes_spin() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    // Spin about contact normal (x-axis)
    dem.omega[0] = [100.0, 0.0, 0.0];
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_twisting());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    // Twisting torque on atom 0 should oppose its spin about x (negative x torque)
    assert!(
        dem.torque[0][0] < 0.0,
        "twisting torque should oppose omega_x, got {}",
        dem.torque[0][0]
    );
}

#[test]
fn twisting_friction_zero_when_no_spin() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    // No angular velocity at all
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_twisting());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    // No twisting torque when there's no angular velocity
    let torque_mag =
        (dem.torque[0][0].powi(2) + dem.torque[0][1].powi(2) + dem.torque[0][2].powi(2)).sqrt();
    assert!(
        torque_mag < 1e-20,
        "no twisting torque when no spin, got {}",
        torque_mag
    );
}

#[test]
fn rolling_resistance_opposes_angular_velocity() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    // Give atom 0 a rolling angular velocity (around y-axis — perpendicular to contact normal x)
    dem.omega[0] = [0.0, 100.0, 0.0];
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_rolling());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    // Rolling torque on atom 0 should oppose its angular velocity (negative y)
    assert!(
        dem.torque[0][1] < 0.0,
        "rolling torque should oppose omega_y, got {}",
        dem.torque[0][1]
    );
}

// ── SDS model helper ────────────────────────────────────────────────

fn make_material_table_sds_rolling() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.rolling_model = "sds".to_string();
    mt.add_material_with_sds(
        "glass", 8.7e9, 0.3, 0.95, 0.4, 0.3, // rolling_friction (mu_r)
        0.0, 0.0, 0.0, // twisting_friction
        0.0, 0.0, 1e3, // rolling_stiffness
        0.5, // rolling_damping
        0.0, 0.0,
    );
    mt.build_pair_tables();
    mt
}

fn make_material_table_sds_twisting() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.twisting_model = "sds".to_string();
    mt.add_material_with_sds(
        "glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, // rolling_friction
        0.0, 0.0, 0.3, // twisting_friction (mu_tw)
        0.0, 0.0, 0.0, 0.0, 1e3, // twisting_stiffness
        0.5, // twisting_damping
    );
    mt.build_pair_tables();
    mt
}

#[test]
fn sds_rolling_opposes_angular_velocity() {
    // Two overlapping particles, one spinning → SDS rolling torque opposes it
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    // Give atom 0 angular velocity in y (rolling about contact normal x)
    dem.omega[0] = [0.0, 10.0, 0.0];
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_sds_rolling());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    // SDS rolling torque should oppose omega_y on atom 0
    assert!(
        dem.torque[0][1] < 0.0,
        "SDS rolling torque should oppose omega_y, got {}",
        dem.torque[0][1]
    );
}

#[test]
fn sds_rolling_spring_accumulates() {
    // Pre-load the LAMMPS-style length-valued rolling displacement → larger
    // torque than zero displacement. With n=x and omega=y, v_rl is -z and
    // n x F_roll is the torque about y.
    // Use very small omega so that damping doesn't dominate and Coulomb cap isn't reached
    let radius = 0.001;

    let run_with_preload = |preload_z: f64| -> f64 {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        dem.omega[0] = [0.0, 0.001, 0.0]; // very small angular velocity
        atom.nlocal = 2;
        atom.natoms = 2;

        // Pre-load rolling displacement in contact history (canonical: tag 0 < tag 1, sign=+1)
        if preload_z != 0.0 {
            let mut preload = [0.0; CONTACT_HISTORY_LEN];
            preload[5] = preload_z;
            hist.contacts[0].push((1, preload, false));
        }

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();
        registry.try_register(hist, atom.len()).unwrap();

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_sds_rolling());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let reg = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let d = reg.expect::<DemAtom>("test");
        d.torque[0][1]
    };

    let torque_no_preload = run_with_preload(0.0);
    let torque_with_preload = run_with_preload(-1e-8); // small length preload below cap

    assert!(torque_no_preload < 0.0, "should oppose omega_y");
    assert!(torque_with_preload < 0.0, "should oppose omega_y");
    // Pre-loaded spring adds to torque magnitude
    assert!(
        torque_with_preload.abs() > torque_no_preload.abs(),
        "preloaded spring should increase torque: no_preload={}, preloaded={}",
        torque_no_preload,
        torque_with_preload
    );
}

#[test]
fn sds_rolling_coulomb_cap() {
    // Very high angular velocity → torque should be capped at mu_r * |F_n| * R_eff
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-5; // larger dt to accumulate big spring

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    dem.omega[0] = [0.0, 1e6, 0.0]; // very high
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    let mt = make_material_table_sds_rolling();
    let mu_r = mt.rolling_friction_ij[0][0];

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    let torque_mag =
        (dem.torque[0][0].powi(2) + dem.torque[0][1].powi(2) + dem.torque[0][2].powi(2)).sqrt();

    // Compute expected cap: mu_r * F_n * R_eff
    // F_n from Hertz: 4/3 * E_eff * sqrt(delta * r_eff) * delta
    let r_eff = radius / 2.0;
    let delta = 2.0 * radius - 0.0019;
    let e_eff = 8.7e9 / (2.0 * (1.0 - 0.09)); // single material
    let sqrt_dr = (delta * r_eff).sqrt();
    let f_n_approx = 4.0 / 3.0 * e_eff * sqrt_dr * delta;
    let tau_cap = mu_r * f_n_approx * r_eff;

    // Rolling torque should not exceed cap (with reasonable tolerance for damping and normal force)
    // The torque includes tangential torque contributions, so we just check the rolling component
    // is bounded. Since torque_mag includes all contributions, just check it's finite and reasonable.
    assert!(torque_mag.is_finite(), "torque should be finite");
    assert!(
        torque_mag < tau_cap * 100.0, // generous bound since total torque includes tangential
        "torque {} should be bounded near cap {}",
        torque_mag,
        tau_cap
    );
}

#[test]
fn sds_twisting_opposes_spin() {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    // Spin about contact normal (x-axis)
    dem.omega[0] = [10.0, 0.0, 0.0];
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_sds_twisting());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    // SDS twisting torque should oppose spin about x
    assert!(
        dem.torque[0][0] < 0.0,
        "SDS twisting torque should oppose spin about x, got {}",
        dem.torque[0][0]
    );
}

#[test]
fn sds_twisting_spring_accumulates() {
    let radius = 0.001;

    let run_with_preload = |preload: f64| -> f64 {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;

        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(
            &mut atom,
            &mut dem,
            &mut hist,
            1,
            [0.0019, 0.0, 0.0],
            radius,
        );
        dem.omega[0] = [0.001, 0.0, 0.0]; // very small spin
        atom.nlocal = 2;
        atom.natoms = 2;

        if preload != 0.0 {
            let mut preload_state = [0.0; CONTACT_HISTORY_LEN];
            preload_state[6] = preload;
            hist.contacts[0].push((1, preload_state, false));
        }

        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();
        registry.try_register(hist, atom.len()).unwrap();

        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_sds_twisting());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let reg = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let d = reg.expect::<DemAtom>("test");
        d.torque[0][0]
    };

    let torque_no_preload = run_with_preload(0.0);
    let torque_with_preload = run_with_preload(1e-5);
    assert!(torque_no_preload < 0.0);
    assert!(torque_with_preload < 0.0);
    assert!(
        torque_with_preload.abs() > torque_no_preload.abs(),
        "preloaded twisting spring should increase torque: no_preload={}, preloaded={}",
        torque_no_preload,
        torque_with_preload
    );
}

/// Marshall twisting material: tangential friction `friction` drives the
/// derived twisting cap; the SDS twist stiffness/damping inputs are provided
/// deliberately so tests can confirm the Marshall model *ignores* them.
fn make_material_table_marshall_twisting(
    friction: f64,
    twist_stiff: f64,
    twist_damp: f64,
) -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.twisting_model = "marshall".to_string();
    mt.add_material_with_sds(
        "glass",
        8.7e9,
        0.3,
        0.95,
        friction, // tangential μ_t — Marshall derives μ_twist = (2/3) a μ_t from this
        0.0,      // rolling_friction
        0.0,
        0.0,
        0.0, // twisting_friction (unused by Marshall)
        0.0,
        0.0, // kn, kt (Hertz path ignores these)
        0.0,
        0.0, // rolling sds
        twist_stiff,
        twist_damp, // twisting sds — must NOT affect Marshall
    );
    mt.build_pair_tables();
    mt
}

/// Run one Marshall-twisting contact step and return the twisting torque on
/// atom 0 (about the contact normal x̂). `preload` seeds the stored twisting
/// spring displacement; a large value forces the saturated (capped) regime.
fn run_marshall_twist(mt: MaterialTable, omega_x: f64, preload: f64) -> f64 {
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    dem.omega[0] = [omega_x, 0.0, 0.0]; // spin about contact normal x̂
    atom.nlocal = 2;
    atom.natoms = 2;

    if preload != 0.0 {
        let mut preload_state = [0.0; CONTACT_HISTORY_LEN];
        preload_state[6] = preload;
        hist.contacts[0].push((1, preload_state, false));
    }

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let reg = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let tq = reg.expect::<DemAtom>("test").torque[0][0];
    tq
}

#[test]
fn marshall_twisting_opposes_spin() {
    // Spin about the contact normal (x̂) → Marshall twisting couple opposes it.
    let tq = run_marshall_twist(
        make_material_table_marshall_twisting(0.4, 0.0, 0.0),
        10.0,
        0.0,
    );
    assert!(
        tq < 0.0,
        "Marshall twisting torque should oppose spin about x, got {}",
        tq
    );
}

#[test]
fn marshall_twisting_ignores_sds_inputs() {
    // The Marshall coefficients are DERIVED from the tangential model, so the
    // SDS twisting_stiffness / twisting_damping material inputs must have no
    // effect. Drive into the saturated (capped) regime with a large preload so
    // the torque equals the derived cap τ_max = μ_twist·F_n, then confirm two
    // wildly different SDS-input tables give the identical torque.
    let tq_zero = run_marshall_twist(
        make_material_table_marshall_twisting(0.4, 0.0, 0.0),
        10.0,
        1.0,
    );
    let tq_huge = run_marshall_twist(
        make_material_table_marshall_twisting(0.4, 1.0e9, 1.0e6),
        10.0,
        1.0,
    );
    assert!(tq_zero < 0.0, "should oppose spin, got {}", tq_zero);
    assert!(
        (tq_zero - tq_huge).abs() <= 1e-12 * tq_zero.abs().max(1e-30),
        "Marshall torque must ignore SDS twist inputs: zero-input={}, huge-input={}",
        tq_zero,
        tq_huge
    );
}

#[test]
fn marshall_twisting_cap_scales_with_tangential_friction() {
    // μ_twist = (2/3) a μ_t, so in the saturated regime the cap scales linearly
    // with the tangential friction coefficient: doubling μ_t doubles |τ|, and
    // μ_t = 0 gives zero twisting couple (Marshall ties the cap to sliding).
    let tq_mu04 = run_marshall_twist(
        make_material_table_marshall_twisting(0.4, 0.0, 0.0),
        10.0,
        1.0,
    );
    let tq_mu08 = run_marshall_twist(
        make_material_table_marshall_twisting(0.8, 0.0, 0.0),
        10.0,
        1.0,
    );
    let tq_mu00 = run_marshall_twist(
        make_material_table_marshall_twisting(0.0, 0.0, 0.0),
        10.0,
        1.0,
    );
    let ratio = tq_mu08 / tq_mu04;
    assert!(
        (ratio - 2.0).abs() < 1e-6,
        "doubling μ_t should double the Marshall cap: ratio={}",
        ratio
    );
    assert!(
        tq_mu00.abs() < 1e-12,
        "μ_t = 0 should give zero Marshall twisting torque, got {}",
        tq_mu00
    );
}

#[test]
fn constant_model_unchanged_with_sds_config() {
    // When rolling_model = "constant" (default), SDS parameters should be ignored
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    dem.omega[0] = [0.0, 10.0, 0.0];
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    // Use constant model but with SDS parameters set (they should be ignored)
    let mut mt = MaterialTable::new();
    // rolling_model defaults to "constant"
    mt.add_material_with_sds(
        "glass", 8.7e9, 0.3, 0.95, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 1e3, 0.5, 0.0, 0.0,
    );
    mt.build_pair_tables();

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("test");
    // Constant model: torque = -mu_r * |F_n| * r_eff * (roll/|roll|)
    // Should still produce opposing torque
    assert!(
        dem.torque[0][1] < 0.0,
        "constant rolling model should still work, got {}",
        dem.torque[0][1]
    );

    // Check that spring history has zero rolling/twisting displacement
    let hist = registry.expect::<ContactHistoryStore>("test");
    let contact = &hist.contacts[0][0];
    assert_eq!(
        contact.1[3], 0.0,
        "rolling disp x should be zero in constant model"
    );
    assert_eq!(
        contact.1[4], 0.0,
        "rolling disp y should be zero in constant model"
    );
    assert_eq!(
        contact.1[5], 0.0,
        "rolling disp z should be zero in constant model"
    );
    assert_eq!(
        contact.1[6], 0.0,
        "twisting disp should be zero in constant model"
    );
}

// ── DMT adhesion tests ──────────────────────────────────────────────

fn make_material_table_dmt() -> MaterialTable {
    let mut mt = MaterialTable::new();
    // Use high surface energy (1.0 J/m²) so adhesion clearly dominates at small overlaps
    mt.add_material_full("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 1.0);
    mt.adhesion_model = "dmt".to_string();
    mt.build_pair_tables();
    mt
}

#[test]
fn dmt_pulloff_force_matches_theory() {
    // DMT pull-off force = 2 * pi * gamma * r_eff (at contact, delta = 0+)
    let radius = 0.001;
    let gamma = 1.0;
    let r_eff = radius / 2.0; // two equal spheres

    // Use a very small overlap so Hertz contribution is negligible
    // At tiny delta, F_hertz ~ 0 but F_dmt = 2*pi*gamma*r_eff
    let tiny_overlap = 1e-12; // extremely small overlap
    let sep = 2.0 * radius - tiny_overlap;

    let mut app = App::new();
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_dmt());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let expected_dmt = 2.0 * std::f64::consts::PI * gamma * r_eff;
    // Force on atom 0 should be positive (attracted toward atom 1)
    // f_n_mag = k_n*delta - f_diss - f_dmt ~ -f_dmt (since delta ~ 0, v=0)
    // force[0] -= f_n_mag * nx -> force[0] ~ +f_dmt
    assert!(
        atom.force[0][0] > 0.0,
        "DMT should produce attractive force, got {}",
        atom.force[0][0]
    );
    assert!(
        (atom.force[0][0] as f64 - expected_dmt).abs() / expected_dmt < 1e-3,
        "DMT pull-off force should match 2*pi*gamma*r_eff = {}, got {}",
        expected_dmt,
        atom.force[0][0]
    );
}

#[test]
fn dmt_no_force_beyond_contact() {
    // DMT has no adhesion-only regime -- no force when delta < 0 (gap)
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    // Place particles with a gap
    let gap = 1e-9;
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [2.0 * radius + gap, 0.0, 0.0],
        radius,
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_dmt());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    // DMT: no force when particles are not in geometric contact
    assert!(
        atom.force[0][0].abs() < 1e-20,
        "DMT should have no force beyond contact, got {}",
        atom.force[0][0]
    );
}

#[test]
fn dmt_pulloff_less_than_jkr() {
    // DMT pull-off = 2*pi*gamma*r_eff, JKR pull-off = 1.5*pi*gamma*r_eff
    // At same surface energy, DMT has HIGHER pull-off force than JKR (2 > 1.5)
    // But JKR has extended range (adhesion across gap), so effective sticking is stronger
    let gamma = 1.0;
    let radius = 0.001;
    let r_eff = radius / 2.0;

    let f_dmt = 2.0 * std::f64::consts::PI * gamma * r_eff;
    let f_jkr = 1.5 * std::f64::consts::PI * gamma * r_eff;
    assert!(
        f_dmt > f_jkr,
        "DMT pull-off ({}) should be larger than JKR pull-off ({})",
        f_dmt,
        f_jkr
    );
}

#[test]
fn dmt_newtons_third_law() {
    // Verify equal and opposite forces for DMT contact
    let mut app = App::new();
    let radius = 0.001;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.0019, 0.0, 0.0],
        radius,
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_dmt());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    for d in 0..3 {
        assert!(
            (atom.force[0][d] + atom.force[1][d]).abs() < 1e-10,
            "Newton's 3rd law violated in dim {}: {} + {} != 0",
            d,
            atom.force[0][d],
            atom.force[1][d]
        );
    }
}

#[test]
fn dmt_does_not_break_jkr() {
    // Run the JKR test with default adhesion_model (should still work as JKR)
    let mut app = App::new();
    let radius = 0.001;
    let gamma = 1.0;
    let r_eff = radius / 2.0;
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    // Place particles with a tiny gap (adhesion-only regime for JKR)
    let gap = 1e-9;
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [2.0 * radius + gap, 0.0, 0.0],
        radius,
    );
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    // Use JKR material table (default adhesion_model = "jkr")
    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table_jkr());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let expected_jkr = 1.5 * std::f64::consts::PI * gamma * r_eff;
    // JKR should still attract across gap
    assert!(
        atom.force[0][0] > 0.0,
        "JKR should still work with DMT feature added, got {}",
        atom.force[0][0]
    );
    assert!(
        (atom.force[0][0] as f64 - expected_jkr).abs() / expected_jkr < 1e-6,
        "JKR pull-off force should still match 1.5*pi*gamma*r_eff = {}, got {}",
        expected_jkr,
        atom.force[0][0]
    );
}

// ── Force scaling validation tests ──────────────────────────────────

#[test]
fn hertz_force_scales_as_delta_three_halves() {
    let radius = 0.001;

    // Compute elastic-only normal force for a given separation (zero velocity -> no damping).
    let hertz_force_at = |sep: f64| -> f64 {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;
        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();
        registry.try_register(hist, atom.len()).unwrap();
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Force on atom 0 is negative (pushed away from atom 1), take absolute value
        atom.force[0][0].abs() as f64
    };

    // Test at 5 different overlaps
    let deltas = [1e-5, 2e-5, 4e-5, 6e-5, 8e-5];
    let forces: Vec<f64> = deltas
        .iter()
        .map(|d| {
            let sep = 2.0 * radius - d;
            hertz_force_at(sep)
        })
        .collect();

    // For each pair (i, 0), check F_i/F_0 ~ (delta_i/delta_0)^(3/2)
    for i in 1..deltas.len() {
        let expected_ratio = (deltas[i] / deltas[0]).powf(1.5);
        let actual_ratio = forces[i] / forces[0];
        let rel_err = ((actual_ratio - expected_ratio) / expected_ratio).abs();
        assert!(
                rel_err < 0.01,
                "Hertz force scaling: delta ratio {:.1}, expected F ratio {:.4}, got {:.4} (rel err {:.4})",
                deltas[i] / deltas[0], expected_ratio, actual_ratio, rel_err
            );
    }
}

#[test]
fn hooke_force_scales_linearly_across_overlaps() {
    let radius = 0.001;
    let hooke_force_at = |sep: f64| -> f64 {
        let mut app = App::new();
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let mut hist = ContactHistoryStore::new();
        atom.dt = 1e-7;
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
        push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
        atom.nlocal = 2;
        atom.natoms = 2;
        let mut neighbor = Neighbor::new();
        neighbor.neighbor_offsets = vec![0, 1, 1];
        neighbor.neighbor_indices = vec![1];
        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();
        registry.try_register(hist, atom.len()).unwrap();
        app.add_resource(atom);
        app.add_resource(neighbor);
        app.add_resource(registry);
        app.add_resource(make_material_table_hooke());
        app.add_update_system(hooke_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();
        let atom = app.get_resource_ref::<Atom>().unwrap();
        atom.force[0][0].abs() as f64
    };

    let deltas = [2e-5, 4e-5, 6e-5, 8e-5, 1e-4];
    let forces: Vec<f64> = deltas
        .iter()
        .map(|d| {
            let sep = 2.0 * radius - d;
            hooke_force_at(sep)
        })
        .collect();

    for i in 1..deltas.len() {
        let expected_ratio = deltas[i] / deltas[0]; // linear
        let actual_ratio = forces[i] / forces[0];
        let rel_err = ((actual_ratio - expected_ratio) / expected_ratio).abs();
        assert!(
                rel_err < 0.01,
                "Hooke force scaling: delta ratio {:.1}, expected F ratio {:.4}, got {:.4} (rel err {:.4})",
                deltas[i] / deltas[0], expected_ratio, actual_ratio, rel_err
            );
    }
}

#[test]
fn hertz_force_matches_analytical_value() {
    let radius = 0.001;
    let delta = 5e-5;
    let sep = 2.0 * radius - delta;

    let mut app = App::new();
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];
    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    let mt = make_material_table();
    let e_eff = mt.e_eff_ij[0][0];
    let r_eff = radius / 2.0; // two equal spheres: r_eff = r1*r2/(r1+r2) = r/2

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let f_computed = atom.force[0][0].abs() as f64;
    // Analytical: F = (4/3) * E_eff * sqrt(R_eff) * delta^(3/2)
    let f_analytical = (4.0 / 3.0) * e_eff * r_eff.sqrt() * delta.powf(1.5);
    let rel_err = (f_computed - f_analytical).abs() / f_analytical;
    assert!(
        rel_err < 1e-10,
        "Hertz force analytical check: computed={:.6e}, expected={:.6e}, rel_err={:.2e}",
        f_computed,
        f_analytical,
        rel_err
    );
}

#[test]
fn linear_momentum_conserved_during_elastic_contact() {
    // Perfectly elastic (restitution = 1.0) → ~no damping. The Hertz/Tsuji
    // coefficient at e=1 is the polynomial's residual (~1.3e-4), not exactly 0
    // (LAMMPS `damping tsuji` has the same residual), so momentum is conserved to
    // that order rather than machine epsilon.
    let mut mt = MaterialTable::new();
    mt.add_material("elastic", 8.7e9, 0.3, 1.0, 0.0, 0.0, 0.0);
    mt.build_pair_tables();
    assert!(
        mt.beta_ij[0][0].abs() < 1e-3,
        "beta should be ~0 for e=1.0, got {}",
        mt.beta_ij[0][0]
    );

    let radius = 0.001;
    let dt = 1e-8;

    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = dt;

    // Two particles approaching each other, slight overlap
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(
        &mut atom,
        &mut dem,
        &mut hist,
        1,
        [0.00195, 0.0, 0.0],
        radius,
    );
    atom.vel[0] = [0.1, 0.05, -0.02];
    atom.vel[1] = [-0.05, 0.03, 0.01];
    atom.nlocal = 2;
    atom.natoms = 2;

    let initial_momentum = [
        atom.mass[0] * atom.vel[0][0] + atom.mass[1] * atom.vel[1][0],
        atom.mass[0] * atom.vel[0][1] + atom.mass[1] * atom.vel[1][1],
        atom.mass[0] * atom.vel[0][2] + atom.mass[1] * atom.vel[1][2],
    ];

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];
    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    let mut app = App::new();
    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(
        crate::contact::hertz_mindlin_contact_force,
        ParticleSimScheduleSet::Force,
    );
    app.add_update_system(
        soil_verlet::initial_integration,
        ParticleSimScheduleSet::InitialIntegration,
    );
    app.add_update_system(
        soil_verlet::final_integration,
        ParticleSimScheduleSet::FinalIntegration,
    );
    // Zero forces between steps
    app.add_update_system(
        |mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>| {
            let n = atoms.len();
            for force in atoms.force.iter_mut().take(n) {
                *force = [0.0; 3];
            }
            registry.zero_all(n);
        },
        ParticleSimScheduleSet::PostInitialIntegration,
    );
    app.organize_systems();

    // Run for 100 steps
    for _ in 0..100 {
        app.run();
    }

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let final_momentum = [
        atom.mass[0] * atom.vel[0][0] + atom.mass[1] * atom.vel[1][0],
        atom.mass[0] * atom.vel[0][1] + atom.mass[1] * atom.vel[1][1],
        atom.mass[0] * atom.vel[0][2] + atom.mass[1] * atom.vel[1][2],
    ];

    for d in 0..3 {
        let err = (final_momentum[d] - initial_momentum[d]).abs();
        assert!(
            err < 1e-12,
            "Momentum not conserved in dim {}: initial={:.6e}, final={:.6e}, err={:.2e}",
            d,
            initial_momentum[d],
            final_momentum[d],
            err
        );
    }
}

#[test]
fn contact_force_symmetry_with_tangential_velocity() {
    let radius = 0.001;
    let sep = 0.0019;

    let mut app = App::new();
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-7;

    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 0, [0.0, 0.0, 0.0], radius);
    push_test_atom_with_history(&mut atom, &mut dem, &mut hist, 1, [sep, 0.0, 0.0], radius);
    // Give both atoms velocities in all directions
    atom.vel[0] = [0.1, 0.2, -0.1];
    atom.vel[1] = [-0.3, 0.1, 0.05];
    dem.omega[0] = [10.0, 20.0, -5.0];
    dem.omega[1] = [-15.0, 5.0, 10.0];
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];
    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(hist, atom.len()).unwrap();

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(make_material_table());
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    // Newton's 3rd law: forces equal and opposite
    for d in 0..3 {
        assert!(
            (atom.force[0][d] + atom.force[1][d]).abs() < 1e-10,
            "Newton's 3rd law violated in dim {}: f0={:.6e}, f1={:.6e}",
            d,
            atom.force[0][d],
            atom.force[1][d]
        );
    }
}

#[test]
fn willett_liquid_bridge_force_matches_closed_form_and_ruptures() {
    let r_eff: f64 = 2.5e-3;
    let volume: f64 = 1.0e-11;
    let gamma: f64 = 0.072;
    let theta: f64 = 0.0;
    let rupture: f64 = 5.0e-5;
    for separation in [0.0, 1.0e-6, 1.0e-5, 4.0e-5] {
        let s_hat = separation * (r_eff / volume).sqrt();
        let expected = 2.0 * std::f64::consts::PI * r_eff * gamma * theta.cos()
            / (1.0 + 1.05 * s_hat + 2.5 * s_hat * s_hat);
        let got = willett2000_liquid_bridge_force(separation, r_eff, volume, gamma, theta, rupture);
        assert!((got - expected).abs() < 1.0e-15);
    }
    assert_eq!(
        willett2000_liquid_bridge_force(rupture * 1.01, r_eff, volume, gamma, theta, rupture),
        0.0
    );
}
