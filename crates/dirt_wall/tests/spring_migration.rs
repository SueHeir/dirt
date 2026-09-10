//! A wall friction spring must survive an atom crossing a rank boundary.
//!
//! This is the crate-level pin for `dirt-wall-springs-atomdata`. It exercises
//! the **production** migration contract — `soil_core::ParticleStore`'s
//! `pack_migrant`/`append_migrant`, which is the whole of what an MPI exchange
//! sends and receives — rather than a harness that copies the store by hand.
//!
//! # What was wrong
//!
//! `Walls::tangential_springs` was a `HashMap<(u8, usize, u32), [f64; 3]>` on
//! the `Walls` **resource**. A resource is rank-local, so SOIL's exchange never
//! saw it: an atom in sustained frictional wall contact that migrated arrived
//! on the receiving rank with no entry for its tag, took `unwrap_or([0.0; 3])`,
//! and lost its whole tangential spring.
//!
//! `examples/mpi_wall_springs/RESULTS.md` measured that on shipped native DIRT
//! with no kernel anywhere: 15 crossings in a 2-rank run, every one in wall
//! contact, each discharging **100 %** of the tangential force (worst
//! `|dF_t|` = 1.196797e-4 N) for 1352 steps. The normal force agreed to exactly
//! zero on all 90 000 compared rows, so the spring history was the only thing
//! being dropped.
//!
//! # Why this test lives in `tests/` and not `src/tests.rs`
//!
//! `crates/dirt_wall/src/tests.rs` does not compile against the pinned
//! `soil_core` revision — 77 errors, all `field nlocal … is private`, `field
//! natoms … is private` and `method try_register is private`, pre-existing and
//! unrelated (`RESULTS.md` §7). An integration test is a separate compilation
//! unit, so `cargo test -p dirt_wall --test spring_migration` runs while that
//! is still broken.

use dirt_test_utils::{ParticleFixture, ParticleSpec};
use grass_app::prelude::*;
use soil_core::{Atom, AtomDataRegistry, ParticleSimScheduleSet, ParticleStore};

use dirt_wall::{wall_contact_force, WallMotion, WallPlane, WallSpringStore, Walls};

const RADIUS: f64 = 1.0e-3;
/// Overlap with the floor, so the contact is live from the first step.
const DELTA: f64 = 1.0e-5;
/// Tangential slip, so the Mindlin spring integrates away from zero.
const SLIP: f64 = 0.05;

/// A single frictional floor plane at `z = 0` with the normal pointing up.
fn floor() -> Walls {
    Walls {
        planes: vec![WallPlane {
            point_x: 0.0,
            point_y: 0.0,
            point_z: 0.0,
            normal_x: 0.0,
            normal_y: 0.0,
            normal_z: 1.0,
            material_index: 0,
            name: Some("floor".into()),
            bound_x_low: f64::NEG_INFINITY,
            bound_x_high: f64::INFINITY,
            bound_y_low: f64::NEG_INFINITY,
            bound_y_high: f64::INFINITY,
            bound_z_low: f64::NEG_INFINITY,
            bound_z_high: f64::INFINITY,
            velocity: [0.0; 3],
            motion: WallMotion::Static,
            origin: [0.0; 3],
            force_accumulator: 0.0,
            temperature: None,
        }],
        active: vec![true],
        cylinders: Vec::new(),
        cylinder_active: Vec::new(),
        spheres: Vec::new(),
        sphere_active: Vec::new(),
        regions: Vec::new(),
        region_active: Vec::new(),
        time: 0.0,
    }
}

/// An app holding `specs` atoms over the floor, with the wall force scheduled
/// and `WallSpringStore` registered exactly as `WallPlugin` registers it.
fn rank(specs: &[(u32, [f64; 3], [f64; 3])]) -> App {
    let mut builder = ParticleFixture::single(ParticleSpec::new(specs[0].0, specs[0].1, RADIUS));
    for spec in &specs[1..] {
        builder.push_particle(ParticleSpec::new(spec.0, spec.1, RADIUS));
    }
    let mut fixture = builder.build();
    for (i, spec) in specs.iter().enumerate() {
        fixture.atom.vel[i] = spec.2.map(|v| v as _);
    }
    fixture.register_atom_data(WallSpringStore::new());

    let mut app = fixture.into_app();
    app.add_resource(floor());
    app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app
}

/// The tangential spring app-row `row` holds against the floor.
fn spring_of(app: &App, row: usize) -> [f64; 3] {
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("AtomDataRegistry");
    let store = registry.expect::<WallSpringStore>("wall springs");
    store.tangential(row, 0, 0)
}

fn force_of(app: &App, row: usize) -> [f64; 3] {
    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    [
        atoms.force[row][0] as f64,
        atoms.force[row][1] as f64,
        atoms.force[row][2] as f64,
    ]
}

/// The whole claim, end to end: load a spring, send the atom through the real
/// migration wire, and show the receiving rank continues the same contact.
#[test]
fn a_loaded_wall_spring_crosses_a_rank_boundary_and_the_contact_continues() {
    let contact = (7u32, [0.0, 0.0, RADIUS - DELTA], [SLIP, 0.0, 0.0]);
    // The receiving rank already owns an atom, and it is nowhere near the
    // floor — so it contributes nothing, and the arriving row lands at index 1.
    let bystander = (99u32, [0.02, 0.02, 0.5], [0.0; 3]);

    // ── the sending rank loads a spring ─────────────────────────────────────
    let mut sender = rank(&[contact]);
    sender.run();
    let sent = spring_of(&sender, 0);
    assert!(
        sent.iter().any(|c| *c != 0.0),
        "the fixture must load a non-zero tangential spring for this to prove \
         anything, got {sent:?}"
    );

    // ── it goes into the exchange message ───────────────────────────────────
    // `pack_migrant` is `Atom::pack_exchange` plus every registered `AtomData`,
    // and that is the entirety of what crosses a rank boundary with an atom.
    let message = {
        let registry = sender
            .get_resource_ref::<AtomDataRegistry>()
            .expect("registry");
        let cell = sender
            .resource_cell(std::any::TypeId::of::<Atom>())
            .expect("Atom");
        let mut boxed = cell.borrow_mut();
        let atoms = boxed.downcast_mut::<Atom>().expect("Atom");
        let mut buf = Vec::new();
        ParticleStore::new(atoms, &registry)
            .pack_migrant(0, &mut buf)
            .expect("a local atom with a valid layout packs");
        buf
    };

    // ── and comes back out on a rank that has never held the tag ────────────
    let mut receiver = rank(&[bystander]);
    let arrived_row = {
        let registry = receiver
            .get_resource_ref::<AtomDataRegistry>()
            .expect("registry");
        let cell = receiver
            .resource_cell(std::any::TypeId::of::<Atom>())
            .expect("Atom");
        let mut boxed = cell.borrow_mut();
        let atoms = boxed.downcast_mut::<Atom>().expect("Atom");
        let before = atoms.nlocal() as usize;
        let consumed = ParticleStore::new(atoms, &registry)
            .append_migrant(&message)
            .expect("the message decodes");
        assert_eq!(
            consumed,
            message.len(),
            "the receiving rank must consume exactly the record the sender wrote"
        );
        assert_eq!(atoms.nlocal() as usize, before + 1);
        before
    };
    assert_eq!(arrived_row, 1, "the migrant lands after the rank's own atom");
    assert_eq!(
        receiver
            .get_resource_ref::<Atom>()
            .expect("Atom")
            .tag[arrived_row],
        contact.0
    );

    let received = spring_of(&receiver, arrived_row);
    assert_eq!(
        received, sent,
        "the tangential spring must arrive bit for bit; before \
         dirt-wall-springs-atomdata it arrived as [0, 0, 0], which cost 100 % of \
         the tangential force for ~1352 steps per crossing"
    );

    // ── and the contact simply continues ────────────────────────────────────
    // Both ranks now hold identical state for this atom, so one more step must
    // advance the spring and apply the force identically. This is the part a
    // "the bits are in the buffer" assertion cannot reach: it is what says the
    // *contact* survived, not merely the number.
    let sender_force_before = force_of(&sender, 0);
    let receiver_force_before = force_of(&receiver, arrived_row);
    sender.run();
    receiver.run();

    assert_eq!(
        spring_of(&receiver, arrived_row),
        spring_of(&sender, 0),
        "after the crossing the spring must keep advancing identically on the \
         receiving rank"
    );
    let sender_step: Vec<f64> = force_of(&sender, 0)
        .iter()
        .zip(sender_force_before)
        .map(|(a, b)| a - b)
        .collect();
    let receiver_step: Vec<f64> = force_of(&receiver, arrived_row)
        .iter()
        .zip(receiver_force_before)
        .map(|(a, b)| a - b)
        .collect();
    assert_eq!(
        receiver_step, sender_step,
        "the wall force on the step after the crossing must be identical on both \
         ranks, to the bit"
    );
    assert!(
        sender_step[0] != 0.0,
        "the comparison is only meaningful if a tangential force was applied"
    );

    println!(
        "wall spring migration: tag {} carried [{:.6e}, {:.6e}, {:.6e}] m across a \
         {} f64 exchange message intact; the next step's wall force matches the \
         sending rank's bit for bit ({:.6e} N tangential)",
        contact.0,
        sent[0],
        sent[1],
        sent[2],
        message.len(),
        sender_step[0],
    );
}

/// A contact that **ends** must still lose its spring. Pruning used to be free:
/// both maps were rebuilt from scratch every step, so an untouched key simply
/// was not re-inserted. The per-atom store has to do it explicitly, and getting
/// it wrong would leave a stale spring to be re-applied when the particle came
/// back — a different physics bug in the opposite direction.
#[test]
fn a_contact_that_ends_still_loses_its_spring() {
    let mut app = rank(&[(7u32, [0.0, 0.0, RADIUS - DELTA], [SLIP, 0.0, 0.0])]);
    app.run();
    assert!(
        spring_of(&app, 0).iter().any(|c| *c != 0.0),
        "the contact must be loaded before it can be ended"
    );

    // Lift the atom clear of the floor and step again.
    {
        let cell = app
            .resource_cell(std::any::TypeId::of::<Atom>())
            .expect("Atom");
        let mut boxed = cell.borrow_mut();
        let atoms = boxed.downcast_mut::<Atom>().expect("Atom");
        atoms.pos[0][2] = (10.0 * RADIUS) as _;
    }
    app.run();

    assert_eq!(
        spring_of(&app, 0),
        [0.0; 3],
        "a contact that ended must not keep its spring"
    );
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("registry");
    let store = registry.expect::<WallSpringStore>("wall springs");
    assert!(
        store.row(0).is_empty(),
        "the ended contact's entry must be pruned, not merely zeroed"
    );
}

/// SOIL moves the rows; this checks the two structural operations it uses do
/// what the atom store does, because a row that drifts off its atom is a
/// silently wrong friction history rather than a crash.
#[test]
fn rows_follow_their_atoms_through_a_sort_and_a_deletion() {
    use soil_core::AtomData;

    let mut store = WallSpringStore::new();
    for row in 0..3 {
        store.set_tangential(row, 0, 0, [row as f64 + 1.0, 0.0, 0.0]);
    }

    // A spatial sort: 0 -> 2, 1 -> 0, 2 -> 1.
    unsafe { store.apply_permutation(&[2, 0, 1], 3) };
    assert_eq!(store.tangential(0, 0, 0), [3.0, 0.0, 0.0]);
    assert_eq!(store.tangential(1, 0, 0), [1.0, 0.0, 0.0]);
    assert_eq!(store.tangential(2, 0, 0), [2.0, 0.0, 0.0]);

    // A deletion: the last row is swapped into the hole.
    unsafe { store.swap_remove(0) };
    assert_eq!(store.len(), 2);
    assert_eq!(store.tangential(0, 0, 0), [2.0, 0.0, 0.0]);
    assert_eq!(store.tangential(1, 0, 0), [1.0, 0.0, 0.0]);
}

/// Several walls at once, through the wire. The row is variable-length and
/// keyed by `(kind, index)`, so an atom touching a plane and a cylinder must
/// arrive with both — and with them the right way round.
#[test]
fn a_row_carries_every_wall_it_touches() {
    use soil_core::AtomData;

    let mut sender = WallSpringStore::new();
    sender.set_tangential(0, 0, 3, [1.0, 2.0, 3.0]);
    sender.set_rolling(0, 0, 3, [4.0, 5.0, 6.0]);
    sender.set_tangential(0, 1, 0, [7.0, 8.0, 9.0]);
    sender.set_tangential(0, 2, 11, [-1.0, -2.0, -3.0]);

    let mut wire = Vec::new();
    sender.pack(0, &mut wire);

    let mut receiver = WallSpringStore::new();
    let used = unsafe { receiver.unpack(&wire) };
    assert_eq!(used, wire.len(), "the decoder consumes exactly the record");
    assert_eq!(receiver.row(0).len(), 3);

    assert_eq!(receiver.tangential(0, 0, 3), [1.0, 2.0, 3.0]);
    assert_eq!(receiver.rolling(0, 0, 3), [4.0, 5.0, 6.0]);
    assert_eq!(receiver.tangential(0, 1, 0), [7.0, 8.0, 9.0]);
    assert_eq!(receiver.tangential(0, 2, 11), [-1.0, -2.0, -3.0]);
    // A wall the atom was not touching stays at the "no history" default.
    assert_eq!(receiver.tangential(0, 3, 0), [0.0; 3]);

    // Arrivals are inactive: the receiving rank's own force pass has to touch
    // them or they are pruned, which is what keeps "a contact that ended loses
    // its spring" true across a crossing too.
    assert!(receiver.row(0).iter().all(|e| !e.active));
    receiver.prune(1);
    assert!(receiver.row(0).is_empty());
}

/// A frictionless wall must not allocate a spring row at all — the crate
/// documents frictionless contacts as byte-for-byte identical to a pure normal
/// contact, and an empty row is what keeps the migration message that size too.
#[test]
fn a_frictionless_wall_stores_nothing() {
    let mut app = rank(&[(7u32, [0.0, 0.0, RADIUS - DELTA], [SLIP, 0.0, 0.0])]);
    {
        let cell = app
            .resource_cell(std::any::TypeId::of::<dirt_atom::MaterialTable>())
            .expect("MaterialTable");
        let mut boxed = cell.borrow_mut();
        let mt = boxed
            .downcast_mut::<dirt_atom::MaterialTable>()
            .expect("MaterialTable");
        mt.friction_ij[0][0] = 0.0;
        mt.rolling_friction_ij[0][0] = 0.0;
    }
    app.run();

    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("registry");
    let store = registry.expect::<WallSpringStore>("wall springs");
    assert!(
        store.row(0).is_empty(),
        "a frictionless wall contact must store no spring history"
    );
    // And it is still a live normal contact, so this is not vacuous.
    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    assert!(atoms.force[0][2] > 0.0 as _, "the normal force must be live");
}
