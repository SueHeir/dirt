# Writing a DIRT Test

The `dirt_test_utils` crate provides helpers to set up `Atom`, `GroupRegistry`,
`CommResource`, `DemAtom`, and `MaterialTable` instances quickly, so a unit test
can focus on the logic under test rather than boilerplate. Each helper produces
a minimal, valid object.

## Quick start

```rust
use dirt_test_utils::{make_atoms, make_group_registry, make_single_comm};

let atom = make_atoms(3); // 3 atoms at (0,0,0), (1,0,0), (2,0,0)
let groups = make_group_registry("all", vec![true, true, true]);
let _comm = make_single_comm();
```

`make_atoms(n)` places atom `i` at `(i, 0, 0)` with radius 0.5, mass 1.0, and
`dt = 0.001`. `make_single_comm` returns a single-process comm requiring no MPI
init. `make_material_table` returns a single ready-to-use "glass" material with
pair tables already built (8.7 GPa, ν = 0.3, restitution 0.95, friction 0.4).

For DEM atoms, `push_dem_test_atom` appends one solid sphere (density
2500 kg/m³; mass from `m = ρ·(4/3)πr³`; `I = 0.4·m·r²`) to both the `Atom` and
`DemAtom` arrays in parallel.

## The full pattern: a force-law test

These helpers cover atom/material **construction** only. A force-law test also
needs an `App` (from `grass_app`) with a neighbor list and the system under test
scheduled into it. That assembly lives in the *consumer* crate — `dirt_test_utils`
deliberately depends only on `soil_core` + `dirt_atom`, so the snippet below is
illustrative. It mirrors the real pattern in `dirt_granular`'s `contact.rs`
tests:

```rust
use dirt_test_utils::{push_dem_test_atom, make_material_table};
use dirt_atom::DemAtom;
use soil_core::{Atom, AtomDataRegistry, Neighbor, ParticleSimScheduleSet};
use grass_app::prelude::*;

let mut app = App::new();
let radius = 0.001;

// 1. Build Atom + DemAtom in parallel with the helper.
let mut atom = Atom::new();
let mut dem = DemAtom::new();
atom.dt = 1e-7;                     // DEM needs a tiny dt; the helper does NOT set it
push_dem_test_atom(&mut atom, &mut dem, 0, [0.0, 0.0, 0.0], radius);
push_dem_test_atom(&mut atom, &mut dem, 1, [0.0019, 0.0, 0.0], radius);

// 2. Manually set the counts — push_dem_test_atom does NOT touch them.
atom.nlocal = 2;
atom.natoms = 2;

// 3. Build a Neighbor list by hand (CSR offsets/indices for the pair 0–1).
let mut neighbor = Neighbor::new();
neighbor.neighbor_offsets = vec![0, 1, 1];
neighbor.neighbor_indices = vec![1];

// 4. Register per-atom data (and any history store) into the registry.
let mut registry = AtomDataRegistry::new();
registry.register(dem);

// 5. Add resources, schedule the system under test, organize, run.
app.add_resource(atom);
app.add_resource(neighbor);
app.add_resource(registry);
app.add_resource(make_material_table());
app.add_update_system(my_force_system, ParticleSimScheduleSet::Force);
app.organize_systems();
app.run();

// 6. Read results back off the App resources and assert.
let atom = app.get_resource_ref::<Atom>().unwrap();
assert!(atom.force[0][0] < 0.0);
```

## Two footguns

- **`nlocal` / `natoms` are NOT set by `push_dem_test_atom`.** It only *appends*
  one atom to the arrays; it leaves `atom.nlocal` and `atom.natoms` untouched.
  Force systems iterate `0..nlocal`, so if you forget to set them (they default
  to 0) **your system silently does nothing** — no panic, just an empty loop and
  a passing-but-meaningless test. Set both after your last `push_*` call.
- **DEM needs a tiny timestep.** `make_atoms` sets `dt = 0.001`, far too large
  for stiff DEM contacts. For any contact/bond test override `atom.dt` to
  roughly `1e-7` (the value the real `dirt_granular`/`dirt_bond` tests use) so a
  single `app.run()` step does not blow the contact up.
- **`AtomDataRegistry` registration order matters.** Every per-atom column your
  system reads (`DemAtom`, a history store, …) must be registered into the
  `AtomDataRegistry` *before* the registry is added to the `App`. Systems fetch
  these columns by type at run time, so a missing registration is a panic when
  the system runs, not a compile error.

## A history-store wrapper: `push_test_atom_with_history`

`push_dem_test_atom` populates `Atom` + `DemAtom` only. A contact or bond test
also needs a per-atom slot in the relevant **history store**
(`ContactHistoryStore` for `dirt_granular`, the bond history for `dirt_bond`).
Because those stores are crate-internal types, `dirt_test_utils` cannot push into
them — so each crate defines a thin local wrapper in its `#[cfg(test)]` module:

```rust
// local to the consumer crate's test module
fn push_test_atom_with_history(
    atom: &mut Atom, dem: &mut DemAtom, history: &mut ContactHistoryStore,
    tag: u32, pos: [f64; 3], radius: f64,
) {
    push_dem_test_atom(atom, dem, tag, pos, radius);
    history.push_atom(); // add the matching per-atom history slot
}
```

This local-wrapper pattern is deliberate: keeping `ContactHistoryStore` out of
`dirt_test_utils` preserves the crate's tight dependency footprint
(`soil_core` + `dirt_atom` only). Reproduce it in any crate that carries
per-contact or per-bond history.

## When you need `make_group_registry`

`make_group_registry(name, mask)` builds a one-group `GroupRegistry` whose
membership is the boolean `mask`. It is what **fix** tests use to target a subset
of atoms (`dirt_fixes` tests rely on it). Pure contact-law tests do *not* need
it — the contact systems iterate the neighbor list, not a group — so do not add
it reflexively to a force-law test.

## Customizing materials beyond the glass baseline

`make_material_table()` returns the neutral `"glass"` material (8.7 GPa, ν = 0.3,
e = 0.95, μ = 0.4) with pair tables already built — the no-special-behaviour
case. To test cohesion, JKR adhesion, rolling resistance, SDS, Hooke, or
twisting, build a local table in your test module, register the materials, and
call `build_pair_tables()` yourself. `dirt_granular/contact.rs` has the canonical
set of these local helpers (`make_material_table_jkr`, `make_material_table_sds_rolling`,
…) — follow that pattern rather than extending the shared helper.

`push_dem_test_atom` hardcodes density at 2500 kg/m³ with no override. If a test
needs a different density, bypass the helper: call `atom.push_test_atom` and
populate the `DemAtom` fields (mass, `inv_inertia`, quaternion, …) directly.

## What this crate does NOT provide

- **No `App` / `Neighbor` builders.** You construct and wire those in the
  consuming crate; the neighbor list in particular is built by hand (raw CSR
  `neighbor_offsets` / `neighbor_indices`).
- **No custom assertions or test macros.** Use plain `assert!` / `assert_eq!`.
- **No scheduling / plugin helpers.** Schedule sets and systems are added
  directly via the `grass_app` / `grass_scheduler` API in the consumer.
