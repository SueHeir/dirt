//! Emit public `SimulationFixture` measurements for the companion contract sweep.

use dirt_atom::DemAtom;
use dirt_test_utils::{ParticleFixture, ParticleSpec, DEFAULT_DEM_TIMESTEP};

fn main() {
    let mut custom = ParticleFixture::single(ParticleSpec::new(7, [0.0; 3], 0.002));
    let second = custom.push_particle(ParticleSpec::new(11, [0.003, 0.0, 0.0], 0.001));
    let third = custom.push_particle(ParticleSpec::new(19, [0.006, 0.0, 0.0], 0.0015));
    custom.add_pair((0, second));
    custom.add_pair((second, third));

    println!("case,atom_rows,nlocal,natoms,dem_rows,csr_offsets,csr_indices,materials,pair_rows,pair_columns,dt");
    emit("default_pair", ParticleFixture::new().build());
    emit("custom_chain", custom.build());
}

fn emit(name: &str, fixture: dirt_test_utils::SimulationFixture) {
    let dem_rows = fixture
        .registry
        .expect::<DemAtom>("fixture validation")
        .radius
        .len();
    let offsets = fixture
        .neighbor
        .neighbor_offsets
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(";");
    let indices = fixture
        .neighbor
        .neighbor_indices
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(";");
    println!(
        "{name},{},{},{},{},{offsets},{indices},{},{},{},{:.17e}",
        fixture.atom.len(),
        fixture.atom.nlocal,
        fixture.atom.natoms,
        dem_rows,
        fixture.materials.names.len(),
        fixture.materials.e_eff_ij.len(),
        fixture.materials.e_eff_ij[0].len(),
        fixture.atom.dt,
    );
    assert_eq!(fixture.atom.dt, DEFAULT_DEM_TIMESTEP);
}
