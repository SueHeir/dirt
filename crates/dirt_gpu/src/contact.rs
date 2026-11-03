//! [`GpuContactPlugin`]: the contact force on the GPU, everything else on the CPU.
//!
//! This is the LAMMPS GPU-package model. Each step uploads positions, velocities and
//! angular velocities, evaluates the contact force, and adds force and torque onto the
//! CPU arrays. Integration, walls, fixes, MPI and output are untouched and do not know.
//!
//! # Contact history
//!
//! Springs live on the GPU between neighbor rebuilds, keyed by partner tag so they
//! survive ghost reordering. At each rebuild they come down into DIRT's own
//! [`ContactHistoryStore`], and go back up after exchange, borders and sorting. Those
//! steps already carry that store correctly, so sorting and migration need nothing
//! further here. The store's canonical (lower-tag) frame is converted on the way through,
//! because each GPU atom keeps its own half of a spring.
//!
//! A restart file written between rebuilds holds the history from the last rebuild.

use std::collections::HashMap;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_granular::tangential::{ContactHistoryStore, CONTACT_HISTORY_LEN};
use dirt_schedule::CONTACT_FORCE;
use grass_app::prelude::*;
use grass_gpu::{GpuDevice, GpuPlugin};
use grass_scheduler::prelude::*;
use soil_core::{
    register_atom_data, Accum, Atom, BondStore, CommState, Optional, ParticleSimScheduleSet,
    ParticlesWith, Read, Real, Write,
};
use soil_gpu::neighbor_slots::{MAX_CONTACTS, SENTINEL};
use soil_gpu::{DeviceColumn, GpuReal, GpuState, Grid};

use crate::granular_force::{
    GranularColumns, GranularConfig, GranularForce, MaterialPairs, PAIR_FIELDS, ROLLING_CONSTANT,
    ROLLING_SDS, TWISTING_CONSTANT, TWISTING_MARSHALL, TWISTING_SDS,
};

/// Cells allowed per atom. Memory stays bounded if one particle flies far.
const CELLS_PER_ATOM: usize = 8;

/// The contact table for every pair of materials.
pub(crate) fn material_pairs(mt: &MaterialTable) -> MaterialPairs {
    let ntypes = mt.e_eff_ij.len();
    let mut rows = Vec::with_capacity(ntypes * ntypes * PAIR_FIELDS);
    for i in 0..ntypes {
        for j in 0..ntypes {
            rows.extend(pair_row(mt, i, j));
        }
    }
    let rolling_sds = mt.rolling_model == "sds";
    MaterialPairs {
        ntypes,
        rows,
        rolling_model: if rolling_sds {
            ROLLING_SDS
        } else {
            ROLLING_CONSTANT
        },
        twisting_model: match mt.twisting_model.as_str() {
            "sds" => TWISTING_SDS,
            "marshall" => TWISTING_MARSHALL,
            _ => TWISTING_CONSTANT,
        },
    }
}

/// One pair's parameters, in `PAIR_FIELDS` order.
pub(crate) fn pair_row(mt: &MaterialTable, i: usize, j: usize) -> [GpuReal; PAIR_FIELDS] {
    [
        mt.e_eff_ij[i][j] as GpuReal,
        mt.beta_ij[i][j] as GpuReal,
        mt.g_eff_ij[i][j] as GpuReal,
        mt.friction_ij[i][j] as GpuReal,
        mt.rolling_friction_ij[i][j] as GpuReal,
        mt.twisting_friction_ij[i][j] as GpuReal,
        mt.rolling_stiffness_ij[i][j] as GpuReal,
        mt.rolling_damping_ij[i][j] as GpuReal,
        mt.twisting_stiffness_ij[i][j] as GpuReal,
        mt.twisting_damping_ij[i][j] as GpuReal,
    ]
}

/// What the kernels leave out, for one material table. Empty means they can run it.
fn unsupported_features(mt: &MaterialTable) -> Vec<String> {
    let mut found = Vec::new();
    let any = |table: &[Vec<f64>]| table.iter().flatten().any(|&value| value != 0.0);
    if any(&mt.cohesion_energy_ij) {
        found.push("cohesion energy".to_string());
    }
    if any(&mt.surface_energy_ij) {
        found.push("surface energy".to_string());
    }
    if !mt.limit_damping {
        found.push("limit_damping = false".to_string());
    }
    if mt.contact_model != "hertz" {
        found.push(format!("contact_model = {:?}", mt.contact_model));
    }
    if mt.tangential_model != "history" {
        found.push(format!("tangential_model = {:?}", mt.tangential_model));
    }
    if mt.liquid_bridge_model != "off" {
        found.push(format!(
            "liquid_bridge_model = {:?}",
            mt.liquid_bridge_model
        ));
    }
    found
}

pub(crate) fn require_supported(mt: &MaterialTable, dem: &DemAtom, nall: usize) {
    let mut found = unsupported_features(mt);
    if dem.body_id[..nall].iter().any(|&body| body != 0.0) {
        found.push("rigid clumps".to_string());
    }
    assert!(
        found.is_empty(),
        "The GPU kernels run Hertz-Mindlin with constant-torque rolling and twisting. This run \
         uses: {}. Use HertzMindlinContactPlugin for it.",
        found.join(", ")
    );
}

struct Resident {
    state: GpuState,
    columns: GranularColumns,
    omega: usize,
}

impl Resident {
    fn new(gpu: &GpuDevice, n: usize, grid: Grid, mt: &MaterialTable, dt: GpuReal) -> Self {
        let mut state = GpuState::new(gpu.clone(), n, grid);
        // Gravity off: this contributes only the contact force, added on the CPU
        // alongside gravity, walls, and other fixes.
        state.set_params(dt, [0.0; 3]);
        let omega = state.add_aux_dof();
        let hook = GranularForce::new(
            &state,
            &grid,
            omega,
            &vec![0.0; n],
            &material_pairs(mt),
            GranularConfig::new(dt),
        );
        let columns = hook.columns();
        state.add_force_hook(Box::new(hook));
        Self {
            state,
            columns,
            omega,
        }
    }
}

/// The GPU contact state, as an App resource.
#[derive(Default)]
pub struct GpuContact {
    resident: Option<Resident>,
    /// Tags of the rows the GPU springs belong to, as of the last upload.
    row_tags: Vec<u32>,
    /// The springs are in the [`ContactHistoryStore`], and the GPU copy is stale.
    history_on_host: bool,
    /// Bonds counted when the bond lists were last sent.
    bond_total: usize,
}

/// Every row's bond partners by tag, for the contact kernel: CSR offsets and tags.
pub(crate) fn bond_lists(bonds: Option<&BondStore>, rows: usize) -> (Vec<u32>, Vec<u32>) {
    let mut start = Vec::with_capacity(rows + 1);
    let mut tags = Vec::new();
    start.push(0);
    for i in 0..rows {
        let list = bonds
            .and_then(|b| b.bonds.get(i))
            .map_or(&[][..], Vec::as_slice);
        tags.extend(list.iter().map(|bond| bond.partner_tag));
        start.push(tags.len() as u32);
    }
    (start, tags)
}

/// Sends the CPU's bonds to the contact kernel. Their force stays on the CPU here, so no end is
/// broken as far as the kernel knows: a broken bond is simply gone from the list.
fn write_cpu_bonds(
    gpu: &GpuDevice,
    columns: &GranularColumns,
    bonds: Option<&BondStore>,
    rows: usize,
) {
    let (start, tags) = bond_lists(bonds, rows);
    let alive = DeviceColumn::from_slice(gpu, &vec![0u32; tags.len().max(1)]);
    columns.write_bonds(&start, &tags, alive);
}

fn cast3(column: &[[Real; 3]], out: &mut Vec<GpuReal>) {
    out.clear();
    out.extend(
        column
            .iter()
            .flat_map(|v| [v[0] as GpuReal, v[1] as GpuReal, v[2] as GpuReal]),
    );
}

/// Sliding and rolling displacement flip between a contact's two ends, since the normal and the
/// relative motion both do. Twist, in slot 6, does not: it integrates relative spin along the
/// normal, and both of those flip.
fn frame_sign(slot: usize, sign: f64) -> f64 {
    if slot == 6 {
        1.0
    } else {
        sign
    }
}

/// Packs the store into GPU slots for the first `nlocal` rows. A pair may be recorded on
/// only one side, or on a ghost's copy of its owner's list, so every list is searched
/// for records that mention a local atom.
fn store_to_slots(
    history: &ContactHistoryStore,
    tag: &[u32],
    nlocal: usize,
    payload_len: usize,
) -> (Vec<u32>, Vec<GpuReal>) {
    let local: HashMap<u32, usize> = tag[..nlocal]
        .iter()
        .enumerate()
        .map(|(i, &t)| (t, i))
        .collect();
    let mut partner = vec![SENTINEL; nlocal * MAX_CONTACTS];
    let mut payload = vec![0.0 as GpuReal; nlocal * MAX_CONTACTS * payload_len];
    let mut used = vec![0usize; nlocal];
    let mut put = |row: usize, other: u32, canonical: &[f64; CONTACT_HISTORY_LEN]| {
        let base = row * MAX_CONTACTS;
        let slots = &mut partner[base..base + used[row]];
        if used[row] == MAX_CONTACTS || slots.contains(&other) {
            return;
        }
        let k = base + used[row];
        partner[k] = other;
        let sign = if tag[row] < other { 1.0 } else { -1.0 };
        for d in 0..payload_len {
            payload[k * payload_len + d] = (frame_sign(d, sign) * canonical[d]) as GpuReal;
        }
        used[row] += 1;
    };
    // A row's own records win over what its partners hold.
    for (row, list) in history.contacts.iter().enumerate().take(nlocal) {
        for (other, canonical, _) in list {
            put(row, *other, canonical);
        }
    }
    for (holder, list) in history.contacts.iter().enumerate().take(tag.len()) {
        for (other, canonical, _) in list {
            if let Some(&row) = local.get(other) {
                put(row, tag[holder], canonical);
            }
        }
    }
    (partner, payload)
}

/// Downloads the springs into the store, so that exchange, sorting and borders carry them.
fn gpu_contact_history_to_host(
    atoms: Res<Atom>,
    particles: ParticlesWith<'_, Write<ContactHistoryStore>>,
    mut contact: ResMut<GpuContact>,
) {
    if contact.history_on_host {
        return;
    }
    contact.history_on_host = true;
    let GpuContact {
        resident, row_tags, ..
    } = &mut *contact;
    let Some(resident) = resident.as_ref() else {
        return;
    };
    let rows = row_tags.len();
    let (partner, payload) = resident.columns.read_slots(rows);
    let payload_len = payload.len() / (rows * MAX_CONTACTS).max(1);

    particles.with(|mut history| {
        let nlocal = atoms.nlocal() as usize;
        if history.contacts.len() < nlocal {
            history.contacts.resize_with(nlocal, Vec::new);
        }
        // Atoms only reorder at this rebuild, later on, but a fix may have removed some.
        let moved = atoms.tag[..rows.min(atoms.len())] != row_tags[..];
        let index: HashMap<u32, usize> = if moved {
            atoms.tag[..nlocal]
                .iter()
                .enumerate()
                .map(|(i, &t)| (t, i))
                .collect()
        } else {
            HashMap::new()
        };
        let mut full_rows = 0;
        for (r, &owner) in row_tags.iter().enumerate() {
            let row = if moved {
                match index.get(&owner) {
                    Some(&row) => row,
                    None => continue,
                }
            } else {
                r
            };
            let list = &mut history.contacts[row];
            list.clear();
            let base = r * MAX_CONTACTS;
            for k in base..base + MAX_CONTACTS {
                let other = partner[k];
                if other == SENTINEL {
                    continue;
                }
                let sign = if owner < other { 1.0 } else { -1.0 };
                let mut canonical = [0.0; CONTACT_HISTORY_LEN];
                for d in 0..payload_len {
                    canonical[d] = frame_sign(d, sign) * payload[k * payload_len + d] as f64;
                }
                list.push((other, canonical, false));
            }
            full_rows += (list.len() == MAX_CONTACTS) as usize;
        }
        if full_rows > 0 {
            eprintln!(
                "WARNING: {full_rows} atoms reached the GPU limit of {MAX_CONTACTS} contacts; \
                 further contacts on them keep no spring history"
            );
        }
    });
}

/// Replaces the CPU contact force.
pub fn gpu_contact_force(
    gpu: Res<GpuDevice>,
    mut atoms: ResMut<Atom>,
    particles: ParticlesWith<
        '_,
        (
            Write<DemAtom>,
            Read<ContactHistoryStore>,
            Optional<Read<BondStore>>,
        ),
    >,
    material_table: Res<MaterialTable>,
    mut contact: ResMut<GpuContact>,
    mut scratch: Local<Vec<GpuReal>>,
) {
    let nall = atoms.len();
    let nlocal = atoms.nlocal() as usize;
    let dt = atoms.dt as GpuReal;
    particles.with(|(mut dem, history, bonds)| {
        let bond_total = bonds
            .as_deref()
            .map_or(0, |b| b.bonds.iter().map(Vec::len).sum::<usize>());
        if contact.history_on_host {
            require_supported(&material_table, &dem, nall);
            cast3(&atoms.pos[..nall], &mut scratch);
            let radius: Vec<GpuReal> = dem.radius[..nall].iter().map(|&r| r as GpuReal).collect();
            let r_max = radius
                .iter()
                .copied()
                .fold(GpuReal::MIN_POSITIVE, GpuReal::max);
            let grid =
                Grid::from_positions_capped(&scratch, 2.0 * r_max, nall.max(1) * CELLS_PER_ATOM);
            let resident = contact
                .resident
                .get_or_insert_with(|| Resident::new(&gpu, nall, grid, &material_table, dt));
            resident.state.set_active(nall, nlocal, grid);
            let inv_mass: Vec<GpuReal> = atoms.mass[..nall]
                .iter()
                .map(|&m| 1.0 / m as GpuReal)
                .collect();
            resident.state.write_inv_mass(&inv_mass);
            let zero = vec![0.0; nall];
            resident.state.set_aux_inv_coeff(resident.omega, &zero);
            resident.columns.write_radius(&radius);
            resident.columns.write_tag(&atoms.tag[..nall]);
            resident.columns.write_type(&atoms.atom_type[..nall]);
            let payload_len = material_pairs(&material_table).payload_len();
            let (partner, payload) =
                store_to_slots(&history, &atoms.tag[..nall], nlocal, payload_len);
            resident.columns.write_slots(nlocal, &partner, &payload);

            write_cpu_bonds(&gpu, &resident.columns, bonds.as_deref(), nall);
            contact.bond_total = bond_total;

            contact.row_tags.clear();
            contact.row_tags.extend_from_slice(&atoms.tag[..nlocal]);
            contact.history_on_host = false;
        }

        // A CPU bond broke since the last rebuild, so those two atoms can touch now.
        if bond_total != contact.bond_total {
            let resident = contact.resident.as_ref().expect("built on the first step");
            write_cpu_bonds(&gpu, &resident.columns, bonds.as_deref(), nall);
            contact.bond_total = bond_total;
        }
        let resident = contact.resident.as_mut().expect("built on the first step");
        cast3(&atoms.pos[..nall], &mut scratch);
        resident.state.write_pos(&scratch);
        cast3(&atoms.vel[..nall], &mut scratch);
        resident.state.write_vel(&scratch);
        scratch.clear();
        scratch.extend(
            dem.omega[..nall]
                .iter()
                .flat_map(|w| [w[0] as GpuReal, w[1] as GpuReal, w[2] as GpuReal]),
        );
        resident.state.set_aux_state(resident.omega, &scratch);
        resident.state.eval_force_once();
        let omega = resident.omega;
        resident.state.sync_force_and_rate(omega);

        // Adds, like the CPU force it replaces, so gravity and walls compose. Locals only:
        // each local atom's sum is complete, so nothing needs sending back from ghosts.
        for (i, f) in resident.state.force()[..3 * nlocal]
            .chunks_exact(3)
            .enumerate()
        {
            for (out, f) in atoms.force[i].iter_mut().zip(f) {
                *out += *f as Accum;
            }
        }
        let rate = resident.state.aux_rate(omega)[..3 * nlocal].chunks_exact(3);
        for (out, t) in dem.torque[..nlocal].iter_mut().zip(rate) {
            for (out, t) in out.iter_mut().zip(t) {
                *out += *t as f64;
            }
        }
    });
}

/// GPU drop-in for `HertzMindlinContactPlugin`. [`GranularGpuPlugins`] is
/// `GranularDefaultPlugins` with this swapped in:
///
/// ```rust,ignore
/// app.add_plugins(CorePlugins).add_plugins(GranularGpuPlugins);
/// ```
///
/// It covers Hertz-Mindlin with constant-torque rolling and twisting friction and any number
/// of materials. It panics at the first step on spring-based rolling or twisting, cohesion,
/// adhesion, liquid bridges, other contact or tangential models, and rigid clumps. It also panics on a machine with no GPU.
/// Bonded particles work: add `DemBondPlugin`, whose force stays on the CPU, and bonded
/// pairs and their neighbors are left out of contact as `BondStore::are_excluded` does.
/// It does not add to the virial, so per-contact stress output needs the CPU force.
pub struct GpuContactPlugin;

impl Plugin for GpuContactPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(GpuPlugin::default());
        register_atom_data!(app, ContactHistoryStore::new());
        app.add_resource(GpuContact {
            history_on_host: true,
            ..Default::default()
        });
        app.add_update_system(
            gpu_contact_history_to_host.run_if(in_state(CommState::FullRebuild)),
            ParticleSimScheduleSet::PreExchange,
        );
        app.add_update_system(
            gpu_contact_force.label(CONTACT_FORCE),
            ParticleSimScheduleSet::Force,
        );
    }
}

/// `GranularDefaultPlugins` with the contact force on the GPU.
pub struct GranularGpuPlugins;

impl PluginGroup for GranularGpuPlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(dirt_atom::DemAtomPlugin)
            .add(dirt_atom::DemAtomInsertPlugin)
            .add(soil_verlet::VelocityVerletPlugin::new())
            .add(GpuContactPlugin)
            .add(dirt_granular::RotationalDynamicsPlugin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dirt_atom::{Elastic, Friction, Material};
    use soil_gpu::neighbor_slots::PAYLOAD;

    fn entry(canonical: f64) -> [f64; CONTACT_HISTORY_LEN] {
        let mut history = [0.0; CONTACT_HISTORY_LEN];
        history[0] = canonical;
        history
    }

    fn table(rolling: f64, names: &[&str]) -> MaterialTable {
        let mut mt = MaterialTable::new();
        for name in names {
            let friction = Friction {
                sliding: 0.5,
                rolling,
                ..Default::default()
            };
            mt.add(Material::new(*name, Elastic::new(1e7, 0.25, 0.5)).with_friction(friction))
                .expect("valid material");
        }
        mt.build_pair_tables();
        mt
    }

    #[test]
    fn plain_hertz_mindlin_is_supported_and_others_are_named() {
        assert_eq!(
            unsupported_features(&table(0.0, &["m"])),
            Vec::<String>::new()
        );
        assert_eq!(
            unsupported_features(&table(0.1, &["a", "b"])),
            Vec::<String>::new()
        );
        let mut hooke = table(0.1, &["m"]);
        hooke.contact_model = "hooke".to_string();
        assert!(unsupported_features(&hooke)[0].contains("hooke"));
    }

    #[test]
    fn each_atom_gets_its_own_half_of_a_spring_recorded_on_one_side() {
        // Rows 0 and 1 are local, tags 5 and 9. Only tag 5 recorded their spring, and a
        // ghost of tag 3 recorded one with tag 9.
        let tag = [5, 9, 3];
        let mut history = ContactHistoryStore::new();
        history.contacts = vec![
            vec![(9, entry(2.0), false)],
            vec![],
            vec![(9, entry(4.0), false)],
        ];
        let (partner, payload) = store_to_slots(&history, &tag, 2, PAYLOAD);

        let at = |row: usize, k: usize| {
            (
                partner[row * MAX_CONTACTS + k],
                payload[(row * MAX_CONTACTS + k) * PAYLOAD],
            )
        };
        // Lower tag holds the canonical sign; the higher tag holds the opposite.
        assert_eq!(at(0, 0), (9, 2.0));
        assert_eq!(at(1, 0), (5, -2.0));
        assert_eq!(at(1, 1), (3, -4.0));
        assert_eq!(partner[2], SENTINEL);
    }

    #[test]
    fn a_rows_own_record_wins_over_a_partners() {
        let tag = [5, 9];
        let mut history = ContactHistoryStore::new();
        history.contacts = vec![vec![(9, entry(2.0), false)], vec![(5, entry(7.0), false)]];
        let (partner, payload) = store_to_slots(&history, &tag, 2, PAYLOAD);
        assert_eq!((partner[0], payload[0]), (9, 2.0));
        assert_eq!(
            (partner[MAX_CONTACTS], payload[MAX_CONTACTS * PAYLOAD]),
            (5, -7.0)
        );
        assert_eq!(partner[1], SENTINEL, "no duplicate slot");
    }
}
