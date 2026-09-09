//! Per-atom wall spring history, carried by SOIL's `AtomData` contract.
//!
//! # Why this is an `AtomData` and not a map on the [`Walls`](crate::Walls) resource
//!
//! Until `dirt-wall-springs-atomdata` the tangential and rolling spring
//! histories lived in two `HashMap<(u8, usize, u32), [f64; 3]>` fields on
//! `Walls`. `Walls` is a **resource**, and a resource is rank-local: SOIL's
//! atom exchange never sees it. A particle in sustained frictional contact
//! with a wall that crossed a rank boundary therefore arrived on the receiving
//! rank with no entry for its tag, took the `unwrap_or([0.0; 3])` default, and
//! **lost its whole tangential spring**.
//!
//! That is not a small error. `examples/mpi_wall_springs/RESULTS.md` measured
//! it on the shipped native path with no kernel anywhere: 15 crossings in a
//! 2-rank run, every one of them in wall contact, each discharging **100 %**
//! of the tangential force (worst `|dF_t|` = 1.196797e-4 N) and taking
//! **1352 steps** to ring back down, ~168 of them above 10 % of the force. The
//! normal force agreed to exactly zero on all 90 000 compared `(step, tag)`
//! rows, so the decomposition itself was sound — the spring history was the
//! only thing being dropped.
//!
//! The particle–particle history never had this defect:
//! `dirt_granular::ContactHistoryStore` is a `soil_core::AtomData` with its own
//! `pack`/`unpack`, so it travels inside the exchange message. This store is
//! the same pattern applied to the wall history, which is what makes the wall
//! case the rule rather than the anomaly.
//!
//! # Keying, and why it needs no cross-rank resolution
//!
//! The old key was `(wall_kind, wall_index, atom_tag)`. Splitting it is what
//! makes the history a per-atom row: the **atom** half becomes the row index,
//! and the **wall** half stays inside the row as `(kind, index)`.
//!
//! `(kind, index)` is already a *global* identifier and needs no handshake.
//! Every rank builds `Walls` in `WallPlugin::build` from the same `[[wall]]`
//! array of the same TOML file, appending to `planes`, `cylinders`, `spheres`
//! and `regions` in parse order, so a given wall has the same `(kind, index)`
//! on every rank. Runtime `active` flags do not shift the numbering — a
//! deactivated wall keeps its slot.
//!
//! # Row lifetime
//!
//! Row `i` belongs to local atom `i`, exactly as `Atom::pos[i]` does. SOIL
//! keeps it there: [`AtomData::apply_permutation`] follows a spatial sort,
//! [`AtomData::swap_remove`] follows a deletion, and
//! [`AtomData::pack`]/[`AtomData::unpack`] carry the row inside the migration
//! message. Nothing in this crate has to re-associate anything by tag.
//!
//! # Pruning
//!
//! The old code rebuilt both maps from scratch every step, so a contact that
//! ended was pruned by simply not being re-inserted. The equivalent here is
//! [`WallSpringStore::begin_step`] (clear every entry's `active` flag) and
//! [`WallSpringStore::prune`] (drop the entries nothing touched), called
//! around the wall force loop. The observable semantics are identical: a
//! contact that ends loses its spring, and a contact with no history starts
//! from zero displacement.

use std::any::Any;

use soil_core::AtomData;

/// Scalars stored per `(atom, wall)` entry in the migration wire format:
/// three for the tangential spring, three for the SDS rolling displacement.
pub const WALL_SPRING_LEN: usize = 6;

/// One wall's spring history for one atom.
///
/// `kind` is `0 = plane, 1 = cylinder, 2 = sphere, 3 = region` and `index` is
/// the position within that list — the same `(wall_kind, wall_index)` the
/// old `Walls` maps used in their key, and global for the reason given in the
/// [module documentation](self).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallSpring {
    /// Wall kind: 0 = plane, 1 = cylinder, 2 = sphere, 3 = region.
    pub kind: u8,
    /// Index within the list selected by [`kind`](Self::kind).
    pub index: u32,
    /// Mindlin tangential spring displacement (m).
    pub tangential: [f64; 3],
    /// SDS rolling spring displacement (m). Zero under the default `constant`
    /// rolling model, which is stateless.
    pub rolling: [f64; 3],
    /// Touched by this step's wall force loop. Cleared by
    /// [`WallSpringStore::begin_step`] and used by
    /// [`WallSpringStore::prune`] to drop contacts that have ended. Never
    /// transmitted: an unpacked entry arrives inactive and is kept only if the
    /// receiving rank's own force loop touches it this step.
    pub active: bool,
}

/// Per-atom wall spring history: one row per atom, one entry per wall the atom
/// is (or was, this step) in contact with.
///
/// Registered in the [`AtomDataRegistry`](soil_core::AtomDataRegistry) by
/// [`WallPlugin`](crate::WallPlugin), which is what makes the history travel
/// with a migrating atom.
#[derive(Default)]
pub struct WallSpringStore {
    /// `rows[i]` is local atom `i`'s wall spring history.
    pub rows: Vec<Vec<WallSpring>>,
}

impl WallSpringStore {
    /// An empty store with no atom rows.
    pub fn new() -> Self {
        WallSpringStore { rows: Vec::new() }
    }

    /// Grow the store to `n` rows, so every atom — including ghosts, which
    /// never hold entries but must keep the row indices aligned — has one.
    pub fn ensure_rows(&mut self, n: usize) {
        if self.rows.len() < n {
            self.rows.resize_with(n, Vec::new);
        }
    }

    /// Mark every entry of the first `nlocal` rows untouched, before a force
    /// pass. Entries the pass does not touch are dropped by [`Self::prune`].
    pub fn begin_step(&mut self, nlocal: usize) {
        let n = nlocal.min(self.rows.len());
        for row in &mut self.rows[..n] {
            for entry in row.iter_mut() {
                entry.active = false;
            }
        }
    }

    /// Drop every entry of the first `nlocal` rows that this step's force pass
    /// did not touch — the contacts that have ended.
    pub fn prune(&mut self, nlocal: usize) {
        let n = nlocal.min(self.rows.len());
        for row in &mut self.rows[..n] {
            row.retain(|entry| entry.active);
        }
    }

    fn find(&self, atom: usize, kind: u8, index: u32) -> Option<&WallSpring> {
        self.rows
            .get(atom)?
            .iter()
            .find(|e| e.kind == kind && e.index == index)
    }

    /// This atom's tangential spring against `(kind, index)`, or zeros when
    /// there is no history — the same default the old `HashMap` lookup took.
    pub fn tangential(&self, atom: usize, kind: u8, index: u32) -> [f64; 3] {
        self.find(atom, kind, index)
            .map(|e| e.tangential)
            .unwrap_or([0.0; 3])
    }

    /// This atom's SDS rolling displacement against `(kind, index)`, or zeros.
    pub fn rolling(&self, atom: usize, kind: u8, index: u32) -> [f64; 3] {
        self.find(atom, kind, index)
            .map(|e| e.rolling)
            .unwrap_or([0.0; 3])
    }

    /// Find or create this atom's entry for `(kind, index)` and mark it
    /// touched this step.
    fn entry_mut(&mut self, atom: usize, kind: u8, index: u32) -> &mut WallSpring {
        self.ensure_rows(atom + 1);
        let row = &mut self.rows[atom];
        let at = match row.iter().position(|e| e.kind == kind && e.index == index) {
            Some(at) => at,
            None => {
                row.push(WallSpring {
                    kind,
                    index,
                    tangential: [0.0; 3],
                    rolling: [0.0; 3],
                    active: false,
                });
                row.len() - 1
            }
        };
        let entry = &mut row[at];
        entry.active = true;
        entry
    }

    /// Store this step's tangential spring for `(atom, kind, index)`.
    pub fn set_tangential(&mut self, atom: usize, kind: u8, index: u32, spring: [f64; 3]) {
        self.entry_mut(atom, kind, index).tangential = spring;
    }

    /// Store this step's SDS rolling displacement for `(atom, kind, index)`.
    pub fn set_rolling(&mut self, atom: usize, kind: u8, index: u32, spring: [f64; 3]) {
        self.entry_mut(atom, kind, index).rolling = spring;
    }

    /// Every entry of one atom's row, in insertion order.
    pub fn row(&self, atom: usize) -> &[WallSpring] {
        self.rows.get(atom).map(|r| r.as_slice()).unwrap_or(&[])
    }

    /// Total number of live entries across the first `nlocal` rows. Used by
    /// tests and diagnostics to assert that a fixture is actually holding
    /// springs.
    pub fn live_entries(&self, nlocal: usize) -> usize {
        self.rows[..nlocal.min(self.rows.len())]
            .iter()
            .map(|r| r.len())
            .sum()
    }
}

impl AtomData for WallSpringStore {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn snapshot(&self) -> Box<dyn AtomData> {
        Box::new(WallSpringStore {
            rows: self.rows.clone(),
        })
    }

    fn len(&self) -> usize {
        self.rows.len()
    }

    unsafe fn push_default(&mut self) {
        self.rows.push(Vec::new());
    }

    unsafe fn truncate(&mut self, n: usize) {
        // Grow first: atoms can be inserted without going through `unpack`.
        self.rows.resize_with(n, Vec::new);
        self.rows.truncate(n);
    }

    unsafe fn swap_remove(&mut self, i: usize) {
        if i < self.rows.len() {
            self.rows.swap_remove(i);
        }
    }

    unsafe fn apply_permutation(&mut self, perm: &[usize], n: usize) {
        let permuted: Vec<Vec<WallSpring>> =
            perm.iter().map(|&p| self.rows[p].clone()).collect();
        self.rows[..n].clone_from_slice(&permuted);
    }

    /// Wire format: `count`, then `count` records of
    /// `[kind, index, tangential x y z, rolling x y z]`.
    ///
    /// `active` is deliberately not transmitted. It is a within-step marker,
    /// and an arriving entry must earn it back from the receiving rank's own
    /// force loop or be pruned — which is exactly the "a contact that ended
    /// loses its spring" rule, applied on the rank that can now see the
    /// contact.
    fn pack(&self, i: usize, buf: &mut Vec<f64>) {
        match self.rows.get(i) {
            Some(row) => {
                buf.push(row.len() as f64);
                for entry in row {
                    buf.push(entry.kind as f64);
                    buf.push(entry.index as f64);
                    buf.extend_from_slice(&entry.tangential);
                    buf.extend_from_slice(&entry.rolling);
                }
            }
            None => buf.push(0.0),
        }
    }

    unsafe fn unpack(&mut self, buf: &[f64]) -> usize {
        let count = buf[0] as usize;
        let mut row = Vec::with_capacity(count);
        let mut pos = 1;
        // 2 identity scalars + WALL_SPRING_LEN spring scalars per entry.
        let stride = 2 + WALL_SPRING_LEN;
        for _ in 0..count {
            let kind = buf[pos] as u8;
            let index = buf[pos + 1] as u32;
            let mut tangential = [0.0; 3];
            let mut rolling = [0.0; 3];
            tangential.copy_from_slice(&buf[pos + 2..pos + 5]);
            rolling.copy_from_slice(&buf[pos + 5..pos + 8]);
            row.push(WallSpring {
                kind,
                index,
                tangential,
                rolling,
                active: false,
            });
            pos += stride;
        }
        self.rows.push(row);
        pos
    }
}
