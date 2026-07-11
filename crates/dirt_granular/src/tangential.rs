//! Standalone Mindlin tangential contact force with spring history.
//!
//! This module provides the tangential friction force as a separate plugin,
//! for use when normal and tangential are registered independently. For the
//! recommended fused implementation, use [`crate::contact::HertzMindlinContactPlugin`].
//!
//! # Physics
//!
//! The Mindlin no-slip tangential force uses an incremental spring displacement:
//!
//! 1. Rotate previous spring into current tangent plane (remove normal component)
//! 2. Integrate: `s += v_t · dt`
//! 3. Optionally rescale the history on normal unloading for the
//!    `mindlin_rescale` variants
//! 4. Cap spring at Coulomb limit: `|k_t s| ≤ μ |F_n|`
//! 5. Tangential force: `F_t = k_t s - γ_t v_t`, capped at `μ |F_n|`
//!
//! where `k_t = 8 G* √(R* δ)` is the tangential stiffness, `G*` is the effective
//! shear modulus, and `γ_t = 2 β √(5/6) √(k_t m_r)` is the tangential damping.
//!
//! Spring displacements are stored in **canonical form** (from the perspective of the
//! particle with the lower tag) so that both particles in a pair see a consistent history
//! regardless of which is `i` vs `j` in the neighbor list.

use std::any::Any;

use soil_core::AtomData;

// ── ContactHistoryStore ─────────────────────────────────────────────────────

/// Number of scalar history values stored for each active contact.
pub const CONTACT_HISTORY_LEN: usize = 23;

/// Fixed-width per-contact history payload; see [`ContactHistoryStore`] for slot meanings.
pub type ContactHistory = [f64; CONTACT_HISTORY_LEN];

/// Per-contact spring displacement history for tangential, rolling, and twisting models.
///
/// Each contact entry is stored in **canonical form** — from the perspective of the
/// particle with the lower tag — so the spring is frame-consistent regardless of
/// neighbor list ordering. A `sign` factor of `+1` or `-1` converts between the
/// canonical frame and the local `(i, j)` frame each timestep.
///
/// # Storage layout
///
/// Each contact stores [`CONTACT_HISTORY_LEN`] `f64` displacement/history values:
///
/// | Indices | Model              | Description                             |
/// |---------|--------------------|-----------------------------------------|
/// | `[0..3]`| Mindlin tangential | Tangential displacement, or elastic force for `/force` variants |
/// | `[3..6]`| SDS rolling        | Rolling spring displacement vector      |
/// | `[6]`  | SDS twisting       | Twisting spring displacement (scalar)   |
/// | `[7]`  | Mindlin rescale    | Previous contact radius `a = sqrt(R* delta)` |
/// | `[8..]`| MDR normal         | LAMMPS-style per-side rigid-flat history |
///
/// Rolling and twisting slots are zero when the constant-torque model is used.

/// Per-contact spring and MDR normal history store.
pub struct ContactHistoryStore {
    /// Per-atom list of `(partner_tag, spring_displacement/history, active_flag)`.
    ///
    /// `active_flag` is reset to `false` before each pair loop and set to `true`
    /// when a contact is touched. Stale entries (broken contacts) are pruned after
    /// the loop completes.
    pub contacts: Vec<Vec<(u32, ContactHistory, bool)>>,
}

impl ContactHistoryStore {
    /// Create an empty contact history store with no pre-allocated atoms.
    pub fn new() -> Self {
        ContactHistoryStore {
            contacts: Vec::new(),
        }
    }
}

impl AtomData for ContactHistoryStore {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn len(&self) -> usize {
        self.contacts.len()
    }

    fn push_default(&mut self) {
        self.contacts.push(Vec::new());
    }

    fn truncate(&mut self, n: usize) {
        // Grow if needed (atoms may have been inserted without going through unpack)
        self.contacts.resize_with(n, Vec::new);
        self.contacts.truncate(n);
    }

    fn swap_remove(&mut self, i: usize) {
        if i < self.contacts.len() {
            self.contacts.swap_remove(i);
        }
    }

    fn apply_permutation(&mut self, perm: &[usize], n: usize) {
        let new_contacts: Vec<Vec<(u32, ContactHistory, bool)>> =
            perm.iter().map(|&p| self.contacts[p].clone()).collect();
        self.contacts[..n].clone_from_slice(&new_contacts);
    }

    fn pack(&self, i: usize, buf: &mut Vec<f64>) {
        if i < self.contacts.len() {
            let list = &self.contacts[i];
            buf.push(list.len() as f64);
            for &(tag, ref s, _) in list {
                buf.push(tag as f64);
                for value in s {
                    buf.push(*value);
                }
            }
        } else {
            buf.push(0.0); // no contacts
        }
    }

    fn unpack(&mut self, buf: &[f64]) -> usize {
        let count = buf[0] as usize;
        let mut list = Vec::with_capacity(count);
        let mut pos = 1;
        let old_seven_slot_format = buf.len() == 1 + count * 8;
        let old_eight_slot_format = buf.len() == 1 + count * 9;
        for _ in 0..count {
            let tag = buf[pos] as u32;
            let mut s = [0.0; CONTACT_HISTORY_LEN];
            if old_seven_slot_format {
                s[..7].copy_from_slice(&buf[pos + 1..pos + 8]);
                pos += 8;
            } else if old_eight_slot_format {
                s[..8].copy_from_slice(&buf[pos + 1..pos + 9]);
                pos += 9;
            } else {
                s[..CONTACT_HISTORY_LEN]
                    .copy_from_slice(&buf[pos + 1..pos + 1 + CONTACT_HISTORY_LEN]);
                pos += 1 + CONTACT_HISTORY_LEN;
            }
            list.push((tag, s, false));
        }
        self.contacts.push(list);
        pos
    }
}
