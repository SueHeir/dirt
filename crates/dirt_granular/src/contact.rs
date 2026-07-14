//! Fused Hertz-Mindlin contact force computation.
//!
//! The public contact API is assembled from responsibility-focused modules. The
//! force loops retain the established expression order while dispatch, shared
//! state, normal laws, adhesion, MDR support, and regression tests stay local to
//! their respective implementations.

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use dirt_atom::{self, DemAtom, MaterialTable};
use dirt_schedule::CONTACT_FORCE;
use soil_core::Neighbor;
use soil_core::{forward_comm_overlap, CommBuffers, CommResource, CommTopology};
use soil_core::{
    register_atom_data, Atom, AtomDataRegistry, BondStore, ParticleSimScheduleSet, VirialStress,
    VirialStressPlugin,
};

use crate::tangential::{ContactHistory, ContactHistoryStore, CONTACT_HISTORY_LEN};
use crate::{LARGE_OVERLAP_WARN_THRESHOLD, MAX_OVERLAP_WARNINGS, SQRT_5_6, TANGENTIAL_EPSILON};

mod adhesion;
mod force;
mod mdr;
mod normal;
mod plugin;
mod shared;

use mdr::mdr_normal_force;
use normal::{mdr_nonadhesive_force, mdr_round_up_negative_epsilon};
use shared::{is_mindlin_force_history, is_mindlin_rescale, zero_contact_history};

pub use adhesion::willett2000_liquid_bridge_force;
pub use force::{
    contact_force_core, hertz_mindlin_contact_force, hooke_contact_force, overlapped_contact_force,
};
pub use plugin::HertzMindlinContactPlugin;
pub use shared::ForcePass;

#[cfg(test)]
mod tests;
