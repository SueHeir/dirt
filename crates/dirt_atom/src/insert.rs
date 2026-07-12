//! Particle insertion: random, rate-based, and file-based from `[[particles.insert]]` config.
//!
//! # Parallel insertion model (born-in-owner)
//!
//! Insertion is **not** done on one rank and scattered. Instead it runs on
//! *every* MPI rank, and every rank seeds its RNG with the **same** `seed`, so
//! all ranks generate the bit-identical candidate stream — the same positions,
//! radii, velocities, and tags, in the same order. Each rank then keeps **only**
//! the candidates whose position falls inside its own subdomain, tested with a
//! half-open interval (`low ≤ pos < high`, via the `owns_position` helper) that matches
//! `Domain::exchange()` ownership exactly.
//!
//! Three properties follow:
//!
//! - **Born in its owner.** Each atom is materialized only on the rank that owns
//!   its position, so the first post-insertion `exchange()` never has to migrate
//!   it — no insertion-time communication is needed.
//! - **Exactly-once.** The half-open convention guarantees every position is
//!   claimed by one rank, never two and never zero (no duplicated or dropped
//!   atoms at subdomain boundaries).
//! - **Deterministic across rank counts.** Because the candidate stream depends
//!   only on `seed` (not on the decomposition), the *global* packing is
//!   identical whether you run on 1 rank or 64 — only the partitioning of that
//!   packing changes. Overlap rejection uses a global spatial hash replicated on
//!   every rank, so packings are reproducible run-to-run and rank-count to
//!   rank-count.
//!
//! File-based insertion follows the same rule: the file is parsed on every rank
//! and filtered to the local subdomain, with tags advanced identically so they
//! stay globally consistent.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use serde::Deserialize;

use grass_scheduler::prelude::CurrentState;
use soil_core::{
    deep_merge, Atom, AtomDataRegistry, CommResource, CommState, Config, Domain, DomainConfig,
    ParticleSimScheduleSet, ParticleStore, Real, Region, RunConfig, RunState, ScheduleSetupSet,
    StageOverrides,
};

use crate::{DemAtom, MaterialTable, RadiusSpec};

// ── Particle insertion ─────────────────────────────────────────────────────

mod configuration;
mod immediate;
mod rate;
mod readers;
mod sampling;
mod setup;

use configuration::*;
use immediate::*;
use readers::*;
use sampling::*;

pub use configuration::{
    ColumnMapping, InsertConfig, ParticlesConfig, RateInsertEntry, RateInsertState,
};
pub use immediate::dem_insert_atoms;
pub use rate::dem_rate_insert;
pub use setup::DemAtomInsertPlugin;

#[cfg(test)]
mod tests;
