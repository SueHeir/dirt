//! Per-particle force-balance quiescence: freeze particles that are in static
//! equilibrium, keep a one-contact-shell buffer around anything still moving.
//! Self-contained in this example — no cell grid.
//!
//! Three per-particle modes, assigned every `check_every` steps:
//!
//! - **ACTIVE** (on): not in force equilibrium — |net force| above a fraction
//!   `force_tol_frac` of the particle's weight (m·g). A free-falling or
//!   being-shoved particle is never balanced, so it can never freeze.
//! - **ADJACENT** (buffer): in contact-range of an ACTIVE particle. Forced fully
//!   on (and integrated) so every contact of an active particle has at least one
//!   computed end; distinctly labelled so it can sleep again once activity moves
//!   away.
//! - **ASLEEP** (off): in equilibrium for `still_checks` consecutive checks and
//!   not adjacent to activity. Velocity/force zeroed (pinned); pairs where *both*
//!   ends are asleep are skipped in the contact loop.
//!
//! Wake-up: every step, each sleeper snapshots the partial force it receives from
//! gravity, walls, and *awake* neighbours; a relative drift > `wake_force_rel`
//! wakes it. Removing a support (blocker wall, discharged neighbour) changes that
//! partial force and trips the wake; the adjacency buffer then propagates the
//! front through the contact network.

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;
use serde::Deserialize;

use crate::contact::{QcStore, MODE_ACTIVE, MODE_ADJACENT, MODE_ASLEEP};

/// TOML `[quiescence]` section.
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields, default)]
pub struct QuiescenceConfig {
    /// Enable the optimization.
    pub region_coherence: bool,
    /// Classification cadence [steps].
    pub check_every: usize,
    /// Consecutive in-equilibrium checks before a particle freezes.
    pub still_checks: usize,
    /// Equilibrium threshold: the unbalanced-force ratio — window-averaged
    /// |net force| divided by the window-summed contact-force scale (a particle's
    /// own weight is used as a floor). Below this the particle is "balanced".
    /// Dimensionless and load-independent (interior force-chain particles and
    /// light surface particles use the same cut). A free-flyer's ratio is ~1.
    pub force_ratio_tol: f64,
    /// Relative force deviation that wakes a sleeping particle.
    pub wake_force_rel: f64,
    /// Thickness of the forced-on ADJACENT buffer around active particles, in
    /// contact-shells (neighbour-list hops). 1 = direct neighbours only; larger
    /// values widen the simulated margin around activity (more "green").
    pub buffer_layers: usize,
}

impl Default for QuiescenceConfig {
    fn default() -> Self {
        QuiescenceConfig {
            region_coherence: false,
            check_every: 25,
            still_checks: 3,
            force_ratio_tol: 0.1,
            wake_force_rel: 0.3,
            buffer_layers: 1,
        }
    }
}

/// Per-particle quiescence state machine (no grid).
pub struct Quiescence {
    pub coh_on: bool,
    check_every: usize,
    still_checks: u16,
    force_ratio_tol: f64,
    wake_force_rel: f64,
    buffer_layers: usize,
    /// Gravitational acceleration magnitude, the contact-scale floor (= weight).
    g: f64,
    // Per-atom classification scratch, recomputed each check (sized to nlocal).
    // The window accumulators themselves live in QcStore so they follow atoms
    // through spatial sorting.
    active_flag: Vec<bool>,
    adj_flag: Vec<bool>,
    /// Scratch list of particles newly added to the buffer in one expansion layer.
    adj_frontier: Vec<usize>,
}

impl Quiescence {
    pub fn new() -> Self {
        Quiescence {
            coh_on: false,
            check_every: 25,
            still_checks: 3,
            force_ratio_tol: 0.1,
            wake_force_rel: 0.3,
            buffer_layers: 1,
            g: 9.81,
            active_flag: Vec::new(),
            adj_flag: Vec::new(),
            adj_frontier: Vec::new(),
        }
    }
}

/// Plugin wiring the quiescence systems around the contact loop.
pub struct QuiescencePlugin;

impl Plugin for QuiescencePlugin {
    fn build(&self, app: &mut App) {
        Config::load::<QuiescenceConfig>(app, "quiescence");
        app.add_resource(Quiescence::new());
        app.add_setup_system(q_setup, ScheduleSetupSet::PostSetup);
        app.add_update_system(q_post_force, ParticleSimScheduleSet::PostForce);
        app.add_update_system(q_zero_sleepers, ParticleSimScheduleSet::PreFinalIntegration);

        // Expose the per-particle quiescence mode (0 = active, 1 = adjacent,
        // 2 = asleep) as a dump/VTP column — color by it in OVITO/ParaView to
        // see frozen regions, the adjacent buffer, and the wake-front. Uses the
        // DumpRegistry extension point, from this (non-core) plugin's build().
        if let Some(dump_reg) = app.get_resource_ref::<DumpRegistry>() {
            dump_reg.register_scalar("mode", |atoms, registry| {
                let nlocal = atoms.nlocal as usize;
                match registry.get::<QcStore>() {
                    Some(store) => (0..nlocal)
                        .map(|i| if i < store.mode.len() { store.mode[i] as f64 } else { 0.0 })
                        .collect(),
                    None => vec![0.0; nlocal],
                }
            });
            // Binary "frozen" flag (1 = asleep/pinned, 0 = active or adjacent) —
            // convenient for thresholding/coloring in ParaView or OVITO.
            dump_reg.register_scalar("frozen", |atoms, registry| {
                let nlocal = atoms.nlocal as usize;
                match registry.get::<QcStore>() {
                    Some(store) => (0..nlocal)
                        .map(|i| {
                            let frozen = i < store.mode.len() && store.mode[i] == MODE_ASLEEP;
                            frozen as u8 as f64
                        })
                        .collect(),
                    None => vec![0.0; nlocal],
                }
            });
        }
    }
}

/// Read config + gravity, push the per-contact knob into the [`QcStore`].
fn q_setup(
    config: Res<QuiescenceConfig>,
    mut q: ResMut<Quiescence>,
    gravity: Res<GravityConfig>,
    registry: Res<AtomDataRegistry>,
    comm: Res<CommResource>,
) {
    q.coh_on = config.region_coherence;
    q.check_every = config.check_every.max(1);
    q.still_checks = config.still_checks.max(1) as u16;
    q.force_ratio_tol = config.force_ratio_tol;
    q.wake_force_rel = config.wake_force_rel;
    q.buffer_layers = config.buffer_layers.max(1);
    q.g = (gravity.gx * gravity.gx + gravity.gy * gravity.gy + gravity.gz * gravity.gz).sqrt();

    if comm.size() > 1 && q.coh_on {
        // The prototype's wake logic is single-rank (no ghost mode exchange).
        if comm.rank() == 0 {
            println!("Quiescence: DISABLED — prototype supports single-rank runs only");
        }
        q.coh_on = false;
    }

    if let Some(mut store) = registry.get_mut::<QcStore>() {
        store.sleep_skip_on = q.coh_on;
    }

    if comm.rank() == 0 && q.coh_on {
        println!(
            "Quiescence: per-particle unbalanced-force ratio, g={:.3}, force_ratio_tol={}, still_checks={}, buffer_layers={}, check_every={}",
            q.g, q.force_ratio_tol, q.still_checks, q.buffer_layers, q.check_every
        );
    }
}

/// Wake sleeping particles on force deviation (every step) and run the
/// equilibrium classification (every `check_every` steps).
fn q_post_force(
    mut atoms: ResMut<Atom>,
    mut q: ResMut<Quiescence>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
) {
    let q = &mut *q;
    if !q.coh_on {
        return;
    }
    let mut store = registry.expect_mut::<QcStore>("q_post_force");
    let mut dem = registry.expect_mut::<DemAtom>("q_post_force");
    let nlocal = atoms.nlocal as usize;
    store.ensure_len(atoms.len());
    if q.active_flag.len() < nlocal {
        q.active_flag.resize(nlocal, false);
        q.adj_flag.resize(nlocal, false);
    }

    // ── 1. Every step: wake sleepers on force drift; for awake particles,
    //    accumulate this window's net force and (weight-floored) load scale. ──
    // A sleeper's force is partial (both-asleep pairs are skipped) but steady
    // while the sleeping set is unchanged; a drift means a neighbour woke/moved
    // or a support vanished → wake.
    let mut sleep_set_changed = false;
    let wake_sq = q.wake_force_rel * q.wake_force_rel;
    for i in 0..nlocal {
        let f = atoms.force[i];
        if store.mode[i] == MODE_ASLEEP {
            if !store.has_base[i] {
                store.f_base[i] = f;
                store.has_base[i] = true;
            } else {
                let b = store.f_base[i];
                let dev_sq =
                    (f[0] - b[0]).powi(2) + (f[1] - b[1]).powi(2) + (f[2] - b[2]).powi(2);
                let bmag_sq = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
                if dev_sq > wake_sq * bmag_sq {
                    store.mode[i] = MODE_ACTIVE;
                    store.has_base[i] = false;
                    store.still_count[i] = 0;
                    atoms.force[i] = [0.0; 3];
                    sleep_set_changed = true;
                }
            }
        } else {
            store.accum_fnet[i][0] += f[0];
            store.accum_fnet[i][1] += f[1];
            store.accum_fnet[i][2] += f[2];
            let weight = atoms.mass[i] * q.g;
            store.accum_scale[i] += store.contact_mag[i].max(weight);
        }
    }

    // ── 2. Classification via the windowed unbalanced-force ratio ──────────
    if run_state.total_cycle % q.check_every == 0 && nlocal > 0 {
        // (a) ACTIVE = awake and unbalanced. The window's mean net force exceeds
        //     force_ratio_tol × its mean load scale. The per-step counts cancel,
        //     so compare |Σ f_net| directly against tol · Σ scale.
        let tol = q.force_ratio_tol;
        for i in 0..nlocal {
            if store.mode[i] == MODE_ASLEEP {
                q.active_flag[i] = false;
                continue;
            }
            let a = store.accum_fnet[i];
            let fnet_sq = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
            let lim = tol * store.accum_scale[i];
            q.active_flag[i] = fnet_sq > lim * lim;
        }

        // (b) ADJACENT buffer: forced-on shell(s) around active particles, grown
        //     outward `buffer_layers` neighbour-hops. Each pass marks neighbours
        //     of the current simulated set (active ∪ adjacent) that aren't in it
        //     yet, staged into a frontier so one pass adds exactly one shell.
        for a in q.adj_flag[..nlocal].iter_mut() {
            *a = false;
        }
        for _layer in 0..q.buffer_layers {
            q.adj_frontier.clear();
            for (i, j) in neighbor.pairs(nlocal) {
                let i_in = q.active_flag[i] || q.adj_flag[i];
                let j_in = j < nlocal && (q.active_flag[j] || q.adj_flag[j]);
                if i_in && j < nlocal && !j_in {
                    q.adj_frontier.push(j);
                }
                if j_in && !i_in {
                    q.adj_frontier.push(i);
                }
            }
            if q.adj_frontier.is_empty() {
                break;
            }
            for k in 0..q.adj_frontier.len() {
                q.adj_flag[q.adj_frontier[k]] = true;
            }
        }

        // (c) Apply, then reset the window accumulators.
        for i in 0..nlocal {
            if q.active_flag[i] {
                store.mode[i] = MODE_ACTIVE;
                store.still_count[i] = 0;
            } else if q.adj_flag[i] {
                if store.mode[i] == MODE_ASLEEP {
                    store.has_base[i] = false;
                    sleep_set_changed = true;
                }
                store.mode[i] = MODE_ADJACENT;
                store.still_count[i] = 0;
            } else if store.mode[i] != MODE_ASLEEP {
                // Awake, balanced, away from activity: accumulate, then freeze.
                store.still_count[i] = store.still_count[i].saturating_add(1);
                if store.still_count[i] >= q.still_checks {
                    store.mode[i] = MODE_ASLEEP;
                    store.has_base[i] = false;
                    atoms.vel[i] = [0.0; 3];
                    dem.omega[i] = [0.0; 3];
                    sleep_set_changed = true;
                } else {
                    store.mode[i] = MODE_ACTIVE;
                }
            }
            store.accum_fnet[i] = [0.0; 3];
            store.accum_scale[i] = 0.0;
        }
    }

    // Sleeping-set changes invalidate every partial-force baseline.
    if sleep_set_changed {
        for i in 0..nlocal {
            store.has_base[i] = false;
        }
    }
}

/// Pin-style integration skip: zero force/velocity (and torque/omega) of
/// sleeping particles so both Verlet kicks and the rotational update are no-ops.
fn q_zero_sleepers(mut atoms: ResMut<Atom>, q: Res<Quiescence>, registry: Res<AtomDataRegistry>) {
    if !q.coh_on {
        return;
    }
    let store = registry.expect_mut::<QcStore>("q_zero_sleepers");
    let mut dem = registry.expect_mut::<DemAtom>("q_zero_sleepers");
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal {
        if i < store.mode.len() && store.mode[i] == MODE_ASLEEP {
            atoms.force[i] = [0.0; 3];
            atoms.vel[i] = [0.0; 3];
            dem.torque[i] = [0.0; 3];
            dem.omega[i] = [0.0; 3];
        }
    }
}
