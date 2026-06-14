//! Hertz-Mindlin contact force with the region-coherence optimization,
//! self-contained in this example (the stock `dirt_granular` crate is untouched).
//!
//! Per-particle modes are assigned by the quiescence module: `ACTIVE`,
//! `ADJACENT`, and `ASLEEP`. `ACTIVE` and `ADJACENT` particles are computed
//! normally — `ADJACENT` is the forced-on buffer ring around active cells, so
//! every contact of an active particle has at least one fully-computed end.
//! Only pairs where *both* particles are `ASLEEP` are skipped entirely (sleepers
//! are not integrated, so forces on them are not needed; a sleeper's force on an
//! awake neighbour is still computed because that pair has one awake member).

use dirt_core::prelude::*;
use dirt_core::dirt_atom::{self, DemAtom, MaterialTable};
use dirt_core::dirt_granular::{LARGE_OVERLAP_WARN_THRESHOLD, MAX_OVERLAP_WARNINGS, SQRT_5_6, TANGENTIAL_EPSILON};
use dirt_core::soil_core::{
    register_atom_data, Atom, AtomData, AtomDataRegistry, Neighbor, ParticleSimScheduleSet,
    VirialStress, VirialStressPlugin,
};
use std::any::Any;

/// Particle region modes set by the quiescence module each step.
/// `ACTIVE` (on) and `ADJACENT` (forced-on buffer next to active cells) are both
/// fully computed and integrated; `ASLEEP` (off) is skipped and pinned.
pub const MODE_ACTIVE: u8 = 0;
pub const MODE_ADJACENT: u8 = 1;
pub const MODE_ASLEEP: u8 = 2;

/// One contact's persistent tangential/rolling/twisting spring history.
#[derive(Clone, Copy, Debug)]
pub struct ContactEntry {
    /// Partner atom tag.
    pub tag: u32,
    /// Spring displacements: tangential [0..3], rolling [3..6], twisting [6].
    /// Canonical form (lower-tag perspective), as in the stock implementation.
    pub spring: [f64; 7],
    /// Touched this step (reset before each pair loop, pruned when false).
    pub active: bool,
}

impl ContactEntry {
    pub fn new(tag: u32) -> Self {
        ContactEntry {
            tag,
            spring: [0.0; 7],
            active: false,
        }
    }
}

/// Contact history + per-particle quiescence state. Registered as [`AtomData`]
/// so all per-atom vectors follow spatial-sort permutations and MPI exchange.
pub struct QcStore {
    /// Per-atom contact entries keyed by partner tag.
    pub contacts: Vec<Vec<ContactEntry>>,
    /// Per-atom region mode (MODE_ACTIVE / MODE_ADJACENT / MODE_ASLEEP).
    pub mode: Vec<u8>,
    /// Baseline force snapshot for sleeping particles (wake detector).
    pub f_base: Vec<[f64; 3]>,
    /// Whether `f_base` is valid for this particle.
    pub has_base: Vec<bool>,
    /// Consecutive quiescence checks this particle has been in force equilibrium
    /// (the sleep hysteresis counter). Per-atom so it survives spatial sorting.
    pub still_count: Vec<u16>,
    /// Per-step Σ|contact force| on each particle — the load scale used to
    /// normalize the unbalanced-force ratio. Recomputed every step, so it is not
    /// permuted/migrated (it is refilled before the quiescence check reads it).
    pub contact_mag: Vec<f64>,
    /// Window accumulators for the unbalanced-force ratio: running sum of the net
    /// force vector and of the (weight-floored) load scale over a check window.
    /// Per-atom and permuted on sort so the window follows each particle.
    pub accum_fnet: Vec<[f64; 3]>,
    pub accum_scale: Vec<f64>,

    // ── per-step knobs, set by the quiescence module ──
    /// Both-asleep pairs are skipped entirely.
    pub sleep_skip_on: bool,

    // ── per-step statistics (reset each step) ──
    pub n_pairs: usize,
    pub n_skipped: usize,
}

impl QcStore {
    pub fn new() -> Self {
        QcStore {
            contacts: Vec::new(),
            mode: Vec::new(),
            f_base: Vec::new(),
            has_base: Vec::new(),
            still_count: Vec::new(),
            contact_mag: Vec::new(),
            accum_fnet: Vec::new(),
            accum_scale: Vec::new(),
            sleep_skip_on: false,
            n_pairs: 0,
            n_skipped: 0,
        }
    }

    /// Grow per-particle vectors to `n` (new entries default to awake/active).
    pub fn ensure_len(&mut self, n: usize) {
        if self.contacts.len() < n {
            self.contacts.resize_with(n, Vec::new);
        }
        if self.mode.len() < n {
            self.mode.resize(n, MODE_ACTIVE);
        }
        if self.f_base.len() < n {
            self.f_base.resize(n, [0.0; 3]);
        }
        if self.has_base.len() < n {
            self.has_base.resize(n, false);
        }
        if self.still_count.len() < n {
            self.still_count.resize(n, 0);
        }
        if self.contact_mag.len() < n {
            self.contact_mag.resize(n, 0.0);
        }
        if self.accum_fnet.len() < n {
            self.accum_fnet.resize(n, [0.0; 3]);
        }
        if self.accum_scale.len() < n {
            self.accum_scale.resize(n, 0.0);
        }
    }
}

impl AtomData for QcStore {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn truncate(&mut self, n: usize) {
        self.contacts.resize_with(n, Vec::new);
        self.contacts.truncate(n);
        self.mode.resize(n, MODE_ACTIVE);
        self.mode.truncate(n);
        self.f_base.resize(n, [0.0; 3]);
        self.f_base.truncate(n);
        self.has_base.resize(n, false);
        self.has_base.truncate(n);
        self.still_count.resize(n, 0);
        self.still_count.truncate(n);
        self.contact_mag.resize(n, 0.0);
        self.contact_mag.truncate(n);
        self.accum_fnet.resize(n, [0.0; 3]);
        self.accum_fnet.truncate(n);
        self.accum_scale.resize(n, 0.0);
        self.accum_scale.truncate(n);
    }

    fn swap_remove(&mut self, i: usize) {
        if i < self.contacts.len() {
            self.contacts.swap_remove(i);
        }
        if i < self.mode.len() {
            self.mode.swap_remove(i);
        }
        if i < self.f_base.len() {
            self.f_base.swap_remove(i);
        }
        if i < self.has_base.len() {
            self.has_base.swap_remove(i);
        }
        if i < self.still_count.len() {
            self.still_count.swap_remove(i);
        }
        if i < self.contact_mag.len() {
            self.contact_mag.swap_remove(i);
        }
        if i < self.accum_fnet.len() {
            self.accum_fnet.swap_remove(i);
        }
        if i < self.accum_scale.len() {
            self.accum_scale.swap_remove(i);
        }
    }

    fn apply_permutation(&mut self, perm: &[usize], n: usize) {
        self.ensure_len(n);
        let new_contacts: Vec<Vec<ContactEntry>> =
            perm.iter().map(|&p| self.contacts[p].clone()).collect();
        self.contacts[..n].clone_from_slice(&new_contacts);
        let new_mode: Vec<u8> = perm.iter().map(|&p| self.mode[p]).collect();
        self.mode[..n].copy_from_slice(&new_mode);
        let new_base: Vec<[f64; 3]> = perm.iter().map(|&p| self.f_base[p]).collect();
        self.f_base[..n].copy_from_slice(&new_base);
        let new_has: Vec<bool> = perm.iter().map(|&p| self.has_base[p]).collect();
        for (k, v) in new_has.iter().enumerate() {
            self.has_base[k] = *v;
        }
        let new_still: Vec<u16> = perm.iter().map(|&p| self.still_count[p]).collect();
        self.still_count[..n].copy_from_slice(&new_still);
        // Window accumulators follow their atom so a sort mid-window doesn't
        // scramble the running sums. (contact_mag is refilled each step, so it
        // is not permuted.)
        let new_fnet: Vec<[f64; 3]> = perm.iter().map(|&p| self.accum_fnet[p]).collect();
        self.accum_fnet[..n].copy_from_slice(&new_fnet);
        let new_scale: Vec<f64> = perm.iter().map(|&p| self.accum_scale[p]).collect();
        self.accum_scale[..n].copy_from_slice(&new_scale);
    }

    fn pack(&self, i: usize, buf: &mut Vec<f64>) {
        // Migrating atoms drop their freeze cache and sleep state: the contact
        // recomputes and the region re-evaluates on the receiving side. Spring
        // history is preserved (same wire format intent as the stock store).
        if i < self.contacts.len() {
            let list = &self.contacts[i];
            buf.push(list.len() as f64);
            for e in list {
                buf.push(e.tag as f64);
                for s in e.spring {
                    buf.push(s);
                }
            }
        } else {
            buf.push(0.0);
        }
    }

    fn unpack(&mut self, buf: &[f64]) -> usize {
        let count = buf[0] as usize;
        let mut list = Vec::with_capacity(count);
        let mut pos = 1;
        for _ in 0..count {
            let mut e = ContactEntry::new(buf[pos] as u32);
            for k in 0..7 {
                e.spring[k] = buf[pos + 1 + k];
            }
            list.push(e);
            pos += 8;
        }
        self.contacts.push(list);
        self.mode.push(MODE_ACTIVE);
        self.f_base.push([0.0; 3]);
        self.has_base.push(false);
        self.still_count.push(0);
        self.contact_mag.push(0.0);
        self.accum_fnet.push([0.0; 3]);
        self.accum_scale.push(0.0);
        pos
    }
}

/// Plugin registering the [`QcStore`] and the quiescent contact force system.
/// Drop-in replacement for `HertzMindlinContactPlugin` (Hertz model only).
pub struct QuiescentContactPlugin;

impl Plugin for QuiescentContactPlugin {
    fn provides(&self) -> Vec<&str> {
        vec!["contact_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["dem_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(VirialStressPlugin);
        register_atom_data!(app, QcStore::new());
        app.add_update_system(
            quiescent_contact_force.label("hertz_mindlin_contact"),
            ParticleSimScheduleSet::Force,
        );
    }
}

/// Hertz-Mindlin contact force with frozen-contact caching and sleeping-pair
/// skipping. Mirrors `dirt_granular::contact::hertz_mindlin_contact_force`
/// (without bond exclusions, which this example does not use).
pub fn quiescent_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    material_table: Res<MaterialTable>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    let newton = neighbor.newton;
    let mut dem = registry.expect_mut::<DemAtom>("quiescent_contact_force");
    let mut store = registry.expect_mut::<QcStore>("quiescent_contact_force");
    let dt = atoms.dt;

    let natoms = atoms.len();
    store.ensure_len(natoms);
    store.n_pairs = 0;
    store.n_skipped = 0;
    for m in store.contact_mag[..natoms].iter_mut() {
        *m = 0.0;
    }

    let nlocal = atoms.nlocal as usize;
    let mut overlap_warnings = 0usize;

    let sleep_skip = store.sleep_skip_on;
    let virial_active = virial.as_ref().map(|v| v.active).unwrap_or(false);
    let vs = if newton { 1.0 } else { 0.5 };

    // Reset active flags before the pair loop. Sleeping atoms keep their
    // entries untouched: their both-asleep contacts are skipped below, and
    // pruning them would discard spring history needed on wake.
    for i in 0..nlocal {
        if sleep_skip && store.mode[i] == MODE_ASLEEP {
            continue;
        }
        for entry in &mut store.contacts[i] {
            entry.active = false;
        }
    }

    for (i, j) in neighbor.pairs(nlocal) {
        // ── Optimization B: skip pairs where both particles sleep ──────────
        if sleep_skip
            && store.mode[i] == MODE_ASLEEP
            && j < nlocal
            && store.mode[j] == MODE_ASLEEP
        {
            store.n_skipped += 1;
            continue;
        }

        // Skip same-body pairs (sub-spheres of the same rigid body)
        if dirt_atom::same_body(&dem, i, j) {
            continue;
        }

        let r1 = dem.radius[i];
        let r2 = dem.radius[j];

        let dx = atoms.pos[j][0] - atoms.pos[i][0];
        let dy = atoms.pos[j][1] - atoms.pos[i][1];
        let dz = atoms.pos[j][2] - atoms.pos[i][2];
        let dist_sq = dx * dx + dy * dy + dz * dz;
        let sum_r = r1 + r2;

        let mat_i = atoms.atom_type[i] as usize;
        let mat_j = atoms.atom_type[j] as usize;
        let surface_energy = material_table.surface_energy_ij[mat_i][mat_j];
        let use_dmt = material_table.adhesion_model == "dmt";

        let r_eff = (r1 * r2) / sum_r;
        let e_eff = material_table.e_eff_ij[mat_i][mat_j];
        let delta_pulloff = if surface_energy > 0.0 && !use_dmt {
            let gamma = surface_energy;
            (std::f64::consts::PI * std::f64::consts::PI * gamma * gamma * r_eff
                / (4.0 * e_eff * e_eff))
                .cbrt()
        } else {
            0.0
        };

        let interaction_r = sum_r + delta_pulloff;
        if dist_sq >= interaction_r * interaction_r {
            continue;
        }

        let distance = dist_sq.sqrt();
        if distance == 0.0 {
            continue;
        }

        let r_min = r1.min(r2);
        let delta = (sum_r - distance).min(0.5 * r_min);

        if delta > 0.0 && distance / sum_r < LARGE_OVERLAP_WARN_THRESHOLD {
            overlap_warnings += 1;
            if overlap_warnings > MAX_OVERLAP_WARNINGS {
                panic!(
                    "Over {} excessive overlaps this step — aborting. \
                     Check timestep or initial configuration.",
                    MAX_OVERLAP_WARNINGS
                );
            }
        }

        if delta <= 0.0 && surface_energy <= 0.0 {
            continue;
        }

        store.n_pairs += 1;

        let tag_i = atoms.tag[i];
        let tag_j = atoms.tag[j];
        let sign: f64 = if tag_i < tag_j { 1.0 } else { -1.0 };

        // Tangential spring-history lookup.
        let entry_idx = store.contacts[i].iter().position(|e| e.tag == tag_j);

        // ── Full Hertz-Mindlin model (mirrors the stock implementation) ─────
        let inv_dist = 1.0 / distance;
        let nx = dx * inv_dist;
        let ny = dy * inv_dist;
        let nz = dz * inv_dist;

        let g_eff = material_table.g_eff_ij[mat_i][mat_j];
        let inv_m_i = if atoms.inv_mass[i] > 0.0 { atoms.inv_mass[i] } else { 1.0 / atoms.mass[i] };
        let inv_m_j = if atoms.inv_mass[j] > 0.0 { atoms.inv_mass[j] } else { 1.0 / atoms.mass[j] };
        let m_r = 1.0 / (inv_m_i + inv_m_j);

        let beta = material_table.beta_ij[mat_i][mat_j];
        let mu = material_table.friction_ij[mat_i][mat_j];
        let mu_r = material_table.rolling_friction_ij[mat_i][mat_j];
        let mu_tw = material_table.twisting_friction_ij[mat_i][mat_j];
        let cohesion_energy = material_table.cohesion_energy_ij[mat_i][mat_j];

        let jkr_adhesion_only = surface_energy > 0.0 && !use_dmt && delta <= 0.0;

        let (s_n, k_n, k_t) = if delta > 0.0 {
            let sdr = (delta * r_eff).sqrt();
            let sn = 2.0 * e_eff * sdr;
            let kn = 4.0 / 3.0 * e_eff * sdr;
            let kt = 8.0 * g_eff * sdr;
            (sn, kn, kt)
        } else {
            (0.0, 0.0, 0.0)
        };

        let omega_ix = dem.omega[i][0];
        let omega_iy = dem.omega[i][1];
        let omega_iz = dem.omega[i][2];
        let omega_jx = dem.omega[j][0];
        let omega_jy = dem.omega[j][1];
        let omega_jz = dem.omega[j][2];

        let r1n_x = r1 * nx;
        let r1n_y = r1 * ny;
        let r1n_z = r1 * nz;
        let vc_ix = atoms.vel[i][0] + (omega_iy * r1n_z - omega_iz * r1n_y);
        let vc_iy = atoms.vel[i][1] + (omega_iz * r1n_x - omega_ix * r1n_z);
        let vc_iz = atoms.vel[i][2] + (omega_ix * r1n_y - omega_iy * r1n_x);

        let r2n_x = r2 * nx;
        let r2n_y = r2 * ny;
        let r2n_z = r2 * nz;
        let vc_jx = atoms.vel[j][0] + (-omega_jy * r2n_z + omega_jz * r2n_y);
        let vc_jy = atoms.vel[j][1] + (-omega_jz * r2n_x + omega_jx * r2n_z);
        let vc_jz = atoms.vel[j][2] + (-omega_jx * r2n_y + omega_jy * r2n_x);

        let vr_x = vc_jx - vc_ix;
        let vr_y = vc_jy - vc_iy;
        let vr_z = vc_jz - vc_iz;

        let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

        // ── Normal force ─────────────────────────────────────────────────
        let f_n_mag = if surface_energy > 0.0 && use_dmt {
            let f_dmt = 2.0 * std::f64::consts::PI * surface_energy * r_eff;
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            k_n * delta - f_diss_n - f_dmt
        } else if surface_energy > 0.0 {
            let f_adhesion = 1.5 * std::f64::consts::PI * surface_energy * r_eff;
            if jkr_adhesion_only {
                -f_adhesion
            } else {
                let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
                k_n * delta - f_diss_n - f_adhesion
            }
        } else if cohesion_energy > 0.0 {
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * r_eff;
            k_n * delta - f_diss_n - f_cohesion
        } else {
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            (k_n * delta - f_diss_n).max(0.0)
        };

        let fn_x = f_n_mag * nx;
        let fn_y = f_n_mag * ny;
        let fn_z = f_n_mag * nz;

        atoms.force[i][0] -= fn_x;
        atoms.force[i][1] -= fn_y;
        atoms.force[i][2] -= fn_z;
        if newton {
            atoms.force[j][0] += fn_x;
            atoms.force[j][1] += fn_y;
            atoms.force[j][2] += fn_z;
        }

        // ── Tangential (skip in JKR adhesion-only regime) ─────────────────
        if jkr_adhesion_only {
            if virial_active {
                if let Some(ref mut v) = virial {
                    v.add_pair(dx, dy, dz, -fn_x * vs, -fn_y * vs, -fn_z * vs);
                }
            }
            continue;
        }

        let vt_x = vr_x - v_n * nx;
        let vt_y = vr_y - v_n * ny;
        let vt_z = vr_z - v_n * nz;

        let stored = match entry_idx {
            Some(idx) => store.contacts[i][idx].spring,
            None => [0.0; 7],
        };

        let mut sx = sign * stored[0];
        let mut sy = sign * stored[1];
        let mut sz = sign * stored[2];
        let s_dot_n = sx * nx + sy * ny + sz * nz;
        sx -= s_dot_n * nx;
        sy -= s_dot_n * ny;
        sz -= s_dot_n * nz;
        sx += vt_x * dt;
        sy += vt_y * dt;
        sz += vt_z * dt;

        let s_mag = (sx * sx + sy * sy + sz * sz).sqrt();
        let f_t_spring_mag = k_t * s_mag;
        let f_t_max = mu * f_n_mag.abs();
        if f_t_spring_mag > f_t_max && f_t_spring_mag > TANGENTIAL_EPSILON {
            let scale = f_t_max / f_t_spring_mag;
            sx *= scale;
            sy *= scale;
            sz *= scale;
        }

        let gamma_t = 2.0 * SQRT_5_6 * beta * (k_t * m_r).sqrt();
        // NOTE: deliberate sign fix vs the stock implementation. With
        // vt = tangential velocity of j relative to i and +ft applied to i,
        // the dashpot must act ALONG vt to dissipate (P = -ft·vt < 0).
        // Stock uses `- gamma_t * vt`, which injects energy during sliding
        // and keeps dense beds fluidized indefinitely (the stock hopper
        // example never reaches its KE settling threshold because of this).
        let mut ft_x = k_t * sx + gamma_t * vt_x;
        let mut ft_y = k_t * sy + gamma_t * vt_y;
        let mut ft_z = k_t * sz + gamma_t * vt_z;

        let f_t_mag = (ft_x * ft_x + ft_y * ft_y + ft_z * ft_z).sqrt();
        if f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
            let scale = f_t_max / f_t_mag;
            ft_x *= scale;
            ft_y *= scale;
            ft_z *= scale;
        }

        atoms.force[i][0] += ft_x;
        atoms.force[i][1] += ft_y;
        atoms.force[i][2] += ft_z;
        if newton {
            atoms.force[j][0] -= ft_x;
            atoms.force[j][1] -= ft_y;
            atoms.force[j][2] -= ft_z;
        }

        // Per-pair torque accumulators (needed for the freeze cache).
        let mut tau_i = [
            r1n_y * ft_z - r1n_z * ft_y,
            r1n_z * ft_x - r1n_x * ft_z,
            r1n_x * ft_y - r1n_y * ft_x,
        ];
        let mut tau_j = [
            r2n_y * ft_z - r2n_z * ft_y,
            r2n_z * ft_x - r2n_x * ft_z,
            r2n_x * ft_y - r2n_y * ft_x,
        ];

        // ── Rolling resistance ─────────────────────────────────────────────
        let or_x = omega_ix - omega_jx;
        let or_y = omega_iy - omega_jy;
        let or_z = omega_iz - omega_jz;
        let or_dot_n = or_x * nx + or_y * ny + or_z * nz;
        let roll_x = or_x - or_dot_n * nx;
        let roll_y = or_y - or_dot_n * ny;
        let roll_z = or_z - or_dot_n * nz;

        let mut roll_disp_x = sign * stored[3];
        let mut roll_disp_y = sign * stored[4];
        let mut roll_disp_z = sign * stored[5];
        let mut twist_disp = sign * stored[6];

        if mu_r > 0.0 {
            let roll_mag = (roll_x * roll_x + roll_y * roll_y + roll_z * roll_z).sqrt();
            let sds_rolling = material_table.rolling_model == "sds";
            if sds_rolling {
                let k_roll = material_table.rolling_stiffness_ij[mat_i][mat_j];
                let gamma_roll = material_table.rolling_damping_ij[mat_i][mat_j];

                let rd_dot_n = roll_disp_x * nx + roll_disp_y * ny + roll_disp_z * nz;
                roll_disp_x -= rd_dot_n * nx;
                roll_disp_y -= rd_dot_n * ny;
                roll_disp_z -= rd_dot_n * nz;
                roll_disp_x += roll_x * dt;
                roll_disp_y += roll_y * dt;
                roll_disp_z += roll_z * dt;

                let mut tr_x = -k_roll * roll_disp_x - gamma_roll * roll_x;
                let mut tr_y = -k_roll * roll_disp_y - gamma_roll * roll_y;
                let mut tr_z = -k_roll * roll_disp_z - gamma_roll * roll_z;
                let tr_mag = (tr_x * tr_x + tr_y * tr_y + tr_z * tr_z).sqrt();
                let tau_max = mu_r * f_n_mag.abs() * r_eff;

                if tr_mag > tau_max && tr_mag > TANGENTIAL_EPSILON {
                    let scale = tau_max / tr_mag;
                    tr_x *= scale;
                    tr_y *= scale;
                    tr_z *= scale;
                    if k_roll > TANGENTIAL_EPSILON {
                        roll_disp_x = (tr_x + gamma_roll * roll_x) / (-k_roll);
                        roll_disp_y = (tr_y + gamma_roll * roll_y) / (-k_roll);
                        roll_disp_z = (tr_z + gamma_roll * roll_z) / (-k_roll);
                    }
                }

                tau_i[0] += tr_x;
                tau_i[1] += tr_y;
                tau_i[2] += tr_z;
                tau_j[0] -= tr_x;
                tau_j[1] -= tr_y;
                tau_j[2] -= tr_z;
            } else if roll_mag > 1e-30 {
                let tau_mag = mu_r * f_n_mag.abs() * r_eff;
                let inv_roll = tau_mag / roll_mag;
                let tr_x = -inv_roll * roll_x;
                let tr_y = -inv_roll * roll_y;
                let tr_z = -inv_roll * roll_z;
                tau_i[0] += tr_x;
                tau_i[1] += tr_y;
                tau_i[2] += tr_z;
                tau_j[0] -= tr_x;
                tau_j[1] -= tr_y;
                tau_j[2] -= tr_z;
            }
        }

        // ── Twisting friction ──────────────────────────────────────────────
        if mu_tw > 0.0 {
            let twist_vel = or_dot_n;
            let sds_twisting = material_table.twisting_model == "sds";
            if sds_twisting {
                let k_twist = material_table.twisting_stiffness_ij[mat_i][mat_j];
                let gamma_twist = material_table.twisting_damping_ij[mat_i][mat_j];

                twist_disp += twist_vel * dt;

                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_tw * f_n_mag.abs() * r_eff;

                if tau_twist.abs() > tau_max {
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                tau_i[0] += tau_twist * nx;
                tau_i[1] += tau_twist * ny;
                tau_i[2] += tau_twist * nz;
                tau_j[0] -= tau_twist * nx;
                tau_j[1] -= tau_twist * ny;
                tau_j[2] -= tau_twist * nz;
            } else if twist_vel.abs() > 1e-30 {
                let tau = mu_tw * f_n_mag.abs() * r_eff;
                let sign_tw = if twist_vel > 0.0 { -1.0 } else { 1.0 };
                tau_i[0] += sign_tw * tau * nx;
                tau_i[1] += sign_tw * tau * ny;
                tau_i[2] += sign_tw * tau * nz;
                tau_j[0] -= sign_tw * tau * nx;
                tau_j[1] -= sign_tw * tau * ny;
                tau_j[2] -= sign_tw * tau * nz;
            }
        }

        dem.torque[i][0] += tau_i[0];
        dem.torque[i][1] += tau_i[1];
        dem.torque[i][2] += tau_i[2];
        if newton {
            dem.torque[j][0] += tau_j[0];
            dem.torque[j][1] += tau_j[1];
            dem.torque[j][2] += tau_j[2];
        }

        // Total pair force on i (normal + tangential), for virial and the
        // unbalanced-force-ratio load scale.
        let fp_x = -fn_x + ft_x;
        let fp_y = -fn_y + ft_y;
        let fp_z = -fn_z + ft_z;

        // Accumulate |contact force| on each end (the load scale). Newton: the
        // force on j has equal magnitude, so add to both; otherwise the (j,i)
        // pair will add j's own contribution.
        let fp_mag = (fp_x * fp_x + fp_y * fp_y + fp_z * fp_z).sqrt();
        store.contact_mag[i] += fp_mag;
        if newton {
            store.contact_mag[j] += fp_mag;
        }

        if virial_active {
            if let Some(ref mut v) = virial {
                v.add_pair(dx, dy, dz, fp_x * vs, fp_y * vs, fp_z * vs);
            }
        }

        // ── Write back tangential/rolling/twisting spring history ─────────
        let new_spring = [
            sign * sx,
            sign * sy,
            sign * sz,
            sign * roll_disp_x,
            sign * roll_disp_y,
            sign * roll_disp_z,
            sign * twist_disp,
        ];
        let idx = match entry_idx {
            Some(idx) => idx,
            None => {
                store.contacts[i].push(ContactEntry::new(tag_j));
                store.contacts[i].len() - 1
            }
        };
        let e = &mut store.contacts[i][idx];
        e.spring = new_spring;
        e.active = true;
    }

    // Prune stale contacts (skip sleeping atoms — their lists are preserved).
    for i in 0..nlocal {
        if sleep_skip && store.mode[i] == MODE_ASLEEP {
            continue;
        }
        store.contacts[i].retain(|e| e.active);
    }
}
