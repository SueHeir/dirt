//! Wall contact-force kernel.

use crate::geometry::Walls;
use dirt_atom::{DemAtom, MaterialTable, SQRT_5_6};
use grass_scheduler::prelude::*;
use soil_core::{Accum, Atom, ParticlesWith, Write};

fn wall_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Mindlin tangential (sliding) friction for a single particle–wall contact,
/// mirroring the particle–particle Hertz–Mindlin model in `dirt_granular`.
///
/// `n` is the unit contact normal pointing from the wall surface toward the
/// particle center; `v_rel` is the particle velocity minus the wall velocity;
/// `old_spring` is the accumulated tangential displacement from the previous
/// step. Returns `(force_on_particle, torque_on_particle, new_spring)`. Returns
/// zeros when `mu <= 0` or `delta <= 0`, so frictionless walls are byte-for-byte
/// unchanged and a purely normal impact yields no tangential force.
#[allow(clippy::too_many_arguments)]
fn wall_tangential_force(
    n: [f64; 3],
    v_rel: [f64; 3],
    omega: [f64; 3],
    radius: f64,
    delta: f64,
    f_n: f64,
    m_r: f64,
    g_eff: f64,
    mu: f64,
    beta: f64,
    dt: f64,
    old_spring: [f64; 3],
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    if mu <= 0.0 || delta <= 0.0 {
        return ([0.0; 3], [0.0; 3], [0.0; 3]);
    }
    // Particle surface velocity relative to the wall at the contact point
    // (contact point is at r_c = -radius * n from the particle center).
    let oxn = wall_cross(omega, n);
    let vs = [
        v_rel[0] - radius * oxn[0],
        v_rel[1] - radius * oxn[1],
        v_rel[2] - radius * oxn[2],
    ];
    let vsn = vs[0] * n[0] + vs[1] * n[1] + vs[2] * n[2];
    let vt = [vs[0] - vsn * n[0], vs[1] - vsn * n[1], vs[2] - vsn * n[2]];

    // Mindlin tangential stiffness; R* = radius (wall is flat / infinite radius).
    let k_t = 8.0 * g_eff * (radius * delta).sqrt();

    // Advance the stored spring: project out the normal component, then
    // accumulate the tangential relative displacement vt*dt.
    let mut s = old_spring;
    let sn = s[0] * n[0] + s[1] * n[1] + s[2] * n[2];
    s = [
        s[0] - sn * n[0] + vt[0] * dt,
        s[1] - sn * n[1] + vt[1] * dt,
        s[2] - sn * n[2] + vt[2] * dt,
    ];

    let ft_max = mu * f_n.abs();
    let s_mag = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    if k_t * s_mag > ft_max && k_t * s_mag > 1e-30 {
        let c = ft_max / (k_t * s_mag);
        s = [s[0] * c, s[1] * c, s[2] * c];
    }

    let gamma_t = 2.0 * SQRT_5_6 * beta * (k_t * m_r).sqrt();
    // Friction opposes the particle's relative tangential surface motion.
    let mut ft = [
        -(k_t * s[0] + gamma_t * vt[0]),
        -(k_t * s[1] + gamma_t * vt[1]),
        -(k_t * s[2] + gamma_t * vt[2]),
    ];
    let ft_mag = (ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2]).sqrt();
    if ft_mag > ft_max && ft_mag > 1e-30 {
        let c = ft_max / ft_mag;
        ft = [ft[0] * c, ft[1] * c, ft[2] * c];
    }

    // Torque about the particle center: r_c × ft with r_c = -radius * n.
    let nxf = wall_cross(n, ft);
    let tau = [-radius * nxf[0], -radius * nxf[1], -radius * nxf[2]];
    (ft, tau, s)
}

/// Rolling-resistance torque for a particle–wall contact, mirroring the
/// particle–particle rolling model in `dirt_granular` with the wall as a
/// zero-spin second body (so the rolling angular velocity is just the
/// particle's spin with its normal/twisting component removed).
///
/// Supports both the `constant` (fixed-magnitude couple opposing the roll) and
/// `sds` (spring–dashpot–slider, Coulomb-capped) models.  SDS follows the
/// LAMMPS pseudo-force convention: it stores a length, computes a force, then
/// converts it to a couple. `r_eff = radius` for a flat wall. Returns
/// `(torque_on_particle, new_rolling_displacement)`; the
/// displacement is only meaningful (and stored) for the `sds` model. Returns
/// zeros when `mu_r <= 0`.
#[allow(clippy::too_many_arguments)]
fn wall_rolling_torque(
    n: [f64; 3],
    omega: [f64; 3],
    radius: f64,
    f_n: f64,
    mu_r: f64,
    k_roll: f64,
    gamma_roll: f64,
    sds: bool,
    dt: f64,
    old_roll_disp: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    if mu_r <= 0.0 {
        return ([0.0; 3], [0.0; 3]);
    }
    // Rolling angular velocity = particle spin minus its normal (twisting) part.
    let odn = omega[0] * n[0] + omega[1] * n[1] + omega[2] * n[2];
    let roll = [
        omega[0] - odn * n[0],
        omega[1] - odn * n[1],
        omega[2] - odn * n[2],
    ];
    let tau_max = mu_r * f_n.abs() * radius;

    if sds {
        // LAMMPS GranularModel: v_rl = R_eff (omega_rel x n).  The history is
        // therefore a displacement in metres and k_roll/gamma_roll have the
        // pseudo-force units N/m and N s/m.
        let vrl = [
            radius * (omega[1] * n[2] - omega[2] * n[1]),
            radius * (omega[2] * n[0] - omega[0] * n[2]),
            radius * (omega[0] * n[1] - omega[1] * n[0]),
        ];
        let mut rd = old_roll_disp;
        let rdn = rd[0] * n[0] + rd[1] * n[1] + rd[2] * n[2];
        rd = [
            rd[0] - rdn * n[0] + vrl[0] * dt,
            rd[1] - rdn * n[1] + vrl[1] * dt,
            rd[2] - rdn * n[2] + vrl[2] * dt,
        ];
        let mut fr = [
            -k_roll * rd[0] - gamma_roll * vrl[0],
            -k_roll * rd[1] - gamma_roll * vrl[1],
            -k_roll * rd[2] - gamma_roll * vrl[2],
        ];
        let fr_mag = (fr[0] * fr[0] + fr[1] * fr[1] + fr[2] * fr[2]).sqrt();
        let fr_max = mu_r * f_n.abs();
        if fr_mag > fr_max && fr_mag > 1e-30 {
            let s = fr_max / fr_mag;
            fr = [fr[0] * s, fr[1] * s, fr[2] * s];
            if k_roll > 1e-30 {
                rd = [
                    (fr[0] + gamma_roll * vrl[0]) / (-k_roll),
                    (fr[1] + gamma_roll * vrl[1]) / (-k_roll),
                    (fr[2] + gamma_roll * vrl[2]) / (-k_roll),
                ];
            }
        }
        let torque = [
            radius * (n[1] * fr[2] - n[2] * fr[1]),
            radius * (n[2] * fr[0] - n[0] * fr[2]),
            radius * (n[0] * fr[1] - n[1] * fr[0]),
        ];
        (torque, rd)
    } else {
        // Constant-torque model: fixed magnitude opposing the rolling direction.
        let roll_mag = (roll[0] * roll[0] + roll[1] * roll[1] + roll[2] * roll[2]).sqrt();
        if roll_mag > 1e-30 {
            let inv = tau_max / roll_mag;
            ([-inv * roll[0], -inv * roll[1], -inv * roll[2]], [0.0; 3])
        } else {
            ([0.0; 3], [0.0; 3])
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn wall_normal_force(
    material_table: &MaterialTable,
    mat_i: usize,
    wall_mat: usize,
    radius: f64,
    delta: f64,
    v_n: f64,
    m_r: f64,
    allow_plane_surface_energy: bool,
) -> f64 {
    let beta = material_table.beta_ij[mat_i][wall_mat];
    let cohesion_energy = material_table.cohesion_energy_ij[mat_i][wall_mat];

    if material_table.contact_model == "hooke" {
        let kn = material_table.kn_ij[mat_i][wall_mat];
        let gamma_n = 2.0 * beta * (kn * m_r).sqrt();
        let f_total = if cohesion_energy > 0.0 {
            let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * radius;
            kn * delta - gamma_n * v_n - f_cohesion
        } else {
            kn * delta - gamma_n * v_n
        };
        if material_table.limit_damping && cohesion_energy <= 0.0 {
            f_total.max(0.0)
        } else {
            f_total
        }
    } else {
        let r_eff = radius;
        let e_eff = material_table.e_eff_ij[mat_i][wall_mat];
        let sdr = (delta * r_eff).sqrt();
        let k_n = 4.0 / 3.0 * e_eff * sdr;
        let s_n = 2.0 * e_eff * sdr;
        let surface_energy = if allow_plane_surface_energy {
            material_table.surface_energy_ij[mat_i][wall_mat]
        } else {
            0.0
        };
        let use_dmt = material_table.adhesion_model == "dmt";

        if surface_energy > 0.0 && use_dmt {
            let f_dmt = 2.0 * std::f64::consts::PI * surface_energy * r_eff;
            let f_diss = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            k_n * delta - f_diss - f_dmt
        } else if surface_energy > 0.0 {
            let f_adhesion = 1.5 * std::f64::consts::PI * surface_energy * r_eff;
            let f_diss = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            k_n * delta - f_diss - f_adhesion
        } else if cohesion_energy > 0.0 {
            let f_diss = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * r_eff;
            k_n * delta - f_diss - f_cohesion
        } else {
            let f_diss = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_total = k_n * delta - f_diss;
            if material_table.limit_damping {
                f_total.max(0.0)
            } else {
                f_total
            }
        }
    }
}

/// Constant twisting-friction torque for a particle-wall contact.
///
/// The local contact normal points from the wall surface toward the particle
/// center. Twisting friction opposes the particle spin component about that
/// normal and uses the same cap as the existing plane-wall model,
/// `tau = mu_tw * |F_n| * R*` with `R* = particle_radius`.
fn wall_twisting_torque(
    n: [f64; 3],
    omega: [f64; 3],
    radius: f64,
    f_n: f64,
    mu_tw: f64,
) -> [f64; 3] {
    if mu_tw <= 0.0 {
        return [0.0; 3];
    }

    let twist = omega[0] * n[0] + omega[1] * n[1] + omega[2] * n[2];
    if twist.abs() <= 1e-30 {
        return [0.0; 3];
    }

    let tau = mu_tw * f_n.abs() * radius;
    let sign_tw = if twist > 0.0 { -1.0 } else { 1.0 };
    [
        sign_tw * tau * n[0],
        sign_tw * tau * n[1],
        sign_tw * tau * n[2],
    ]
}

/// System that applies wall–particle contact forces for every registered wall.
///
/// For each local atom in contact with a wall it evaluates the Hertzian normal
/// force with viscous damping, the tangential/rolling/twisting friction springs
/// (whose per-contact history is carried in [`Walls`]), and any adhesion model,
/// accumulating the result into the atom's force and torque. Tangential and
/// rolling spring histories are rebuilt each step so that contacts which ended
/// are pruned automatically.
pub fn wall_contact_force(
    mut atoms: ResMut<Atom>,
    mut walls: ResMut<Walls>,
    particles: ParticlesWith<'_, Write<DemAtom>>,
    material_table: Res<MaterialTable>,
) {
    particles.with(|mut dem| {
        let nlocal = atoms.nlocal as usize;
        let dt = atoms.dt;

        // Take the tangential-spring history out so the per-wall-list immutable
        // borrows below don't conflict with mutating it; rebuild it fresh this step
        // (contacts that ended are pruned by not being re-inserted).
        let old_springs = std::mem::take(&mut walls.tangential_springs);
        let mut new_springs: std::collections::HashMap<(u8, usize, u32), [f64; 3]> =
            std::collections::HashMap::new();
        let old_rolling = std::mem::take(&mut walls.rolling_springs);
        let mut new_rolling: std::collections::HashMap<(u8, usize, u32), [f64; 3]> =
            std::collections::HashMap::new();

        // Collect per-wall forces to accumulate after the loop
        let nwalls = walls.planes.len();
        let mut wall_forces = vec![0.0f64; nwalls];

        for (wall_idx, wall) in walls.planes.iter().enumerate() {
            if !walls.active[wall_idx] {
                continue;
            }

            let wall_mat = wall.material_index;

            for i in 0..nlocal {
                let px = atoms.pos[i][0] as f64;
                let py = atoms.pos[i][1] as f64;
                let pz = atoms.pos[i][2] as f64;

                // Check if atom is within the wall's bounding region
                if !wall.in_bounds(px, py, pz) {
                    continue;
                }

                // Vector from wall point to atom position
                let dx = px - wall.point_x;
                let dy = py - wall.point_y;
                let dz = pz - wall.point_z;

                // Signed distance from atom to wall plane (positive = on normal side)
                let distance = dx * wall.normal_x + dy * wall.normal_y + dz * wall.normal_z;

                // Only apply force when atom center is on the normal side of the wall
                if distance <= 0.0 {
                    continue;
                }

                let radius = dem.radius[i];
                let mat_i = atoms.atom_type[i] as usize;

                let surface_energy = material_table.surface_energy_ij[mat_i][wall_mat];
                let use_hertz = material_table.contact_model != "hooke";
                let use_dmt = material_table.adhesion_model == "dmt";

                // JKR pull-off distance for extended interaction range
                // DMT: no extended range (particles separate at delta = 0)
                let delta_pulloff = if use_hertz && surface_energy > 0.0 && !use_dmt {
                    let r_eff = radius;
                    let e_eff = material_table.e_eff_ij[mat_i][wall_mat];
                    let gamma = surface_energy;
                    (std::f64::consts::PI * std::f64::consts::PI * gamma * gamma * r_eff
                        / (4.0 * e_eff * e_eff))
                        .cbrt()
                } else {
                    0.0
                };

                let delta = (radius - distance).min(0.5 * radius);

                // Skip if no contact and no JKR adhesion range
                if delta <= 0.0 && (!use_hertz || surface_energy <= 0.0) {
                    continue;
                }
                if delta < -delta_pulloff {
                    continue;
                }

                // JKR adhesion-only regime; DMT has no adhesion-only regime
                let jkr_adhesion_only =
                    use_hertz && surface_energy > 0.0 && !use_dmt && delta <= 0.0;

                // Wall has infinite mass → m_reduced = m_particle
                let m_r = atoms.mass[i] as f64;

                // Relative velocity along wall normal (subtract wall velocity)
                let v_rel_x = atoms.vel[i][0] as f64 - wall.velocity[0];
                let v_rel_y = atoms.vel[i][1] as f64 - wall.velocity[1];
                let v_rel_z = atoms.vel[i][2] as f64 - wall.velocity[2];
                let v_n =
                    v_rel_x * wall.normal_x + v_rel_y * wall.normal_y + v_rel_z * wall.normal_z;

                let f_net = if jkr_adhesion_only {
                    let f_adhesion = 1.5 * std::f64::consts::PI * surface_energy * radius;
                    -f_adhesion
                } else {
                    wall_normal_force(
                        &material_table,
                        mat_i,
                        wall_mat,
                        radius,
                        delta,
                        v_n,
                        m_r,
                        true,
                    )
                };

                // Force direction: along wall normal (pushes atom away from wall)
                atoms.force[i][0] += (f_net * wall.normal_x) as Accum;
                atoms.force[i][1] += (f_net * wall.normal_y) as Accum;
                atoms.force[i][2] += (f_net * wall.normal_z) as Accum;

                // Twisting friction torque (wall-particle)
                if delta > 0.0 {
                    let mu_tw = material_table.twisting_friction_ij[mat_i][wall_mat];
                    let tt = wall_twisting_torque(
                        [wall.normal_x, wall.normal_y, wall.normal_z],
                        dem.omega[i],
                        radius,
                        f_net,
                        mu_tw,
                    );
                    dem.torque[i][0] += tt[0];
                    dem.torque[i][1] += tt[1];
                    dem.torque[i][2] += tt[2];
                }

                // Tangential (Mindlin) sliding friction.
                let mu = material_table.friction_ij[mat_i][wall_mat];
                if mu > 0.0 && delta > 0.0 {
                    let beta = material_table.beta_ij[mat_i][wall_mat];
                    let g_eff = material_table.g_eff_ij[mat_i][wall_mat];
                    let n = [wall.normal_x, wall.normal_y, wall.normal_z];
                    let v_rel = [v_rel_x, v_rel_y, v_rel_z];
                    let key = (0u8, wall_idx, atoms.tag[i]);
                    let old = old_springs.get(&key).copied().unwrap_or([0.0; 3]);
                    let (ft, tau, ns) = wall_tangential_force(
                        n,
                        v_rel,
                        dem.omega[i],
                        radius,
                        delta,
                        f_net,
                        m_r,
                        g_eff,
                        mu,
                        beta,
                        dt,
                        old,
                    );
                    atoms.force[i][0] += ft[0] as Accum;
                    atoms.force[i][1] += ft[1] as Accum;
                    atoms.force[i][2] += ft[2] as Accum;
                    dem.torque[i][0] += tau[0];
                    dem.torque[i][1] += tau[1];
                    dem.torque[i][2] += tau[2];
                    new_springs.insert(key, ns);
                }

                // Rolling-resistance torque.
                let mu_r = material_table.rolling_friction_ij[mat_i][wall_mat];
                if mu_r > 0.0 && delta > 0.0 {
                    let sds = material_table.rolling_model == "sds";
                    let k_roll = material_table.rolling_stiffness_ij[mat_i][wall_mat];
                    let gamma_roll = material_table.rolling_damping_ij[mat_i][wall_mat];
                    let key = (0u8, wall_idx, atoms.tag[i]);
                    let old_rd = old_rolling.get(&key).copied().unwrap_or([0.0; 3]);
                    let (tr, new_rd) = wall_rolling_torque(
                        [wall.normal_x, wall.normal_y, wall.normal_z],
                        dem.omega[i],
                        radius,
                        f_net,
                        mu_r,
                        k_roll,
                        gamma_roll,
                        sds,
                        dt,
                        old_rd,
                    );
                    dem.torque[i][0] += tr[0];
                    dem.torque[i][1] += tr[1];
                    dem.torque[i][2] += tr[2];
                    if sds {
                        new_rolling.insert(key, new_rd);
                    }
                }

                // Accumulate wall force for servo control
                wall_forces[wall_idx] += f_net;
            }
        }

        // Write accumulated forces back to walls
        for (idx, &f) in wall_forces.iter().enumerate() {
            walls.planes[idx].force_accumulator += f;
        }

        // ── Cylinder walls ──────────────────────────────────────────────────
        let ncyl = walls.cylinders.len();
        let mut cyl_forces = vec![0.0f64; ncyl];
        for (cyl_idx, cyl) in walls.cylinders.iter().enumerate() {
            if !walls.cylinder_active[cyl_idx] {
                continue;
            }
            let wall_mat = cyl.material_index;
            for i in 0..nlocal {
                let pos = [
                    atoms.pos[i][0] as f64,
                    atoms.pos[i][1] as f64,
                    atoms.pos[i][2] as f64,
                ];
                let radius = dem.radius[i];
                let mat_i = atoms.atom_type[i] as usize;

                // Decompose position into axial and radial components.
                // `axial` = coordinate along the cylinder axis.
                // `d0, d1` = displacement from cylinder center in the 2D cross-section plane.
                // For axis=Z: axial=z, d0=x-cx, d1=y-cy.
                let (axial, d0, d1) = match cyl.axis {
                    0 => (pos[0], pos[1] - cyl.center[0], pos[2] - cyl.center[1]),
                    1 => (pos[1], pos[0] - cyl.center[0], pos[2] - cyl.center[1]),
                    _ => (pos[2], pos[0] - cyl.center[0], pos[1] - cyl.center[1]),
                };

                // Check axial bounds
                if axial < cyl.lo || axial > cyl.hi {
                    continue;
                }

                let radial_dist = (d0 * d0 + d1 * d1).sqrt();
                if radial_dist < 1e-30 {
                    continue;
                }

                // Compute overlap (delta) and 3D contact normal (nx, ny, nz).
                // The gap is the distance from the particle center to the wall surface.
                // The normal always points from the wall surface toward the particle center.
                let inv_r = 1.0 / radial_dist;
                let (delta, nx, ny, nz) = if cyl.inside {
                    // Inside: gap = cylinder_radius - radial_distance
                    let gap = cyl.radius - radial_dist;
                    let delta = radius - gap;
                    // Normal points inward (toward axis), pushing particle away from wall
                    let (n0, n1) = (-d0 * inv_r, -d1 * inv_r);
                    let (nx, ny, nz) = match cyl.axis {
                        0 => (0.0, n0, n1),
                        1 => (n0, 0.0, n1),
                        _ => (n0, n1, 0.0),
                    };
                    (delta, nx, ny, nz)
                } else {
                    // Outside: gap = radial_distance - cylinder_radius
                    let gap = radial_dist - cyl.radius;
                    let delta = radius - gap;
                    // Normal points outward (away from axis), pushing particle away from wall
                    let (n0, n1) = (d0 * inv_r, d1 * inv_r);
                    let (nx, ny, nz) = match cyl.axis {
                        0 => (0.0, n0, n1),
                        1 => (n0, 0.0, n1),
                        _ => (n0, n1, 0.0),
                    };
                    (delta, nx, ny, nz)
                };

                if delta <= 0.0 {
                    continue;
                }
                let delta = delta.min(0.5 * radius);

                let m_r = atoms.mass[i] as f64;
                let v_n = atoms.vel[i][0] as f64 * nx
                    + atoms.vel[i][1] as f64 * ny
                    + atoms.vel[i][2] as f64 * nz;
                let f_net = wall_normal_force(
                    &material_table,
                    mat_i,
                    wall_mat,
                    radius,
                    delta,
                    v_n,
                    m_r,
                    false,
                );

                atoms.force[i][0] += (f_net * nx) as Accum;
                atoms.force[i][1] += (f_net * ny) as Accum;
                atoms.force[i][2] += (f_net * nz) as Accum;

                // Twisting friction torque (cylinder wall is static).
                let mu_tw = material_table.twisting_friction_ij[mat_i][wall_mat];
                if mu_tw > 0.0 {
                    let tt = wall_twisting_torque([nx, ny, nz], dem.omega[i], radius, f_net, mu_tw);
                    dem.torque[i][0] += tt[0];
                    dem.torque[i][1] += tt[1];
                    dem.torque[i][2] += tt[2];
                }

                // Tangential (Mindlin) sliding friction (cylinder wall is static).
                let mu = material_table.friction_ij[mat_i][wall_mat];
                if mu > 0.0 {
                    let beta = material_table.beta_ij[mat_i][wall_mat];
                    let g_eff = material_table.g_eff_ij[mat_i][wall_mat];
                    let key = (1u8, cyl_idx, atoms.tag[i]);
                    let old = old_springs.get(&key).copied().unwrap_or([0.0; 3]);
                    let (ft, tau, ns) = wall_tangential_force(
                        [nx, ny, nz],
                        [
                            atoms.vel[i][0] as f64,
                            atoms.vel[i][1] as f64,
                            atoms.vel[i][2] as f64,
                        ],
                        dem.omega[i],
                        radius,
                        delta,
                        f_net,
                        m_r,
                        g_eff,
                        mu,
                        beta,
                        dt,
                        old,
                    );
                    atoms.force[i][0] += ft[0] as Accum;
                    atoms.force[i][1] += ft[1] as Accum;
                    atoms.force[i][2] += ft[2] as Accum;
                    dem.torque[i][0] += tau[0];
                    dem.torque[i][1] += tau[1];
                    dem.torque[i][2] += tau[2];
                    new_springs.insert(key, ns);
                }

                // Rolling-resistance torque (cylinder wall is static).
                let mu_r = material_table.rolling_friction_ij[mat_i][wall_mat];
                if mu_r > 0.0 {
                    let sds = material_table.rolling_model == "sds";
                    let k_roll = material_table.rolling_stiffness_ij[mat_i][wall_mat];
                    let gamma_roll = material_table.rolling_damping_ij[mat_i][wall_mat];
                    let key = (1u8, cyl_idx, atoms.tag[i]);
                    let old_rd = old_rolling.get(&key).copied().unwrap_or([0.0; 3]);
                    let (tr, new_rd) = wall_rolling_torque(
                        [nx, ny, nz],
                        dem.omega[i],
                        radius,
                        f_net,
                        mu_r,
                        k_roll,
                        gamma_roll,
                        sds,
                        dt,
                        old_rd,
                    );
                    dem.torque[i][0] += tr[0];
                    dem.torque[i][1] += tr[1];
                    dem.torque[i][2] += tr[2];
                    if sds {
                        new_rolling.insert(key, new_rd);
                    }
                }

                cyl_forces[cyl_idx] += f_net;
            }
        }
        for (idx, &f) in cyl_forces.iter().enumerate() {
            walls.cylinders[idx].force_accumulator += f;
        }

        // ── Sphere walls ────────────────────────────────────────────────────
        let nsph = walls.spheres.len();
        let mut sph_forces = vec![0.0f64; nsph];
        for (sph_idx, sph) in walls.spheres.iter().enumerate() {
            if !walls.sphere_active[sph_idx] {
                continue;
            }
            let wall_mat = sph.material_index;
            for i in 0..nlocal {
                let pos = [
                    atoms.pos[i][0] as f64,
                    atoms.pos[i][1] as f64,
                    atoms.pos[i][2] as f64,
                ];
                let radius = dem.radius[i];
                let mat_i = atoms.atom_type[i] as usize;

                let dx = pos[0] - sph.center[0];
                let dy = pos[1] - sph.center[1];
                let dz = pos[2] - sph.center[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < 1e-30 {
                    continue;
                }

                let inv_dist = 1.0 / dist;
                let (nx, ny, nz, delta) = if sph.inside {
                    let gap = sph.radius - dist;
                    let delta = radius - gap;
                    // Normal points inward (toward center)
                    (-dx * inv_dist, -dy * inv_dist, -dz * inv_dist, delta)
                } else {
                    let gap = dist - sph.radius;
                    let delta = radius - gap;
                    // Normal points outward (away from center)
                    (dx * inv_dist, dy * inv_dist, dz * inv_dist, delta)
                };

                if delta <= 0.0 {
                    continue;
                }
                let delta = delta.min(0.5 * radius);

                let m_r = atoms.mass[i] as f64;
                let v_n = atoms.vel[i][0] as f64 * nx
                    + atoms.vel[i][1] as f64 * ny
                    + atoms.vel[i][2] as f64 * nz;
                let f_net = wall_normal_force(
                    &material_table,
                    mat_i,
                    wall_mat,
                    radius,
                    delta,
                    v_n,
                    m_r,
                    false,
                );

                atoms.force[i][0] += (f_net * nx) as Accum;
                atoms.force[i][1] += (f_net * ny) as Accum;
                atoms.force[i][2] += (f_net * nz) as Accum;

                // Twisting friction torque (sphere wall is static).
                let mu_tw = material_table.twisting_friction_ij[mat_i][wall_mat];
                if mu_tw > 0.0 {
                    let tt = wall_twisting_torque([nx, ny, nz], dem.omega[i], radius, f_net, mu_tw);
                    dem.torque[i][0] += tt[0];
                    dem.torque[i][1] += tt[1];
                    dem.torque[i][2] += tt[2];
                }

                // Tangential (Mindlin) sliding friction (sphere wall is static).
                let mu = material_table.friction_ij[mat_i][wall_mat];
                if mu > 0.0 {
                    let beta = material_table.beta_ij[mat_i][wall_mat];
                    let g_eff = material_table.g_eff_ij[mat_i][wall_mat];
                    let key = (2u8, sph_idx, atoms.tag[i]);
                    let old = old_springs.get(&key).copied().unwrap_or([0.0; 3]);
                    let (ft, tau, ns) = wall_tangential_force(
                        [nx, ny, nz],
                        [
                            atoms.vel[i][0] as f64,
                            atoms.vel[i][1] as f64,
                            atoms.vel[i][2] as f64,
                        ],
                        dem.omega[i],
                        radius,
                        delta,
                        f_net,
                        m_r,
                        g_eff,
                        mu,
                        beta,
                        dt,
                        old,
                    );
                    atoms.force[i][0] += ft[0] as Accum;
                    atoms.force[i][1] += ft[1] as Accum;
                    atoms.force[i][2] += ft[2] as Accum;
                    dem.torque[i][0] += tau[0];
                    dem.torque[i][1] += tau[1];
                    dem.torque[i][2] += tau[2];
                    new_springs.insert(key, ns);
                }

                // Rolling-resistance torque (sphere wall is static).
                let mu_r = material_table.rolling_friction_ij[mat_i][wall_mat];
                if mu_r > 0.0 {
                    let sds = material_table.rolling_model == "sds";
                    let k_roll = material_table.rolling_stiffness_ij[mat_i][wall_mat];
                    let gamma_roll = material_table.rolling_damping_ij[mat_i][wall_mat];
                    let key = (2u8, sph_idx, atoms.tag[i]);
                    let old_rd = old_rolling.get(&key).copied().unwrap_or([0.0; 3]);
                    let (tr, new_rd) = wall_rolling_torque(
                        [nx, ny, nz],
                        dem.omega[i],
                        radius,
                        f_net,
                        mu_r,
                        k_roll,
                        gamma_roll,
                        sds,
                        dt,
                        old_rd,
                    );
                    dem.torque[i][0] += tr[0];
                    dem.torque[i][1] += tr[1];
                    dem.torque[i][2] += tr[2];
                    if sds {
                        new_rolling.insert(key, new_rd);
                    }
                }

                sph_forces[sph_idx] += f_net;
            }
        }
        for (idx, &f) in sph_forces.iter().enumerate() {
            walls.spheres[idx].force_accumulator += f;
        }

        // ── Region walls ────────────────────────────────────────────────────
        let nreg = walls.regions.len();
        let mut reg_forces = vec![0.0f64; nreg];
        for (reg_idx, reg) in walls.regions.iter().enumerate() {
            if !walls.region_active[reg_idx] {
                continue;
            }
            let wall_mat = reg.material_index;
            for i in 0..nlocal {
                let pos = [
                    atoms.pos[i][0] as f64,
                    atoms.pos[i][1] as f64,
                    atoms.pos[i][2] as f64,
                ];
                let radius = dem.radius[i];
                let mat_i = atoms.atom_type[i] as usize;

                let sr = reg.region.closest_point_on_surface(&pos);

                // Compute gap: distance from particle surface to region surface.
                // If inside=true, particles live inside the region, so the wall
                // surface is the region boundary and the gap shrinks as the particle
                // approaches the boundary from inside.
                // sr.distance is positive outside, negative inside.
                let gap = if reg.inside {
                    // Particle inside: gap = |signed_dist| when signed_dist < 0 (inside)
                    // and gap = -signed_dist (positive when inside, wall distance shrinks to 0)
                    -sr.distance
                } else {
                    // Particle outside: gap = signed_dist (positive when outside)
                    sr.distance
                };

                let delta = radius - gap;
                if delta <= 0.0 {
                    continue;
                }
                let delta = delta.min(0.5 * radius);

                // Normal direction: points from wall surface toward particle center
                let (nx, ny, nz) = if reg.inside {
                    // Inside: normal points inward (toward particle, away from surface)
                    (-sr.normal[0], -sr.normal[1], -sr.normal[2])
                } else {
                    // Outside: normal is already outward (toward particle)
                    (sr.normal[0], sr.normal[1], sr.normal[2])
                };

                let m_r = atoms.mass[i] as f64;
                let v_n = atoms.vel[i][0] as f64 * nx
                    + atoms.vel[i][1] as f64 * ny
                    + atoms.vel[i][2] as f64 * nz;
                let f_net = wall_normal_force(
                    &material_table,
                    mat_i,
                    wall_mat,
                    radius,
                    delta,
                    v_n,
                    m_r,
                    false,
                );

                atoms.force[i][0] += (f_net * nx) as Accum;
                atoms.force[i][1] += (f_net * ny) as Accum;
                atoms.force[i][2] += (f_net * nz) as Accum;

                // Twisting friction torque (region wall is static).
                let mu_tw = material_table.twisting_friction_ij[mat_i][wall_mat];
                if mu_tw > 0.0 {
                    let tt = wall_twisting_torque([nx, ny, nz], dem.omega[i], radius, f_net, mu_tw);
                    dem.torque[i][0] += tt[0];
                    dem.torque[i][1] += tt[1];
                    dem.torque[i][2] += tt[2];
                }

                // Tangential (Mindlin) sliding friction (region wall is static).
                let mu = material_table.friction_ij[mat_i][wall_mat];
                if mu > 0.0 {
                    let beta = material_table.beta_ij[mat_i][wall_mat];
                    let g_eff = material_table.g_eff_ij[mat_i][wall_mat];
                    let key = (3u8, reg_idx, atoms.tag[i]);
                    let old = old_springs.get(&key).copied().unwrap_or([0.0; 3]);
                    let (ft, tau, ns) = wall_tangential_force(
                        [nx, ny, nz],
                        [
                            atoms.vel[i][0] as f64,
                            atoms.vel[i][1] as f64,
                            atoms.vel[i][2] as f64,
                        ],
                        dem.omega[i],
                        radius,
                        delta,
                        f_net,
                        m_r,
                        g_eff,
                        mu,
                        beta,
                        dt,
                        old,
                    );
                    atoms.force[i][0] += ft[0] as Accum;
                    atoms.force[i][1] += ft[1] as Accum;
                    atoms.force[i][2] += ft[2] as Accum;
                    dem.torque[i][0] += tau[0];
                    dem.torque[i][1] += tau[1];
                    dem.torque[i][2] += tau[2];
                    new_springs.insert(key, ns);
                }

                // Rolling-resistance torque (region wall is static).
                let mu_r = material_table.rolling_friction_ij[mat_i][wall_mat];
                if mu_r > 0.0 {
                    let sds = material_table.rolling_model == "sds";
                    let k_roll = material_table.rolling_stiffness_ij[mat_i][wall_mat];
                    let gamma_roll = material_table.rolling_damping_ij[mat_i][wall_mat];
                    let key = (3u8, reg_idx, atoms.tag[i]);
                    let old_rd = old_rolling.get(&key).copied().unwrap_or([0.0; 3]);
                    let (tr, new_rd) = wall_rolling_torque(
                        [nx, ny, nz],
                        dem.omega[i],
                        radius,
                        f_net,
                        mu_r,
                        k_roll,
                        gamma_roll,
                        sds,
                        dt,
                        old_rd,
                    );
                    dem.torque[i][0] += tr[0];
                    dem.torque[i][1] += tr[1];
                    dem.torque[i][2] += tr[2];
                    if sds {
                        new_rolling.insert(key, new_rd);
                    }
                }

                reg_forces[reg_idx] += f_net;
            }
        }
        for (idx, &f) in reg_forces.iter().enumerate() {
            walls.regions[idx].force_accumulator += f;
        }

        // Persist the rebuilt tangential- and rolling-spring history for next step.
        walls.tangential_springs = new_springs;
        walls.rolling_springs = new_rolling;
    });
}
