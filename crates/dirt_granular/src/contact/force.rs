use super::*;

/// Standard system: compute the full Hertz-Mindlin contact force in one pass.
pub fn hertz_mindlin_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    particles: ParticlesWith<
        '_,
        (
            Write<DemAtom>,
            Write<ContactHistoryStore>,
            Optional<Read<BondStore>>,
        ),
    >,
    material_table: Res<MaterialTable>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    particles.with(|(mut dem, mut history, bonds)| {
        contact_force_core_views(
            &mut atoms,
            &neighbor,
            &mut dem,
            &mut history,
            bonds.as_deref(),
            &material_table,
            virial.as_deref_mut(),
            ForcePass::All,
        );
    });
}

/// Overlapped Hertz-Mindlin force (roadmap step 4): compute the interior pairs
/// (`j < nlocal`, no ghosts needed) *while the ghost halo is in flight*, then the
/// boundary pairs (`j >= nlocal`) once it lands. The interior force runs as the
/// overlap closure of [`forward_comm_overlap`], so on a multi-rank run its compute
/// hides the MPI latency of the ghost exchange. Bit-identical to
/// [`hertz_mindlin_contact_force`] (the interior/boundary split is exact — see
/// `interior_boundary_split_matches_single_pass`).
pub fn overlapped_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    particles: ParticlesWith<
        '_,
        (
            Write<DemAtom>,
            Write<ContactHistoryStore>,
            Optional<Read<BondStore>>,
        ),
    >,
    material_table: Res<MaterialTable>,
    comm: Res<CommResource>,
    topo: Res<CommTopology>,
    mut buffers: ResMut<CommBuffers>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    let mut pool = std::mem::take(&mut buffers.forward_scratch);
    {
        // Interior pairs need no fresh ghosts — run them during the in-flight halo.
        let mut interior = |a: &mut Atom| {
            contact_force_core_particles(
                a,
                &neighbor,
                &particles,
                &material_table,
                None,
                ForcePass::Interior,
            );
        };
        forward_comm_overlap(
            &mut atoms,
            &registry,
            &topo,
            &**comm,
            &mut pool,
            &mut interior,
        );
    }
    buffers.forward_scratch = pool;
    // Boundary pairs, now that the halo has landed.
    contact_force_core_particles(
        &mut atoms,
        &neighbor,
        &particles,
        &material_table,
        virial.as_deref_mut(),
        ForcePass::Boundary,
    );
}

/// Core force computation, parameterised by [`ForcePass`] so the interior and
/// boundary pairs can be computed in separate passes (with the halo exchange
/// between) for comm/compute overlap. `ForcePass::All` is the single-pass force.
pub fn contact_force_core(
    atoms: &mut Atom,
    neighbor: &Neighbor,
    registry: &AtomDataRegistry,
    material_table: &MaterialTable,
    virial: Option<&mut VirialStress>,
    pass: ForcePass,
) {
    let mut dem = registry.expect_mut::<DemAtom>("contact_force_core test helper");
    let mut history = registry.expect_mut::<ContactHistoryStore>("contact_force_core test helper");
    let bonds = registry.get::<BondStore>();
    contact_force_core_views(
        atoms,
        neighbor,
        &mut dem,
        &mut history,
        bonds.as_deref(),
        material_table,
        virial,
        pass,
    );
}

/// Typed scheduled-contact implementation. The legacy registry helper above is
/// retained solely for direct kernel unit tests; scheduled systems use this
/// query-backed path so required extensions are validated at preparation.
fn contact_force_core_particles(
    atoms: &mut Atom,
    neighbor: &Neighbor,
    particles: &ParticlesWith<
        '_,
        (
            Write<DemAtom>,
            Write<ContactHistoryStore>,
            Optional<Read<BondStore>>,
        ),
    >,
    material_table: &MaterialTable,
    virial: Option<&mut VirialStress>,
    pass: ForcePass,
) {
    particles.with(|(mut dem, mut history, bonds)| {
        contact_force_core_views(
            atoms,
            neighbor,
            &mut dem,
            &mut history,
            bonds.as_deref(),
            material_table,
            virial,
            pass,
        );
    });
}

/// Typed core used by scheduled contact systems after `ParticlesWith` has
/// validated and borrowed the exact extension columns they require.
pub fn contact_force_core_views(
    atoms: &mut Atom,
    neighbor: &Neighbor,
    dem: &mut DemAtom,
    history: &mut ContactHistoryStore,
    bond_store: Option<&BondStore>,
    material_table: &MaterialTable,
    mut virial: Option<&mut VirialStress>,
    pass: ForcePass,
) {
    let newton = neighbor.newton;
    let dt = atoms.dt;

    let natoms = atoms.len();
    if history.contacts.len() < natoms {
        history.contacts.resize_with(natoms, Vec::new);
    }

    let nlocal = atoms.nlocal as usize;
    let mut overlap_warnings = 0usize;

    // Reset all active flags before pair loop (skipped on the Boundary pass, which
    // continues the Interior pass's history instead of clearing it).
    if pass != ForcePass::Boundary {
        for i in 0..nlocal {
            for entry in &mut history.contacts[i] {
                entry.2 = false;
            }
        }
    }

    for (i, j) in neighbor.pairs(nlocal) {
        // Interior/boundary split (step 4): boundary pairs touch a ghost (j >= nlocal).
        let is_boundary_pair = j >= nlocal;
        match pass {
            ForcePass::Interior if is_boundary_pair => continue,
            ForcePass::Boundary if !is_boundary_pair => continue,
            _ => {}
        }
        if let Some(ref bonds) = bond_store {
            if bonds.are_excluded(i, j, &atoms.tag) {
                continue;
            }
        }

        // Skip same-body pairs (sub-spheres of the same rigid body don't interact)
        if dirt_atom::same_body(&dem, i, j) {
            continue;
        }

        let r1 = dem.radius[i];
        let r2 = dem.radius[j];

        let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
        let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
        let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        let sum_r = r1 + r2;

        let mat_i = atoms.atom_type[i] as usize;
        let mat_j = atoms.atom_type[j] as usize;
        let surface_energy = material_table.surface_energy_ij[mat_i][mat_j];

        let use_dmt = material_table.adhesion_model == "dmt";

        // JKR: compute pull-off distance for extended interaction range
        // DMT: no extended range (particles separate at delta = 0)
        // Effective radius: R* = R1 R2 / (R1 + R2)
        let r_eff = (r1 * r2) / sum_r;
        // Effective Young's modulus: 1/E* = (1-ν1²)/E1 + (1-ν2²)/E2
        let e_eff = material_table.e_eff_ij[mat_i][mat_j];
        let liquid_bridge_active = material_table.liquid_bridge_model == "willett2000";
        let liquid_volume = material_table.liquid_bridge_volume_ij[mat_i][mat_j];
        let liquid_surface_tension = material_table.liquid_surface_tension_ij[mat_i][mat_j];
        let liquid_contact_angle = material_table.liquid_contact_angle_ij[mat_i][mat_j];
        let liquid_rupture_distance = material_table.liquid_rupture_distance_ij[mat_i][mat_j];
        let liquid_range =
            if liquid_bridge_active && liquid_volume > 0.0 && liquid_surface_tension > 0.0 {
                if liquid_rupture_distance > 0.0 {
                    liquid_rupture_distance
                } else {
                    (1.0 + 0.5 * liquid_contact_angle) * liquid_volume.cbrt()
                }
            } else {
                0.0
            };
        // JKR pull-off distance: particles interact beyond geometric contact
        let delta_pulloff = if surface_energy > 0.0 && !use_dmt {
            let gamma = surface_energy;
            (std::f64::consts::PI * std::f64::consts::PI * gamma * gamma * r_eff
                / (4.0 * e_eff * e_eff))
                .cbrt()
        } else {
            0.0
        };

        // Check contact: geometric touch or within JKR adhesion range
        let interaction_r = sum_r + delta_pulloff.max(liquid_range);
        if dist_sq >= interaction_r * interaction_r {
            continue;
        }

        let distance = dist_sq.sqrt();

        if distance == 0.0 {
            #[cfg(debug_assertions)]
            eprintln!(
                "WARNING: zero separation between tags {} {}",
                atoms.tag[i], atoms.tag[j]
            );
            continue;
        }

        // delta > 0 means geometric overlap, delta < 0 means gap
        // Cap at half the smaller radius to keep the Hertz model numerically valid.
        let r_min = r1.min(r2);
        let delta = (sum_r - distance).min(0.5 * r_min);

        if delta > 0.0 && distance / sum_r < LARGE_OVERLAP_WARN_THRESHOLD {
            overlap_warnings += 1;
            #[cfg(debug_assertions)]
            eprintln!(
                "WARNING: large overlap tags {} {} ratio {:.3}",
                atoms.tag[i],
                atoms.tag[j],
                distance / sum_r
            );
            if overlap_warnings > MAX_OVERLAP_WARNINGS {
                panic!(
                    "Over {} excessive overlaps this step — aborting. \
                     Check timestep or initial configuration.",
                    MAX_OVERLAP_WARNINGS
                );
            }
            // Cap overlap at half the smaller radius to keep Hertz model valid,
            // but still compute the repulsive force (skipping would remove all
            // repulsion and cause runaway penetration).
        }

        let separation = (distance - sum_r).max(0.0);

        // For non-JKR/non-liquid, skip if no geometric overlap
        if delta <= 0.0 && surface_energy <= 0.0 && liquid_range <= 0.0 {
            continue;
        }

        // ── Shared quantities (computed once) ────────────────────────────
        let inv_dist = 1.0 / distance;
        let nx = dx * inv_dist;
        let ny = dy * inv_dist;
        let nz = dz * inv_dist;

        // Effective shear modulus: 1/G* = (2-ν1)/G1 + (2-ν2)/G2
        let g_eff = material_table.g_eff_ij[mat_i][mat_j];

        // Reduced mass: m_r = 1 / (1/m1 + 1/m2)
        // For clump sub-spheres inv_mass is 0 (body-integrated); use real mass.
        let inv_m_i = if atoms.inv_mass[i] as f64 > 0.0 {
            atoms.inv_mass[i] as f64
        } else {
            1.0 / atoms.mass[i] as f64
        };
        let inv_m_j = if atoms.inv_mass[j] as f64 > 0.0 {
            atoms.inv_mass[j] as f64
        } else {
            1.0 / atoms.mass[j] as f64
        };
        let m_r = 1.0 / (inv_m_i + inv_m_j);

        let beta = material_table.beta_ij[mat_i][mat_j];
        let mu = material_table.friction_ij[mat_i][mat_j];
        let mu_r = material_table.rolling_friction_ij[mat_i][mat_j];
        let mu_tw = material_table.twisting_friction_ij[mat_i][mat_j];
        let cohesion_energy = material_table.cohesion_energy_ij[mat_i][mat_j];
        let use_mdr = material_table.contact_model == "mdr";

        // JKR adhesion-only regime: gap exists but within pull-off distance
        // DMT has no adhesion-only regime (no force beyond contact)
        let jkr_adhesion_only = surface_energy > 0.0 && !use_dmt && delta <= 0.0;
        let f_liquid_bridge = if liquid_range > 0.0 {
            willett2000_liquid_bridge_force(
                separation,
                r_eff,
                liquid_volume,
                liquid_surface_tension,
                liquid_contact_angle,
                liquid_rupture_distance,
            )
        } else {
            0.0
        };
        let bridge_only = delta <= 0.0 && !jkr_adhesion_only && f_liquid_bridge > 0.0;

        // Hertz stiffness parameters (only meaningful when δ > 0)
        // S_n = 2 E* √(R* δ)  — normal stiffness parameter (used in damping)
        // k_n = 4/3 E* √(R* δ) — normal spring constant
        // k_t = 8 G* √(R* δ)  — tangential spring constant (Mindlin)
        let (s_n, k_n, k_t, contact_radius) = if delta > 0.0 {
            let sdr = (delta * r_eff).sqrt();
            let sn = 2.0 * e_eff * sdr;
            let kn = 4.0 / 3.0 * e_eff * sdr;
            let kt = 8.0 * g_eff * sdr;
            (sn, kn, kt, sdr)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };

        // Full relative velocity (including angular contributions)
        let omega_ix = dem.omega[i][0];
        let omega_iy = dem.omega[i][1];
        let omega_iz = dem.omega[i][2];
        let omega_jx = dem.omega[j][0];
        let omega_jy = dem.omega[j][1];
        let omega_jz = dem.omega[j][2];

        // v_contact_i = vel_i + omega_i × (r1 * n)
        let r1n_x = r1 * nx;
        let r1n_y = r1 * ny;
        let r1n_z = r1 * nz;
        let vc_ix = atoms.vel[i][0] as f64 + (omega_iy * r1n_z - omega_iz * r1n_y);
        let vc_iy = atoms.vel[i][1] as f64 + (omega_iz * r1n_x - omega_ix * r1n_z);
        let vc_iz = atoms.vel[i][2] as f64 + (omega_ix * r1n_y - omega_iy * r1n_x);

        // v_contact_j = vel_j + omega_j × (-r2 * n)
        let r2n_x = r2 * nx;
        let r2n_y = r2 * ny;
        let r2n_z = r2 * nz;
        let vc_jx = atoms.vel[j][0] as f64 + (-omega_jy * r2n_z + omega_jz * r2n_y);
        let vc_jy = atoms.vel[j][1] as f64 + (-omega_jz * r2n_x + omega_jx * r2n_z);
        let vc_jz = atoms.vel[j][2] as f64 + (-omega_jx * r2n_y + omega_jy * r2n_x);

        let vr_x = vc_jx - vc_ix;
        let vr_y = vc_jy - vc_iy;
        let vr_z = vc_jz - vc_iz;

        let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

        let tag_i = atoms.tag[i];
        let tag_j = atoms.tag[j];
        let sign: f64 = if tag_i < tag_j { 1.0 } else { -1.0 };

        // Look up existing spring/history (single search, reused for write-back).
        let entry_idx = history.contacts[i].iter().position(|(t, _, _)| *t == tag_j);
        let mut stored = match entry_idx {
            Some(idx) => history.contacts[i][idx].1,
            None => zero_contact_history(),
        };

        // ── Normal force ─────────────────────────────────────────────────
        // F_n > 0 → repulsive (along contact normal from i to j)
        // F_n < 0 → attractive (adhesion/cohesion pulls particles together)
        let (mut f_n_mag, k_t) = if use_mdr {
            let (f, _k_mdr, kt_mdr) = mdr_normal_force(
                delta.max(0.0),
                r1,
                r2,
                e_eff,
                0.5 * (material_table.poisson_ratio[mat_i] + material_table.poisson_ratio[mat_j]),
                material_table.mdr_yield_stress_ij[mat_i][mat_j],
                surface_energy,
                material_table.mdr_damping_ij[mat_i][mat_j],
                m_r,
                v_n,
                &mut stored,
            );
            (f, kt_mdr)
        } else if surface_energy > 0.0 && use_dmt {
            // DMT: Hertz contact + constant adhesive force F_dmt = 2π γ R*
            let f_dmt = 2.0 * std::f64::consts::PI * surface_energy * r_eff;
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            (k_n * delta - f_diss_n - f_dmt, k_t)
        } else if surface_energy > 0.0 {
            // JKR: adhesion force F_adh = 3/2 π γ R* (simplified explicit model)
            let f_adhesion = 1.5 * std::f64::consts::PI * surface_energy * r_eff;
            if jkr_adhesion_only {
                // Gap regime (δ ≤ 0): pure adhesion, no Hertz contact or damping
                (-f_adhesion, k_t)
            } else {
                // Contact regime (δ > 0): Hertz repulsion + damping − adhesion
                let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
                (k_n * delta - f_diss_n - f_adhesion, k_t)
            }
        } else if cohesion_energy > 0.0 {
            // SJKR: cohesion proportional to contact area A = π δ R*
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * r_eff;
            (k_n * delta - f_diss_n - f_cohesion, k_t) // can go negative (attractive)
        } else {
            // Standard Hertz repulsion + viscoelastic damping. With
            // `limit_damping` (default) the total is clamped to ≥ 0 so damping
            // can never pull particles together; with it disabled the damping may
            // go net-attractive near separation, matching LAMMPS's default
            // `pair granular` (no tensile cutoff) — required for exact cross-code
            // COR at low restitution (see bench_hertz_rebound).
            let f_diss_n = 2.0 * beta * SQRT_5_6 * (s_n * m_r).sqrt() * v_n;
            let f_total = k_n * delta - f_diss_n;
            if material_table.limit_damping {
                (f_total.max(0.0), k_t)
            } else {
                (f_total, k_t)
            }
        };
        f_n_mag -= f_liquid_bridge;

        let fn_x = f_n_mag * nx;
        let fn_y = f_n_mag * ny;
        let fn_z = f_n_mag * nz;

        atoms.force[i][0] -= fn_x as soil_core::Accum;
        atoms.force[i][1] -= fn_y as soil_core::Accum;
        atoms.force[i][2] -= fn_z as soil_core::Accum;
        if newton {
            atoms.force[j][0] += fn_x as soil_core::Accum;
            atoms.force[j][1] += fn_y as soil_core::Accum;
            atoms.force[j][2] += fn_z as soil_core::Accum;
        }

        // ── Tangential force (skip in JKR adhesion-only regime) ──────────
        // No tangential friction when particles are not in geometric contact
        if jkr_adhesion_only || bridge_only {
            // No tangential, rolling, or spring history in gap-only attraction.
            // Virial contribution from normal only
            if let Some(ref mut v) = virial {
                if v.active {
                    let vs = if newton { 1.0 } else { 0.5 };
                    v.add_pair(dx, dy, dz, -fn_x * vs, -fn_y * vs, -fn_z * vs);
                }
            }
            continue;
        }

        let vt_x = vr_x - v_n * nx;
        let vt_y = vr_y - v_n * ny;
        let vt_z = vr_z - v_n * nz;

        // Tangential spring displacement (history model). For the history-free
        // `linear_nohistory` model the spring is identically zero, so the force
        // collapses to the velocity-Coulomb law
        //   F_t = -min(μ |F_n|, γ_t |v_t|) t̂ ,   t̂ = v_t / |v_t|
        // (LAMMPS pair_granular `tangential linear_nohistory`, and the classic
        // `pair gran/hooke`) — the force depends only on the instantaneous
        // relative tangential velocity, with NO accumulated displacement.
        let tangential_model = material_table.tangential_model.as_str();
        let nohistory = tangential_model == "linear_nohistory";
        let mindlin_rescale = is_mindlin_rescale(tangential_model);
        let mindlin_force_history = is_mindlin_force_history(tangential_model);
        let f_t_max = mu * f_n_mag.abs();
        let (sx, sy, sz) = if nohistory {
            (0.0, 0.0, 0.0)
        } else if mindlin_force_history {
            // LAMMPS `mindlin_rescale/force` stores the elastic tangential force
            // itself as history. On normal unloading it scales that force by the
            // contact-radius ratio a_n/a_{n-1} before adding the new increment.
            let mut fx = sign * stored[0];
            let mut fy = sign * stored[1];
            let mut fz = sign * stored[2];
            let prev_a = stored[7];
            if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                let scale = contact_radius / prev_a;
                fx *= scale;
                fy *= scale;
                fz *= scale;
            }
            let f_dot_n = fx * nx + fy * ny + fz * nz;
            fx -= f_dot_n * nx;
            fy -= f_dot_n * ny;
            fz -= f_dot_n * nz;
            fx += k_t * vt_x * dt;
            fy += k_t * vt_y * dt;
            fz += k_t * vt_z * dt;
            (fx, fy, fz)
        } else {
            // Convert stored spring from canonical form to local (i,j) frame
            let mut sx = sign * stored[0];
            let mut sy = sign * stored[1];
            let mut sz = sign * stored[2];
            let prev_a = stored[7];
            if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                let scale = contact_radius / prev_a;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }
            // Rotate spring into current tangent plane (remove normal component)
            let s_dot_n = sx * nx + sy * ny + sz * nz;
            sx -= s_dot_n * nx;
            sy -= s_dot_n * ny;
            sz -= s_dot_n * nz;
            // Integrate tangential velocity into spring displacement
            sx += vt_x * dt;
            sy += vt_y * dt;
            sz += vt_z * dt;

            // Coulomb cap on spring: |k_t s| ≤ μ |F_n|
            let s_mag = (sx * sx + sy * sy + sz * sz).sqrt();
            let f_t_spring_mag = k_t * s_mag;
            if f_t_spring_mag > f_t_max && f_t_spring_mag > TANGENTIAL_EPSILON {
                let scale = f_t_max / f_t_spring_mag;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }
            (sx, sy, sz)
        };

        // Tangential damping coefficient: γ_t = 2 β √(5/6) √(k_t m_r)
        let gamma_t = 2.0 * SQRT_5_6 * beta * (k_t * m_r).sqrt();
        let mut ft_x = (if mindlin_force_history { sx } else { k_t * sx }) + gamma_t * vt_x;
        let mut ft_y = (if mindlin_force_history { sy } else { k_t * sy }) + gamma_t * vt_y;
        let mut ft_z = (if mindlin_force_history { sz } else { k_t * sz }) + gamma_t * vt_z;

        // Coulomb cap on total tangential force
        let f_t_mag = (ft_x * ft_x + ft_y * ft_y + ft_z * ft_z).sqrt();
        if f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
            let scale = f_t_max / f_t_mag;
            ft_x *= scale;
            ft_y *= scale;
            ft_z *= scale;
        }

        let (sx, sy, sz) =
            if mindlin_force_history && f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
                (
                    ft_x - gamma_t * vt_x,
                    ft_y - gamma_t * vt_y,
                    ft_z - gamma_t * vt_z,
                )
            } else {
                (sx, sy, sz)
            };

        // Torques: τ_i = (r1 * n) × f_t, τ_j = (-r2 * n) × (-f_t) = (r2 * n) × f_t
        let ti_x = r1n_y * ft_z - r1n_z * ft_y;
        let ti_y = r1n_z * ft_x - r1n_x * ft_z;
        let ti_z = r1n_x * ft_y - r1n_y * ft_x;
        let tj_x = r2n_y * ft_z - r2n_z * ft_y;
        let tj_y = r2n_z * ft_x - r2n_x * ft_z;
        let tj_z = r2n_x * ft_y - r2n_y * ft_x;

        atoms.force[i][0] += ft_x as soil_core::Accum;
        atoms.force[i][1] += ft_y as soil_core::Accum;
        atoms.force[i][2] += ft_z as soil_core::Accum;
        if newton {
            atoms.force[j][0] -= ft_x as soil_core::Accum;
            atoms.force[j][1] -= ft_y as soil_core::Accum;
            atoms.force[j][2] -= ft_z as soil_core::Accum;
        }
        dem.torque[i][0] += ti_x;
        dem.torque[i][1] += ti_y;
        dem.torque[i][2] += ti_z;
        if newton {
            dem.torque[j][0] += tj_x;
            dem.torque[j][1] += tj_y;
            dem.torque[j][2] += tj_z;
        }

        // ── Rolling resistance torque ───────────────────────────────────
        // Relative angular velocity (rolling component)
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
                // SDS rolling: spring-dashpot-slider model
                let k_roll = material_table.rolling_stiffness_ij[mat_i][mat_j];
                let gamma_roll = material_table.rolling_damping_ij[mat_i][mat_j];

                // Update rolling displacement: remove normal component, integrate
                let rd_dot_n = roll_disp_x * nx + roll_disp_y * ny + roll_disp_z * nz;
                roll_disp_x -= rd_dot_n * nx;
                roll_disp_y -= rd_dot_n * ny;
                roll_disp_z -= rd_dot_n * nz;
                roll_disp_x += roll_x * dt;
                roll_disp_y += roll_y * dt;
                roll_disp_z += roll_z * dt;

                // Spring + dashpot torque
                let mut tr_x = -k_roll * roll_disp_x - gamma_roll * roll_x;
                let mut tr_y = -k_roll * roll_disp_y - gamma_roll * roll_y;
                let mut tr_z = -k_roll * roll_disp_z - gamma_roll * roll_z;
                let tr_mag = (tr_x * tr_x + tr_y * tr_y + tr_z * tr_z).sqrt();
                let tau_max = mu_r * f_n_mag.abs() * r_eff;

                if tr_mag > tau_max && tr_mag > TANGENTIAL_EPSILON {
                    // Cap and rescale spring displacement
                    let scale = tau_max / tr_mag;
                    tr_x *= scale;
                    tr_y *= scale;
                    tr_z *= scale;
                    // Rescale spring: δ = (τ + γ·ω) / (-k)
                    if k_roll > TANGENTIAL_EPSILON {
                        roll_disp_x = (tr_x + gamma_roll * roll_x) / (-k_roll);
                        roll_disp_y = (tr_y + gamma_roll * roll_y) / (-k_roll);
                        roll_disp_z = (tr_z + gamma_roll * roll_z) / (-k_roll);
                    }
                }

                dem.torque[i][0] += tr_x;
                dem.torque[i][1] += tr_y;
                dem.torque[i][2] += tr_z;
                if newton {
                    dem.torque[j][0] -= tr_x;
                    dem.torque[j][1] -= tr_y;
                    dem.torque[j][2] -= tr_z;
                }
            } else if roll_mag > 1e-30 {
                // Constant torque model (existing behavior)
                let tau_mag = mu_r * f_n_mag.abs() * r_eff;
                let inv_roll = tau_mag / roll_mag;
                let tr_x = -inv_roll * roll_x;
                let tr_y = -inv_roll * roll_y;
                let tr_z = -inv_roll * roll_z;
                dem.torque[i][0] += tr_x;
                dem.torque[i][1] += tr_y;
                dem.torque[i][2] += tr_z;
                if newton {
                    dem.torque[j][0] -= tr_x;
                    dem.torque[j][1] -= tr_y;
                    dem.torque[j][2] -= tr_z;
                }
            }
        }

        // ── Twisting friction torque ─────────────────────────────────────
        // Three selectable models (material_table.twisting_model):
        //   "constant" / "sds"  — user-supplied coefficients (gated on μ_tw > 0);
        //   "marshall"           — coefficients DERIVED from the active tangential
        //                          (Mindlin) model, no separate twist inputs.
        if material_table.twisting_model == "marshall" {
            // Marshall (2009) twisting, per LAMMPS pair_granular `twisting marshall`
            // (doc/src/pair_granular.rst §twisting, Marshall2009 eqs 32-33). The
            // twisting stiffness/damping/friction are expressed in terms of the
            // tangential (sliding) coefficients and the Hertz contact radius
            // a = √(R* δ):
            //     k_twist = ½ k_t a²,  γ_twist = ½ γ_t a²,  μ_twist = (2/3) a μ_t
            // with k_t, γ_t the tangential spring/damping computed above and μ_t
            // the tangential friction coefficient. Below the cap the couple is
            // the spring–dashpot τ = −k_twist ξ − γ_twist Ω; it is then truncated
            // to |τ| ≤ μ_twist F_n and the angular displacement rescaled to the
            // critical value (identical bookkeeping to the SDS slider).
            if delta > 0.0 {
                let twist_vel = or_dot_n;
                let a_sq = delta * r_eff; // a² = (√(R* δ))²
                let a = a_sq.sqrt(); // Hertz contact radius
                let k_twist = 0.5 * k_t * a_sq;
                let gamma_twist = 0.5 * gamma_t * a_sq;
                let mu_twist = (2.0 / 3.0) * a * mu; // μ = tangential friction coeff

                twist_disp += twist_vel * dt;

                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_twist * f_n_mag.abs();
                if tau_twist.abs() > tau_max {
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                let tt_x = tau_twist * nx;
                let tt_y = tau_twist * ny;
                let tt_z = tau_twist * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            }
        } else if mu_tw > 0.0 {
            let twist_vel = or_dot_n; // twisting component of relative angular velocity
            let sds_twisting = material_table.twisting_model == "sds";
            if sds_twisting {
                // SDS twisting: spring-dashpot-slider model
                let k_twist = material_table.twisting_stiffness_ij[mat_i][mat_j];
                let gamma_twist = material_table.twisting_damping_ij[mat_i][mat_j];

                // Update twisting displacement
                twist_disp += twist_vel * dt;

                // Spring + dashpot torque (scalar, along contact normal)
                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_tw * f_n_mag.abs() * r_eff;

                if tau_twist.abs() > tau_max {
                    // Cap and rescale spring
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                let tt_x = tau_twist * nx;
                let tt_y = tau_twist * ny;
                let tt_z = tau_twist * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            } else if twist_vel.abs() > 1e-30 {
                // Constant torque model (existing behavior)
                let tau = mu_tw * f_n_mag.abs() * r_eff;
                let sign_tw = if twist_vel > 0.0 { -1.0 } else { 1.0 };
                let tt_x = sign_tw * tau * nx;
                let tt_y = sign_tw * tau * ny;
                let tt_z = sign_tw * tau * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            }
        }

        // Virial: force on i from j = (-fn + ft)
        // When newton=false, each pair is visited twice so halve virial contribution
        if let Some(ref mut v) = virial {
            if v.active {
                let vs = if newton { 1.0 } else { 0.5 };
                let vfx = (-fn_x + ft_x) * vs;
                let vfy = (-fn_y + ft_y) * vs;
                let vfz = (-fn_z + ft_z) * vs;
                v.add_pair(dx, dy, dz, vfx, vfy, vfz);
            }
        }

        // Store updated spring back (canonical form) and mark active
        let mut new_spring = stored;
        new_spring[0] = sign * sx;
        new_spring[1] = sign * sy;
        new_spring[2] = sign * sz;
        new_spring[3] = sign * roll_disp_x;
        new_spring[4] = sign * roll_disp_y;
        new_spring[5] = sign * roll_disp_z;
        new_spring[6] = sign * twist_disp;
        new_spring[7] = contact_radius;
        match entry_idx {
            Some(idx) => {
                history.contacts[i][idx].1 = new_spring;
                history.contacts[i][idx].2 = true;
            }
            None => history.contacts[i].push((tag_j, new_spring, true)),
        }
    }

    // Prune stale contacts (skipped on the Interior pass; the Boundary pass prunes
    // once after both passes have marked their active contacts).
    if pass != ForcePass::Interior {
        for i in 0..nlocal {
            history.contacts[i].retain(|(_, _, active)| *active);
        }
    }

    // Debug: check total force + torque on all atoms (local + ghost).
    // In a correct Newton's 3rd law implementation, the sum of all forces
    // from pair interactions must be zero (each pair contributes +F to one atom
    // and -F to the other). A nonzero sum means a pair was counted asymmetrically.
    // Skip this check when newton=false (forces only written to i).
    #[cfg(debug_assertions)]
    if newton {
        let total = atoms.len();
        let mut sum_fx = 0.0;
        let mut sum_fy = 0.0;
        let mut sum_fz = 0.0;
        for i in 0..total {
            sum_fx += atoms.force[i][0] as f64;
            sum_fy += atoms.force[i][1] as f64;
            sum_fz += atoms.force[i][2] as f64;
        }
        let sum_f = (sum_fx * sum_fx + sum_fy * sum_fy + sum_fz * sum_fz).sqrt();
        if sum_f > 1e-6 {
            eprintln!(
                "WARNING: nonzero net force after contact: |F|={:.6e} ({:.6e},{:.6e},{:.6e})",
                sum_f, sum_fx, sum_fy, sum_fz
            );
        }
    }
}

/// Hooke (linear spring) contact force — alternative to Hertz-Mindlin.
///
/// Normal: `f_n = kn * delta`, tangential uses `kt` directly.
/// Damping: `gamma = 2 * beta * sqrt(kn_ij * m_r)`.
/// All other features (friction, rolling, twisting, cohesion, JKR) reused.
pub fn hooke_contact_force(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    particles: ParticlesWith<
        '_,
        (
            Write<DemAtom>,
            Write<ContactHistoryStore>,
            Optional<Read<BondStore>>,
        ),
    >,
    material_table: Res<MaterialTable>,
    mut virial: Option<ResMut<VirialStress>>,
) {
    particles.with(|(mut dem, mut history, bond_store)| {
        let newton = neighbor.newton;
        let dt = atoms.dt;

        while history.contacts.len() < atoms.len() {
            history.contacts.push(Vec::new());
        }

        let nlocal = atoms.nlocal as usize;
        let mut overlap_warnings = 0usize;

        for i in 0..nlocal {
            for entry in &mut history.contacts[i] {
                entry.2 = false;
            }
        }

        for (i, j) in neighbor.pairs(nlocal) {
            if let Some(ref bonds) = bond_store {
                if bonds.are_excluded(i, j, &atoms.tag) {
                    continue;
                }
            }

            // Skip same-body pairs (sub-spheres of the same rigid body don't interact)
            if dirt_atom::same_body(&dem, i, j) {
                continue;
            }

            let r1 = dem.radius[i];
            let r2 = dem.radius[j];

            let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
            let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
            let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
            let dist_sq = dx * dx + dy * dy + dz * dz;
            let sum_r = r1 + r2;

            if dist_sq >= sum_r * sum_r {
                continue;
            }

            let distance = dist_sq.sqrt();
            if distance == 0.0 {
                continue;
            }

            let r_min = r1.min(r2);
            let delta = (sum_r - distance).min(0.5 * r_min);
            if delta <= 0.0 {
                continue;
            }

            if distance / sum_r < LARGE_OVERLAP_WARN_THRESHOLD {
                overlap_warnings += 1;
                if overlap_warnings > MAX_OVERLAP_WARNINGS {
                    panic!(
                        "Over {} excessive overlaps this step — aborting.",
                        MAX_OVERLAP_WARNINGS
                    );
                }
                // Still compute force (don't skip) — removing repulsion causes runaway.
            }

            let inv_dist = 1.0 / distance;
            let nx = dx * inv_dist;
            let ny = dy * inv_dist;
            let nz = dz * inv_dist;

            let mat_i = atoms.atom_type[i] as usize;
            let mat_j = atoms.atom_type[j] as usize;
            let r_eff = (r1 * r2) / sum_r;
            // For clump sub-spheres inv_mass is 0 (body-integrated); use real mass.
            let inv_m_i = if atoms.inv_mass[i] as f64 > 0.0 {
                atoms.inv_mass[i] as f64
            } else {
                1.0 / atoms.mass[i] as f64
            };
            let inv_m_j = if atoms.inv_mass[j] as f64 > 0.0 {
                atoms.inv_mass[j] as f64
            } else {
                1.0 / atoms.mass[j] as f64
            };
            let m_r = 1.0 / (inv_m_i + inv_m_j);
            let beta = material_table.beta_ij[mat_i][mat_j];
            let mu = material_table.friction_ij[mat_i][mat_j];
            let mu_r = material_table.rolling_friction_ij[mat_i][mat_j];
            let mu_tw = material_table.twisting_friction_ij[mat_i][mat_j];
            let cohesion_energy = material_table.cohesion_energy_ij[mat_i][mat_j];

            let kn = material_table.kn_ij[mat_i][mat_j];
            let kt = material_table.kt_ij[mat_i][mat_j];
            let contact_radius = (r_eff * delta).sqrt();

            // Hooke normal: f_n = kn * delta
            // Damping: gamma_n = 2 * beta * sqrt(kn * m_r)
            let gamma_n = 2.0 * beta * (kn * m_r).sqrt();

            // Relative velocity
            let omega_ix = dem.omega[i][0];
            let omega_iy = dem.omega[i][1];
            let omega_iz = dem.omega[i][2];
            let omega_jx = dem.omega[j][0];
            let omega_jy = dem.omega[j][1];
            let omega_jz = dem.omega[j][2];

            let r1n_x = r1 * nx;
            let r1n_y = r1 * ny;
            let r1n_z = r1 * nz;
            let vc_ix = atoms.vel[i][0] as f64 + (omega_iy * r1n_z - omega_iz * r1n_y);
            let vc_iy = atoms.vel[i][1] as f64 + (omega_iz * r1n_x - omega_ix * r1n_z);
            let vc_iz = atoms.vel[i][2] as f64 + (omega_ix * r1n_y - omega_iy * r1n_x);

            let r2n_x = r2 * nx;
            let r2n_y = r2 * ny;
            let r2n_z = r2 * nz;
            let vc_jx = atoms.vel[j][0] as f64 + (-omega_jy * r2n_z + omega_jz * r2n_y);
            let vc_jy = atoms.vel[j][1] as f64 + (-omega_jz * r2n_x + omega_jx * r2n_z);
            let vc_jz = atoms.vel[j][2] as f64 + (-omega_jx * r2n_y + omega_jy * r2n_x);

            let vr_x = vc_jx - vc_ix;
            let vr_y = vc_jy - vc_iy;
            let vr_z = vc_jz - vc_iz;
            let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

            // Normal force
            let f_n_mag = if cohesion_energy > 0.0 {
                let f_cohesion = cohesion_energy * std::f64::consts::PI * delta * r_eff;
                kn * delta - gamma_n * v_n - f_cohesion
            } else {
                // See the Hertz path: `limit_damping` (default) clamps to repulsive-
                // only; disabling it matches LAMMPS's default (no tensile cutoff).
                let f_total = kn * delta - gamma_n * v_n;
                if material_table.limit_damping {
                    f_total.max(0.0)
                } else {
                    f_total
                }
            };

            let fn_x = f_n_mag * nx;
            let fn_y = f_n_mag * ny;
            let fn_z = f_n_mag * nz;

            atoms.force[i][0] -= fn_x as soil_core::Accum;
            atoms.force[i][1] -= fn_y as soil_core::Accum;
            atoms.force[i][2] -= fn_z as soil_core::Accum;
            if newton {
                atoms.force[j][0] += fn_x as soil_core::Accum;
                atoms.force[j][1] += fn_y as soil_core::Accum;
                atoms.force[j][2] += fn_z as soil_core::Accum;
            }

            // Tangential force
            let vt_x = vr_x - v_n * nx;
            let vt_y = vr_y - v_n * ny;
            let vt_z = vr_z - v_n * nz;

            let tag_i = atoms.tag[i];
            let tag_j = atoms.tag[j];
            let sign: f64 = if tag_i < tag_j { 1.0 } else { -1.0 };

            let entry_idx = history.contacts[i].iter().position(|(t, _, _)| *t == tag_j);
            let stored = match entry_idx {
                Some(idx) => history.contacts[i][idx].1,
                None => zero_contact_history(),
            };

            // History-free `linear_nohistory` tangential model → zero spring (see the
            // Hertz path above); the force reduces to velocity-Coulomb with no
            // accumulated displacement. "history" keeps the incremental Hooke spring.
            let tangential_model = material_table.tangential_model.as_str();
            let nohistory = tangential_model == "linear_nohistory";
            let mindlin_rescale = is_mindlin_rescale(tangential_model);
            let mindlin_force_history = is_mindlin_force_history(tangential_model);
            let f_t_max = mu * f_n_mag.abs();
            let (sx, sy, sz) = if nohistory {
                (0.0, 0.0, 0.0)
            } else if mindlin_force_history {
                let mut fx = sign * stored[0];
                let mut fy = sign * stored[1];
                let mut fz = sign * stored[2];
                let prev_a = stored[7];
                if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                    let scale = contact_radius / prev_a;
                    fx *= scale;
                    fy *= scale;
                    fz *= scale;
                }
                let f_dot_n = fx * nx + fy * ny + fz * nz;
                fx -= f_dot_n * nx;
                fy -= f_dot_n * ny;
                fz -= f_dot_n * nz;
                fx += kt * vt_x * dt;
                fy += kt * vt_y * dt;
                fz += kt * vt_z * dt;
                (fx, fy, fz)
            } else {
                let mut sx = sign * stored[0];
                let mut sy = sign * stored[1];
                let mut sz = sign * stored[2];
                let prev_a = stored[7];
                if mindlin_rescale && prev_a > TANGENTIAL_EPSILON && contact_radius < prev_a {
                    let scale = contact_radius / prev_a;
                    sx *= scale;
                    sy *= scale;
                    sz *= scale;
                }
                let s_dot_n = sx * nx + sy * ny + sz * nz;
                sx -= s_dot_n * nx;
                sy -= s_dot_n * ny;
                sz -= s_dot_n * nz;
                sx += vt_x * dt;
                sy += vt_y * dt;
                sz += vt_z * dt;

                let s_mag = (sx * sx + sy * sy + sz * sz).sqrt();
                let f_t_spring_mag = kt * s_mag;
                if f_t_spring_mag > f_t_max && f_t_spring_mag > TANGENTIAL_EPSILON {
                    let scale = f_t_max / f_t_spring_mag;
                    sx *= scale;
                    sy *= scale;
                    sz *= scale;
                }
                (sx, sy, sz)
            };

            let gamma_t = 2.0 * SQRT_5_6 * beta * (kt * m_r).sqrt();
            let mut ft_x = (if mindlin_force_history { sx } else { kt * sx }) + gamma_t * vt_x;
            let mut ft_y = (if mindlin_force_history { sy } else { kt * sy }) + gamma_t * vt_y;
            let mut ft_z = (if mindlin_force_history { sz } else { kt * sz }) + gamma_t * vt_z;

            let f_t_mag = (ft_x * ft_x + ft_y * ft_y + ft_z * ft_z).sqrt();
            if f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
                let scale = f_t_max / f_t_mag;
                ft_x *= scale;
                ft_y *= scale;
                ft_z *= scale;
            }

            let (sx, sy, sz) =
                if mindlin_force_history && f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON {
                    (
                        ft_x - gamma_t * vt_x,
                        ft_y - gamma_t * vt_y,
                        ft_z - gamma_t * vt_z,
                    )
                } else {
                    (sx, sy, sz)
                };

            // Torques
            let ti_x = r1n_y * ft_z - r1n_z * ft_y;
            let ti_y = r1n_z * ft_x - r1n_x * ft_z;
            let ti_z = r1n_x * ft_y - r1n_y * ft_x;
            let tj_x = r2n_y * ft_z - r2n_z * ft_y;
            let tj_y = r2n_z * ft_x - r2n_x * ft_z;
            let tj_z = r2n_x * ft_y - r2n_y * ft_x;

            atoms.force[i][0] += ft_x as soil_core::Accum;
            atoms.force[i][1] += ft_y as soil_core::Accum;
            atoms.force[i][2] += ft_z as soil_core::Accum;
            if newton {
                atoms.force[j][0] -= ft_x as soil_core::Accum;
                atoms.force[j][1] -= ft_y as soil_core::Accum;
                atoms.force[j][2] -= ft_z as soil_core::Accum;
            }
            dem.torque[i][0] += ti_x;
            dem.torque[i][1] += ti_y;
            dem.torque[i][2] += ti_z;
            if newton {
                dem.torque[j][0] += tj_x;
                dem.torque[j][1] += tj_y;
                dem.torque[j][2] += tj_z;
            }

            // Rolling/twisting relative angular velocity
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

            // Rolling resistance
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

                    dem.torque[i][0] += tr_x;
                    dem.torque[i][1] += tr_y;
                    dem.torque[i][2] += tr_z;
                    if newton {
                        dem.torque[j][0] -= tr_x;
                        dem.torque[j][1] -= tr_y;
                        dem.torque[j][2] -= tr_z;
                    }
                } else if roll_mag > 1e-30 {
                    let tau_mag = mu_r * f_n_mag.abs() * r_eff;
                    let inv_roll = tau_mag / roll_mag;
                    let tr_x = -inv_roll * roll_x;
                    let tr_y = -inv_roll * roll_y;
                    let tr_z = -inv_roll * roll_z;
                    dem.torque[i][0] += tr_x;
                    dem.torque[i][1] += tr_y;
                    dem.torque[i][2] += tr_z;
                    if newton {
                        dem.torque[j][0] -= tr_x;
                        dem.torque[j][1] -= tr_y;
                        dem.torque[j][2] -= tr_z;
                    }
                }
            }

            // Twisting friction (see the Hertz-Mindlin path for model semantics).
            if material_table.twisting_model == "marshall" {
                // Marshall (2009) derived-coefficient twisting on the linear (Hooke)
                // tangential model: k_twist = ½ k_t a², γ_twist = ½ γ_t a²,
                // μ_twist = (2/3) a μ_t with a = √(R* δ) and k_t = kt, γ_t the Hooke
                // tangential spring/damping computed above.
                let twist_vel = or_dot_n;
                let a_sq = delta * r_eff;
                let a = a_sq.sqrt();
                let k_twist = 0.5 * kt * a_sq;
                let gamma_twist = 0.5 * gamma_t * a_sq;
                let mu_twist = (2.0 / 3.0) * a * mu;

                twist_disp += twist_vel * dt;

                let mut tau_twist = -k_twist * twist_disp - gamma_twist * twist_vel;
                let tau_max = mu_twist * f_n_mag.abs();
                if tau_twist.abs() > tau_max {
                    tau_twist = tau_twist.signum() * tau_max;
                    if k_twist > TANGENTIAL_EPSILON {
                        twist_disp = (tau_twist + gamma_twist * twist_vel) / (-k_twist);
                    }
                }

                let tt_x = tau_twist * nx;
                let tt_y = tau_twist * ny;
                let tt_z = tau_twist * nz;
                dem.torque[i][0] += tt_x;
                dem.torque[i][1] += tt_y;
                dem.torque[i][2] += tt_z;
                if newton {
                    dem.torque[j][0] -= tt_x;
                    dem.torque[j][1] -= tt_y;
                    dem.torque[j][2] -= tt_z;
                }
            } else if mu_tw > 0.0 {
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

                    let tt_x = tau_twist * nx;
                    let tt_y = tau_twist * ny;
                    let tt_z = tau_twist * nz;
                    dem.torque[i][0] += tt_x;
                    dem.torque[i][1] += tt_y;
                    dem.torque[i][2] += tt_z;
                    if newton {
                        dem.torque[j][0] -= tt_x;
                        dem.torque[j][1] -= tt_y;
                        dem.torque[j][2] -= tt_z;
                    }
                } else if twist_vel.abs() > 1e-30 {
                    let tau = mu_tw * f_n_mag.abs() * r_eff;
                    let sign_tw = if twist_vel > 0.0 { -1.0 } else { 1.0 };
                    let tt_x = sign_tw * tau * nx;
                    let tt_y = sign_tw * tau * ny;
                    let tt_z = sign_tw * tau * nz;
                    dem.torque[i][0] += tt_x;
                    dem.torque[i][1] += tt_y;
                    dem.torque[i][2] += tt_z;
                    if newton {
                        dem.torque[j][0] -= tt_x;
                        dem.torque[j][1] -= tt_y;
                        dem.torque[j][2] -= tt_z;
                    }
                }
            }

            // Virial
            if let Some(ref mut v) = virial {
                if v.active {
                    let vs = if newton { 1.0 } else { 0.5 };
                    let vfx = (-fn_x + ft_x) * vs;
                    let vfy = (-fn_y + ft_y) * vs;
                    let vfz = (-fn_z + ft_z) * vs;
                    v.add_pair(dx, dy, dz, vfx, vfy, vfz);
                }
            }

            let mut new_spring = stored;
            new_spring[0] = sign * sx;
            new_spring[1] = sign * sy;
            new_spring[2] = sign * sz;
            new_spring[3] = sign * roll_disp_x;
            new_spring[4] = sign * roll_disp_y;
            new_spring[5] = sign * roll_disp_z;
            new_spring[6] = sign * twist_disp;
            new_spring[7] = contact_radius;
            match entry_idx {
                Some(idx) => {
                    history.contacts[i][idx].1 = new_spring;
                    history.contacts[i][idx].2 = true;
                }
                None => history.contacts[i].push((tag_j, new_spring, true)),
            }
        }

        for i in 0..nlocal {
            history.contacts[i].retain(|(_, _, active)| *active);
        }
    });
}
