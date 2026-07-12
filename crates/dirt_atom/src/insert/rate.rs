use super::*;
// ── Update system: rate-based insertion ─────────────────────────────────────

/// Update system for rate-based particle insertion during the simulation run.
///
/// Checks each registered [`RateInsertEntry`] against the current timestep, interval,
/// start/end bounds, and total limit. Uses a [`SpatialHash`] for O(1) overlap detection
/// when placing new particles. Runs in `ParticleSimScheduleSet::PreInitialIntegration`.
#[allow(clippy::too_many_arguments)]
pub fn dem_rate_insert(
    comm: Res<CommResource>,
    domain: Res<Domain>,
    mut atom: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    mut rate_state: ResMut<RateInsertState>,
    mut comm_state: ResMut<CurrentState<CommState>>,
) {
    // Rate insertion runs on EVERY rank (born-in-owner): each rank generates the
    // identical candidate stream from a step-derived seed and stores only the
    // candidates that fall inside its own subdomain, so a new atom is born inside
    // its owner and never needs a multi-hop exchange.
    if rate_state.entries.is_empty() {
        return;
    }

    let step = run_state.total_cycle;
    let mut any_to_insert = false;

    // Quick check if any entry needs insertion this step (before stripping ghosts)

    // Quick check if any entry needs insertion this step
    for entry in rate_state.entries.iter() {
        let interval = entry.config.rate_interval.unwrap_or(1);
        let start = entry.config.rate_start.unwrap_or(0);
        if step < start {
            continue;
        }
        if let Some(end) = entry.config.rate_end {
            if step > end {
                continue;
            }
        }
        if let Some(limit) = entry.config.rate_limit {
            if entry.total_inserted >= limit {
                continue;
            }
        }
        let steps_since_start = step - start;
        if interval == 0 || steps_since_start % interval == 0 {
            any_to_insert = true;
            break;
        }
    }

    if !any_to_insert {
        return;
    }

    // Strip ghost atoms before inserting new local atoms.
    // New atoms are appended at atom.len(), which must equal nlocal so that
    // the subsequent borders() truncate_to_nlocal() doesn't discard them.
    if atom.nghost > 0 {
        // Keep the core ghost suffix and every DIRT extension in one
        // transaction.  In particular, a plugin registered after setup has a
        // row here too; truncating Atom and AtomDataRegistry separately would
        // leave a panic/error window between the two mutations.
        ParticleStore::new(&mut atom, &registry)
            .discard_ghosts()
            .expect("rate insertion requires an aligned local/ghost particle layout");
    }

    // Base tag must be globally consistent across ranks. (No all_reduce_max in
    // the backend, so reduce -max with min and negate.) Tags advance only for
    // accepted candidates, matching immediate insertion. Acceptance is globally
    // replicated, so this remains rank-count invariant without another
    // collective and rejected attempts cannot leave gaps in the accepted stream.
    let local_max_tag = if atom.tag.is_empty() {
        -1.0
    } else {
        atom.get_max_tag() as f64
    };
    let base_tag = (-comm.all_reduce_min_f64(-local_max_tag)) as i64 + 1;
    let mut tag_cursor: u32 = base_tag.max(0) as u32;

    for entry_idx in 0..rate_state.entries.len() {
        let interval = rate_state.entries[entry_idx]
            .config
            .rate_interval
            .unwrap_or(1);
        let start = rate_state.entries[entry_idx].config.rate_start.unwrap_or(0);
        let (rate, _, _) = validate_rate_insert_config(
            &rate_state.entries[entry_idx].config,
            "rate-based [[particles.insert]]",
        )
        .expect("rate insertion was validated during fallible plugin preflight");

        if step < start {
            continue;
        }
        if let Some(end) = rate_state.entries[entry_idx].config.rate_end {
            if step > end {
                continue;
            }
        }
        if let Some(limit) = rate_state.entries[entry_idx].config.rate_limit {
            if rate_state.entries[entry_idx].total_inserted >= limit {
                continue;
            }
        }
        let steps_since_start = step - start;
        if interval > 0 && steps_since_start % interval != 0 {
            continue;
        }

        // How many to insert this step
        let mut to_insert = rate;
        if let Some(limit) = rate_state.entries[entry_idx].config.rate_limit {
            let remaining = limit - rate_state.entries[entry_idx].total_inserted;
            to_insert = to_insert.min(remaining);
        }

        let prepared = &rate_state.entries[entry_idx].prepared;

        // Seed the candidate stream from (config seed, step, entry) so it is
        // identical on every rank yet varies between insertion events.
        let stream_seed = prepared.seed
            ^ (step as u64).wrapping_mul(0x9E3779B97F4A7C15)
            ^ (entry_idx as u64).wrapping_mul(0xD1B54A32D192ED03);

        // Replicated overlap scratch: ONLY the positions/radii accepted THIS step
        // (the global set of new atoms). It is identical on every rank because the
        // candidate stream is seeded identically and the scratch is the same on
        // all ranks. This is what keeps accept/reject — and therefore the RNG
        // advancement and the global accept count `inserted` — in lock-step across
        // ranks, so the collective borders()/exchange() triggered below stay
        // synchronized.
        //
        // NOTE: existing local atoms are deliberately NOT added to the scratch.
        // They differ per rank, so including them would make accept/reject (and
        // the RNG stream) diverge across ranks and desync the collectives. New
        // atoms may therefore be born overlapping already-present particles; the
        // contact model resolves that initial overlap via repulsion. Rate-insert
        // regions are normally placed in free space (e.g. above a settled bed),
        // so this is rare in practice.
        let mut candidates =
            CandidateGenerator::new(prepared, &domain, stream_seed, to_insert as usize);

        let mut inserted = 0u32; // accepted globally
        let mut local_inserted = 0u32; // stored on this rank
        let mut attempts = 0u32;
        let max_attempts = to_insert * 100;

        while inserted < to_insert && attempts < max_attempts {
            attempts += 1;

            let Some(candidate) = candidates.next(prepared) else {
                continue;
            };
            let tag = tag_cursor;
            tag_cursor = tag_cursor.wrapping_add(1);

            // Store only if this rank owns the position.
            if owns_position(&domain, &candidate.pos) {
                insert_single_particle(&mut atom, &registry, candidate.particle(prepared, tag));
                local_inserted += 1;
            }
            inserted += 1;
        }

        rate_state.entries[entry_idx].total_inserted += inserted;
        let _ = local_inserted;

        if inserted > 0 {
            // Force full ghost rebuild on EVERY rank if any rank inserted, so the
            // collective borders()/exchange() stay in lock-step. `inserted` is the
            // global accept count and is identical on all ranks (replicated stream).
            comm_state.0 = CommState::FullRebuild;
        }
        if inserted > 0 && attempts >= max_attempts && comm.rank() == 0 {
            eprintln!(
                "WARNING: Rate insertion at step {} only placed {}/{} particles (max attempts reached)",
                step, inserted, to_insert
            );
        }
    }
}
