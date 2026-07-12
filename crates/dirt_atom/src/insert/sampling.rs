use super::*;

/// Samples one candidate point, treating SOIL's bounded rejection-sampling
/// exhaustion as a rejected insertion attempt.
///
/// A composite region can be structurally valid while its overlap has such a
/// small measure that SOIL exhausts its finite inner rejection budget. No
/// particle can be constructed in that case. The outer insertion loop already
/// has a deterministic attempt limit, so retrying there is the same physical
/// outcome as an overlap rejection: successful draws keep the pre-existing RNG
/// stream and MPI ownership rules, while an exhausted draw cannot panic or
/// terminate just one simulation rank.

pub(super) fn try_sample_insertion_point(
    region: &Region,
    rng: &mut impl rand::Rng,
) -> Option<[f64; 3]> {
    region.random_point_inside(rng).ok()
}

// ── SpatialHash for O(1) overlap checking ───────────────────────────────────

/// Grid-based spatial hash for fast overlap detection during particle insertion.
///
/// Divides space into cubic cells of size `cell_size` (typically ~2× max particle diameter).
/// Overlap queries check the 3×3×3 neighborhood of the candidate cell, ensuring all
/// potential overlaps are found without a full O(N²) scan.
/// Periodic boundary info for minimum-image overlap checks during insertion.
pub(super) struct PeriodicBox {
    pub(super) is_periodic: [bool; 3],
    pub(super) box_size: [f64; 3],
}

impl PeriodicBox {
    pub(super) fn from_domain(domain: &Domain) -> Self {
        PeriodicBox {
            is_periodic: domain.periodic_flags(),
            box_size: domain.size,
        }
    }

    /// Compute minimum-image squared distance between two positions.
    pub(super) fn min_image_dist_sq(&self, a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let mut dist_sq = 0.0;
        for d in 0..3 {
            let mut delta = a[d] - b[d];
            if self.is_periodic[d] {
                let half = 0.5 * self.box_size[d];
                if delta > half {
                    delta -= self.box_size[d];
                } else if delta < -half {
                    delta += self.box_size[d];
                }
            }
            dist_sq += delta * delta;
        }
        dist_sq
    }
}

pub(super) struct SpatialHash {
    cell_size: f64,
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl SpatialHash {
    pub(super) fn new(cell_size: f64) -> Self {
        SpatialHash {
            cell_size,
            cells: HashMap::new(),
        }
    }

    pub(super) fn cell_key(&self, pos: &[f64; 3]) -> (i64, i64, i64) {
        (
            (pos[0] / self.cell_size).floor() as i64,
            (pos[1] / self.cell_size).floor() as i64,
            (pos[2] / self.cell_size).floor() as i64,
        )
    }

    pub(super) fn insert(&mut self, idx: usize, pos: &[f64; 3]) {
        let key = self.cell_key(pos);
        self.cells.entry(key).or_default().push(idx);
    }

    pub(super) fn has_overlap(
        &self,
        pos: &[f64; 3],
        radius: f64,
        positions: &[[f64; 3]],
        radii: &[f64],
        pbc: &PeriodicBox,
    ) -> bool {
        // Collect all cell keys to check: 3x3x3 neighborhood plus periodic images
        let key = self.cell_key(pos);
        let min_dist_check = radius * 2.2; // conservative check radius

        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    let neighbor_key = (key.0 + di, key.1 + dj, key.2 + dk);
                    if let Some(indices) = self.cells.get(&neighbor_key) {
                        for &idx in indices {
                            let dist_sq = pbc.min_image_dist_sq(pos, &positions[idx]);
                            let min_dist = (radius + radii[idx]) * 1.1;
                            if dist_sq <= min_dist * min_dist {
                                return true;
                            }
                        }
                    }
                }
            }
        }

        // For periodic axes where the box is small (< 3 cell sizes), the standard
        // 3x3x3 neighborhood may miss periodic images. Do a brute-force check
        // against all atoms using minimum-image distances.
        let needs_pbc_check = (0..3)
            .any(|d| pbc.is_periodic[d] && pbc.box_size[d] < 3.0 * self.cell_size + min_dist_check);
        if needs_pbc_check {
            for idx in 0..positions.len() {
                let dist_sq = pbc.min_image_dist_sq(pos, &positions[idx]);
                let min_dist = (radius + radii[idx]) * 1.1;
                if dist_sq <= min_dist * min_dist {
                    return true;
                }
            }
        }

        false
    }
}

/// One globally replicated accepted candidate.  Every rank constructs this
/// record before ownership filtering, so RNG consumption and tag assignment are
/// independent of the domain decomposition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct InsertCandidate {
    pub(super) pos: [f64; 3],
    pub(super) radius: f64,
    pub(super) vel: [f64; 3],
}

impl InsertCandidate {
    pub(super) fn particle(self, prepared: &PreparedInsert, tag: u32) -> DemParticle {
        DemParticle {
            pos: self.pos,
            vel: self.vel,
            radius: self.radius,
            cutoff_padding: prepared.cutoff_padding,
            density: prepared.density,
            mat_idx: prepared.mat_idx,
            tag,
        }
    }
}

/// Shared deterministic candidate stream for both random insertion schedules.
/// A call consumes the same random draws on every rank and returns `None` for a
/// rejected attempt; callers deliberately keep their distinct tag policies.
pub(super) struct CandidateGenerator {
    rng: StdRng,
    spatial_hash: SpatialHash,
    positions: Vec<[f64; 3]>,
    radii: Vec<f64>,
    pbc: PeriodicBox,
}

impl CandidateGenerator {
    pub(super) fn new(
        prepared: &PreparedInsert,
        domain: &Domain,
        seed: u64,
        capacity: usize,
    ) -> Self {
        Self {
            rng: StdRng::seed_from_u64(seed),
            spatial_hash: SpatialHash::new((2.0 * prepared.max_radius * 1.1).max(1e-10)),
            positions: Vec::with_capacity(capacity),
            radii: Vec::with_capacity(capacity),
            pbc: PeriodicBox::from_domain(domain),
        }
    }

    pub(super) fn next(&mut self, prepared: &PreparedInsert) -> Option<InsertCandidate> {
        let pos = try_sample_insertion_point(&prepared.region, &mut self.rng)?;
        let radius = prepared
            .radius
            .try_sample(&mut self.rng)
            .expect("PreparedInsert only contains a validated radius specification");
        if self
            .spatial_hash
            .has_overlap(&pos, radius, &self.positions, &self.radii, &self.pbc)
        {
            return None;
        }
        let mut vel = prepared.velocity_base;
        if let Some(normal) = &prepared.velocity_normal {
            vel[0] += normal.sample(&mut self.rng);
            vel[1] += normal.sample(&mut self.rng);
            vel[2] += normal.sample(&mut self.rng);
        }
        let index = self.positions.len();
        self.spatial_hash.insert(index, &pos);
        self.positions.push(pos);
        self.radii.push(radius);
        Some(InsertCandidate { pos, radius, vel })
    }
}
