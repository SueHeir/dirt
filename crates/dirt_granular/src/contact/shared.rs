use super::*;
pub(super) fn is_mindlin_rescale(model: &str) -> bool {
    matches!(
        model,
        "mindlin_rescale" | "mindlin_rescale_force" | "mindlin_rescale/force"
    )
}

pub(super) fn is_mindlin_force_history(model: &str) -> bool {
    matches!(model, "mindlin_rescale_force" | "mindlin_rescale/force")
}

pub(super) fn zero_contact_history() -> ContactHistory {
    [0.0; CONTACT_HISTORY_LEN]
}

/// First slot of the MDR per-side block for the `i` side; see `mdr.rs`.
const MDR_SIDE0: usize = 9;
/// First slot of the MDR per-side block for the `j` side.
const MDR_SIDE1: usize = 16;
/// Length of one MDR per-side block.
const MDR_SIDE_LEN: usize = 7;

/// Find the stored history for pair `(i, j)`.
///
/// The half neighbor list holds each pair once, but which atom is `i` can
/// change at a rebuild: atoms change bins, a spatial sort reorders indices,
/// or a pair crossing a sub-domain face moves to the other rank. The history
/// lives in the list of whichever atom was `i` when it was written, so look
/// there too:
///
/// - `contacts[i]` (same orientation as last step): returned with its index
///   for in-place write-back.
/// - `contacts[j]` (pair flipped): removed from `j`'s list and returned with
///   index `None`, so the caller pushes it onto `contacts[i]`. For a ghost `j`
///   this is the owner's list as of the last border exchange, which is when
///   the flip happened; the owner's copy is pruned there because the pair is
///   no longer listed on that rank.
///
/// Slots 0..8 are stored canonically (lower-tag frame) and need no change;
/// the MDR per-side blocks are in the `(i, j)` frame, so they are swapped.
pub(super) fn find_contact_history(
    history: &mut ContactHistoryStore,
    i: usize,
    j: usize,
    tag_i: u32,
    tag_j: u32,
) -> (Option<usize>, ContactHistory) {
    if let Some(idx) = history.contacts[i].iter().position(|(t, _, _)| *t == tag_j) {
        return (Some(idx), history.contacts[i][idx].1);
    }
    if let Some(list) = history.contacts.get_mut(j) {
        if let Some(k) = list.iter().position(|(t, _, _)| *t == tag_i) {
            let (_, mut stored, _) = list.swap_remove(k);
            for s in 0..MDR_SIDE_LEN {
                stored.swap(MDR_SIDE0 + s, MDR_SIDE1 + s);
            }
            return (None, stored);
        }
    }
    (None, zero_contact_history())
}

/// Which pairs to process during an overlapped force calculation.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ForcePass {
    /// Process every pair in one pass.
    All,
    /// Process local-local pairs while communication is in flight.
    Interior,
    /// Process ghost-boundary pairs after communication completes.
    Boundary,
}
