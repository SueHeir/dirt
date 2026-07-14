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
