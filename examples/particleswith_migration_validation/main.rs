//! Entry point for the ParticlesWith migration compatibility evidence.
//!
//! The measured before/after comparison is intentionally driven by `sweep.py`:
//! it runs the representative contact/wall, bond, and clump cases on pinned
//! source revisions and checks their emitted CSVs byte-for-byte.

fn main() {
    println!("Run `python3 examples/particleswith_migration_validation/sweep.py graph`.");
}
