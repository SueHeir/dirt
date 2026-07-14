//! Typed scheduler seams shared by DIRT plugins.
//!
//! The key type makes DIRT-owned producer/consumer boundaries explicit while
//! preserving the stable names emitted by scheduler diagnostics.

#![deny(missing_docs)]

use grass_scheduler::prelude::{SystemKey, SystemLabel};

macro_rules! seam {
    ($marker:ident, $key:ident, $name:literal, $doc:literal) => {
        #[doc = $doc]
        pub struct $marker;
        impl SystemLabel for $marker {
            const NAME: &'static str = $name;
        }
        #[doc = $doc]
        pub const $key: SystemKey<$marker> = SystemKey::new();
    };
}

seam!(
    ContactForce,
    CONTACT_FORCE,
    "hertz_mindlin_contact",
    "Particle-contact force provider."
);
seam!(
    WallContact,
    WALL_CONTACT,
    "wall_contact",
    "Wall-contact force provider."
);
seam!(
    BondForce,
    BOND_FORCE,
    "dem_bond_force",
    "Bond-force provider."
);
seam!(
    ClumpInsert,
    CLUMP_INSERT,
    "clump_insert_atoms",
    "Clump insertion setup seam."
);
seam!(
    ClumpGhostCutoff,
    CLUMP_GHOST_CUTOFF,
    "extend_ghost_cutoff_for_clumps",
    "Clump ghost-cutoff setup seam."
);
seam!(
    ClumpSnap,
    CLUMP_SNAP,
    "snap_subspheres_to_body_com",
    "Clump pre-exchange snap seam."
);
seam!(
    ClumpExchange,
    CLUMP_EXCHANGE,
    "exchange_bodies",
    "Clump body-exchange seam."
);
seam!(
    ClumpRestore,
    CLUMP_RESTORE,
    "restore_subsphere_positions",
    "Clump post-exchange restore seam."
);
seam!(
    ClumpRemap,
    CLUMP_REMAP,
    "remap_bodies_on_box_resize",
    "Clump box-resize remap seam."
);
seam!(
    ClumpInitialIntegration,
    CLUMP_INITIAL_INTEGRATION,
    "integrate_bodies_initial",
    "Clump initial-integration seam."
);
seam!(
    ClumpPbc,
    CLUMP_PBC,
    "pbc_multisphere_bodies",
    "Clump periodic-boundary seam."
);
seam!(
    ClumpForceAggregation,
    CLUMP_FORCE_AGGREGATION,
    "aggregate_clump_forces",
    "Clump force-aggregation seam."
);
seam!(
    ClumpFinalIntegration,
    CLUMP_FINAL_INTEGRATION,
    "integrate_bodies_final",
    "Clump final-integration seam."
);
seam!(
    ClumpPreExchangeUpdate,
    CLUMP_PRE_EXCHANGE_UPDATE,
    "update_clump_positions_pre_exchange",
    "Clump pre-exchange update seam."
);
seam!(
    ClumpPositionUpdate,
    CLUMP_POSITION_UPDATE,
    "update_clump_positions",
    "Clump position-update seam."
);
seam!(
    ClumpLostAtomCheck,
    CLUMP_LOST_ATOM_CHECK,
    "check_lost_clump_atoms",
    "Clump lost-atom check seam."
);
seam!(
    ContactAnalysis,
    CONTACT_ANALYSIS,
    "contact_analysis",
    "Contact-analysis seam."
);
seam!(
    MeasurePlaneElapsed,
    MEASURE_PLANE_ELAPSED,
    "measure_plane_accumulate_elapsed",
    "Measurement-plane elapsed-time seam."
);
seam!(
    MeasurePlaneReport,
    MEASURE_PLANE_REPORT,
    "measure_plane_report",
    "Measurement-plane reporting seam."
);
seam!(
    MeasurePlaneCrossings,
    MEASURE_PLANE_CROSSINGS,
    "measure_plane_detect_crossings",
    "Measurement-plane crossing seam."
);
seam!(
    AutoBond,
    AUTO_BOND,
    "auto_bond_touching",
    "Automatic bond-creation seam."
);
seam!(
    LoadBonds,
    LOAD_BONDS,
    "load_bonds_from_file",
    "Bond-file loading seam."
);
seam!(
    BondGhostCutoff,
    BOND_GHOST_CUTOFF,
    "extend_ghost_cutoff_for_bonds",
    "Bond ghost-cutoff seam."
);
seam!(
    BondBreakageInit,
    BOND_BREAKAGE_INIT,
    "init_breakage",
    "Bond breakage-initialization seam."
);
seam!(
    BondPlasticityInit,
    BOND_PLASTICITY_INIT,
    "init_plasticity",
    "Bond plasticity-initialization seam."
);
seam!(
    AddForce,
    ADD_FORCE,
    "dirt_fixes::addforce",
    "Add-force fix seam."
);
seam!(
    SetForce,
    SET_FORCE,
    "dirt_fixes::setforce",
    "Set-force fix seam."
);
