//! contact_analysis_demo — coordination number, fabric tensor, and per-contact CSV.
//!
//! Drops a cloud of spheres into a closed box (the `hello_bed` setup) and adds
//! [`ContactAnalysisPlugin`] so the settled packing is characterized as it
//! forms. The `[contact_analysis]` config section turns on coordination number,
//! rattler detection, the fabric tensor, and a periodic per-contact CSV dump.
//!
//! ```bash
//! cargo run --release --example contact_analysis_demo --no-default-features -- \
//!     examples/contact_analysis_demo/config.toml
//! ```
//!
//! # Reading the output
//!
//! **Thermo columns** (printed every `thermo` steps to stdout / the log):
//!
//! - `coord_avg`, `coord_max`, `coord_min` — coordination-number statistics. As
//!   the cloud settles, `coord_avg` climbs from ~0 toward the random-close-pack
//!   value (~6 for frictional spheres).
//! - `n_rattlers`, `rattler_fraction` — particles with fewer than 4 contacts
//!   (mechanically under-constrained in 3D).
//! - `fabric_xx`, `fabric_yy`, `fabric_zz`, `fabric_xy`, `fabric_xz`,
//!   `fabric_yz`, `contacts` — the six independent fabric-tensor components and
//!   total contact count. A gravity-settled bed is anisotropic, so
//!   `fabric_zz` typically differs from `fabric_xx`/`fabric_yy`.
//!
//! **Per-contact CSV** (because `interval > 0`): files are written to
//! `<output_dir>/contact/contact_<step>_rank<rank>.csv` with columns
//! `i_tag, j_tag, overlap, cx, cy, cz, nx, ny, nz` — one row per contacting
//! pair, recorded once per physical contact regardless of the neighbor-list
//! Newton setting. Load one in Python with `numpy.loadtxt(..., delimiter=",",
//! skiprows=1)` (or pandas) to plot the force-chain geometry or build a contact
//! network.

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins) // includes PrintPlugin (required before ContactAnalysisPlugin)
        .add_plugins(GranularDefaultPlugins) // Hertz-Mindlin contact (label "hertz_mindlin_contact")
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        // ContactAnalysisPlugin must come AFTER CorePlugins (DumpRegistry) and
        // the contact plugin (the "hertz_mindlin_contact" ordering anchor).
        .add_plugins(ContactAnalysisPlugin);

    app.start();
}
