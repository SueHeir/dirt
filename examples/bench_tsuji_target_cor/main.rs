//! Command-line access to DIRT's physical target-COR Hertz--Tsuji calibration.
//!
//! The maintained Python campaign uses this executable to obtain the raw Tsuji
//! material input, then independently validates the mapping against its own ODE
//! integration, DIRT impacts, and LAMMPS impacts.

use dirt_core::dirt_atom::{hertz_beta_for_cor, hertz_tsuji_raw_for_target_cor};

fn main() {
    let mut args = std::env::args().skip(1);
    if args.next().as_deref() != Some("calibrate") {
        eprintln!("usage: bench_tsuji_target_cor calibrate <target-cor>");
        std::process::exit(2);
    }
    let target = args
        .next()
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or_else(|| {
            eprintln!("target COR must be a finite number in [0, 1]");
            std::process::exit(2);
        });
    let Some(raw) = hertz_tsuji_raw_for_target_cor(target) else {
        eprintln!("target COR is outside the Hertz--Tsuji calibratable range");
        std::process::exit(2);
    };
    // CSV is intentionally machine-readable for the maintained sweep driver.
    println!("{target:.9},{raw:.9},{:.9}", hertz_beta_for_cor(raw));
}
