use std::env;
use std::error::Error;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

use dirt_clump::{
    compute_inertia_tensor_montecarlo, compute_inertia_tensor_montecarlo_seeded, ClumpSphereConfig,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Config {
    density: f64,
    radius: f64,
    tolerance_rel: f64,
    repeat_count: usize,
    sample_counts: Vec<usize>,
    seeds: Vec<u64>,
}

fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1.0e-30)
}

fn max_diag_rel_err(tensor: [[f64; 3]; 3], reference_i: f64) -> f64 {
    (0..3)
        .map(|d| rel_err(tensor[d][d], reference_i))
        .fold(0.0, f64::max)
}

fn default_data_dir(config_path: &Path) -> PathBuf {
    config_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("data")
}

fn main() -> Result<(), Box<dyn Error>> {
    let config_path = env::args()
        .nth(1)
        .unwrap_or_else(|| "examples/bench_clump_inertia_sampler/config.toml".to_string());
    let config_path = PathBuf::from(config_path);
    let config: Config = toml::from_str(&fs::read_to_string(&config_path)?)?;

    if config.sample_counts.is_empty() || config.seeds.is_empty() || config.repeat_count == 0 {
        return Err("sample_counts, seeds, and repeat_count must be non-empty".into());
    }

    let data_dir = default_data_dir(&config_path);
    fs::create_dir_all(&data_dir)?;
    let csv_path = data_dir.join("inertia_sampler.csv");
    let mut csv = File::create(&csv_path)?;

    let spheres = [ClumpSphereConfig {
        offset: [0.0, 0.0, 0.0],
        radius: config.radius,
    }];
    let expected_mass = config.density * (4.0 / 3.0) * PI * config.radius.powi(3);
    let expected_i = 0.4 * expected_mass * config.radius * config.radius;
    let check_samples = *config
        .sample_counts
        .last()
        .expect("sample_counts is non-empty");

    writeln!(
        csv,
        "mode,sample_count,seed,repeat,mass,ixx,iyy,izz,mass_rel_err,max_diag_rel_err,bitwise_repeat"
    )?;

    let default_first = compute_inertia_tensor_montecarlo(&spheres, config.density, check_samples);
    for repeat in 0..config.repeat_count {
        let (mass, tensor) =
            compute_inertia_tensor_montecarlo(&spheres, config.density, check_samples);
        let bitwise_repeat = mass.to_bits() == default_first.0.to_bits()
            && (0..3)
                .all(|a| (0..3).all(|b| tensor[a][b].to_bits() == default_first.1[a][b].to_bits()));
        writeln!(
            csv,
            "default_repeat,{check_samples},,{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}",
            repeat,
            mass,
            tensor[0][0],
            tensor[1][1],
            tensor[2][2],
            rel_err(mass, expected_mass),
            max_diag_rel_err(tensor, expected_i),
            bitwise_repeat
        )?;
    }

    let seed = config.seeds[0];
    let seeded_first =
        compute_inertia_tensor_montecarlo_seeded(&spheres, config.density, check_samples, seed);
    for repeat in 0..config.repeat_count {
        let (mass, tensor) =
            compute_inertia_tensor_montecarlo_seeded(&spheres, config.density, check_samples, seed);
        let bitwise_repeat = mass.to_bits() == seeded_first.0.to_bits()
            && (0..3)
                .all(|a| (0..3).all(|b| tensor[a][b].to_bits() == seeded_first.1[a][b].to_bits()));
        writeln!(
            csv,
            "explicit_seed_repeat,{check_samples},{seed},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}",
            repeat,
            mass,
            tensor[0][0],
            tensor[1][1],
            tensor[2][2],
            rel_err(mass, expected_mass),
            max_diag_rel_err(tensor, expected_i),
            bitwise_repeat
        )?;
    }

    for &sample_count in &config.sample_counts {
        for &seed in &config.seeds {
            let (mass, tensor) = compute_inertia_tensor_montecarlo_seeded(
                &spheres,
                config.density,
                sample_count,
                seed,
            );
            writeln!(
                csv,
                "seed_spread,{sample_count},{seed},,{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},true",
                mass,
                tensor[0][0],
                tensor[1][1],
                tensor[2][2],
                rel_err(mass, expected_mass),
                max_diag_rel_err(tensor, expected_i),
            )?;
        }
    }

    println!("expected_mass={expected_mass:.12e}");
    println!("expected_i={expected_i:.12e}");
    println!("tolerance_rel={:.6e}", config.tolerance_rel);
    println!("wrote {}", csv_path.display());
    Ok(())
}
