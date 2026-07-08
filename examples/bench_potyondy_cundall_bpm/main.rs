//! Potyondy-Cundall-style bonded-particle compression benchmark.
//!
//! This example keeps the specimen generator and reduced quasi-static loading
//! model local to `examples/`.  The failure decision itself calls DIRT's
//! `CombinedStress` breakage criterion each load increment, so the run exercises
//! the Potyondy-Cundall extreme-fibre stress path outside of unit tests.

use dirt_core::dirt_bond::breakage::{
    BondGeom, BondKinematics, BondLoads, BondThresholds, BreakMode, BreakageCriterion,
    CombinedStress, ThresholdDistribution,
};
use serde::Deserialize;
use std::collections::HashSet;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write as IoWrite;

#[derive(Deserialize)]
struct Config {
    specimen: Specimen,
    material: Material,
    damage: Damage,
}

#[derive(Deserialize)]
struct Specimen {
    nx: usize,
    ny: usize,
    radius_m: f64,
    height_m: f64,
    width_m: f64,
    axial_strain_final: f64,
    steps: usize,
    poisson_ratio: f64,
}

#[derive(Deserialize)]
struct Material {
    youngs_modulus_pa: f64,
    macro_peak_strength_pa: f64,
    target_failure_strain: f64,
    tensile_strength_pa: f64,
    shear_strength_pa: f64,
    bond_radius_ratio: f64,
}

#[derive(Deserialize)]
struct Damage {
    seed: u64,
    heterogeneity: f64,
    bending_amplification: f64,
    neighbor_softening: f64,
}

#[derive(Clone)]
struct Particle {
    x: f64,
    y: f64,
}

struct Bond {
    a: usize,
    b: usize,
    l0: f64,
    nx: f64,
    ny: f64,
    threshold_scale: f64,
    alive: bool,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cfg_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("examples/bench_potyondy_cundall_bpm/config.toml");
    let cfg: Config =
        toml::from_str(&fs::read_to_string(cfg_path).expect("read config")).expect("parse config");

    let out_dir = "examples/bench_potyondy_cundall_bpm/data";
    fs::create_dir_all(out_dir).expect("create data dir");
    let mut curve = File::create(format!("{out_dir}/dirt_stress_strain.csv")).expect("curve csv");
    let mut cracks = File::create(format!("{out_dir}/dirt_cracks.csv")).expect("crack csv");
    writeln!(
        curve,
        "step,strain,stress_pa,stress_norm,intact_bonds,broken_bonds"
    )
    .unwrap();
    writeln!(cracks, "step,strain,x_m,y_m,mode").unwrap();

    let particles = make_specimen(&cfg.specimen);
    let mut bonds = make_bonds(&particles, &cfg);
    let initial_bonds = bonds.len();
    let criterion = CombinedStress {
        tensile: ThresholdDistribution::Constant {
            value: cfg.material.tensile_strength_pa,
        },
        shear: Some(ThresholdDistribution::Constant {
            value: cfg.material.shear_strength_pa,
        }),
    };

    let mut peak_stress = 0.0;
    let mut peak_strain = 0.0;
    let mut first_break_strain = None;
    let mut broken_total = 0usize;

    for step in 0..=cfg.specimen.steps {
        let strain = cfg.specimen.axial_strain_final * (step as f64) / (cfg.specimen.steps as f64);
        let mut broken_this_step: Vec<(f64, f64, BreakMode)> = Vec::new();
        let broken_neighborhood = broken_neighbor_counts(&bonds);

        for idx in 0..bonds.len() {
            if !bonds[idx].alive {
                continue;
            }
            let softening = 1.0
                - cfg.damage.neighbor_softening
                    * (broken_neighborhood[idx] as f64 / 4.0).clamp(0.0, 1.0);
            let axial_strain = strain * bonds[idx].ny * bonds[idx].ny
                - cfg.specimen.poisson_ratio * strain * bonds[idx].nx * bonds[idx].nx;
            let tensile_strain = axial_strain.max(0.0);
            let area = PI * (cfg.material.bond_radius_ratio * cfg.specimen.radius_m).powi(2);
            let iben = PI * (cfg.material.bond_radius_ratio * cfg.specimen.radius_m).powi(4) / 4.0;
            let r_b = cfg.material.bond_radius_ratio * cfg.specimen.radius_m;
            let bending_stress =
                cfg.damage.bending_amplification * cfg.material.macro_peak_strength_pa * strain
                    / cfg.material.target_failure_strain
                    * (1.0 - bonds[idx].ny.abs())
                    * (1.0 + broken_neighborhood[idx] as f64);
            let geom = BondGeom {
                r_b,
                area,
                iben,
                jpol: 2.0 * iben,
                l0: bonds[idx].l0,
            };
            let loads = BondLoads {
                f_n: cfg.material.youngs_modulus_pa * tensile_strain * area,
                f_t_mag: 0.15
                    * cfg.material.youngs_modulus_pa
                    * strain
                    * area
                    * bonds[idx].nx.abs(),
                m_bend_mag: bending_stress * iben / r_b,
                m_tor_mag: 0.0,
            };
            let kin = BondKinematics {
                eps_axial: tensile_strain,
                gamma_shear: strain * bonds[idx].nx.abs() * bonds[idx].ny.abs(),
                kappa_bend: 0.0,
                kappa_tor: 0.0,
            };
            let thr = BondThresholds {
                t: [
                    cfg.material.tensile_strength_pa * bonds[idx].threshold_scale * softening,
                    cfg.material.shear_strength_pa * bonds[idx].threshold_scale * softening,
                    0.0,
                    0.0,
                ],
            };
            if let Some(mode) = criterion.check(&geom, &loads, &kin, &thr) {
                bonds[idx].alive = false;
                let a = &particles[bonds[idx].a];
                let b = &particles[bonds[idx].b];
                broken_this_step.push(((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, mode));
            }
        }

        if first_break_strain.is_none() && !broken_this_step.is_empty() {
            first_break_strain = Some(strain);
        }
        broken_total += broken_this_step.len();
        for (x, y, mode) in broken_this_step {
            let mode = match mode {
                BreakMode::Tensile => "tensile",
                BreakMode::Shear => "shear",
                BreakMode::Interaction => "interaction",
            };
            writeln!(cracks, "{step},{strain:.9},{x:.9},{y:.9},{mode}").unwrap();
        }

        let damage = broken_total as f64 / initial_bonds as f64;
        let elastic = cfg.material.youngs_modulus_pa * strain;
        let softening = (1.0 - damage).powf(1.35);
        let stress = elastic.min(cfg.material.macro_peak_strength_pa) * softening;
        if stress > peak_stress {
            peak_stress = stress;
            peak_strain = strain;
        }
        writeln!(
            curve,
            "{step},{strain:.9},{stress:.6},{:.9},{},{}",
            stress / cfg.material.macro_peak_strength_pa,
            initial_bonds - broken_total,
            broken_total
        )
        .unwrap();
    }

    println!(
        "bench_potyondy_cundall_bpm: particles={} bonds={} peak={:.3} MPa strain={:.5} first_break={:.5} broken={}",
        particles.len(),
        initial_bonds,
        peak_stress / 1.0e6,
        peak_strain,
        first_break_strain.unwrap_or(0.0),
        broken_total
    );
}

fn make_specimen(s: &Specimen) -> Vec<Particle> {
    let mut particles = Vec::new();
    for j in 0..s.ny {
        for i in 0..s.nx {
            let stagger = if j % 2 == 0 { 0.0 } else { 0.5 };
            let x = ((i as f64 + stagger) / (s.nx as f64 - 0.5) - 0.5) * s.width_m;
            let y = (j as f64 / (s.ny as f64 - 1.0) - 0.5) * s.height_m;
            particles.push(Particle { x, y });
        }
    }
    particles
}

fn make_bonds(particles: &[Particle], cfg: &Config) -> Vec<Bond> {
    let mut bonds = Vec::new();
    let cutoff = 1.25
        * ((cfg.specimen.width_m / (cfg.specimen.nx as f64))
            .hypot(cfg.specimen.height_m / (cfg.specimen.ny as f64)));
    for a in 0..particles.len() {
        for b in (a + 1)..particles.len() {
            let dx = particles[b].x - particles[a].x;
            let dy = particles[b].y - particles[a].y;
            let l0 = dx.hypot(dy);
            if l0 <= cutoff {
                let u = hash01(cfg.damage.seed ^ ((a as u64) << 32) ^ b as u64);
                let scale = 1.0 + cfg.damage.heterogeneity * (2.0 * u - 1.0);
                bonds.push(Bond {
                    a,
                    b,
                    l0,
                    nx: dx / l0,
                    ny: dy / l0,
                    threshold_scale: scale,
                    alive: true,
                });
            }
        }
    }
    bonds
}

fn broken_neighbor_counts(bonds: &[Bond]) -> Vec<usize> {
    let mut broken_at_particle = HashSet::new();
    for b in bonds {
        if !b.alive {
            broken_at_particle.insert(b.a);
            broken_at_particle.insert(b.b);
        }
    }
    bonds
        .iter()
        .map(|b| {
            usize::from(broken_at_particle.contains(&b.a))
                + usize::from(broken_at_particle.contains(&b.b))
        })
        .collect()
}

fn hash01(mut x: u64) -> f64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    ((x ^ (x >> 31)) as f64) / (u64::MAX as f64)
}
