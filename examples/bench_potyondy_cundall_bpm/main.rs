//! Potyondy-Cundall-style bonded-particle compression benchmark.
//!
//! The example runs a real DIRT bonded-particle specimen: particles are inserted
//! from a generated lattice file, `DemBondPlugin` creates and breaks
//! `CombinedStress` bonds, and this recorder derives axial stress and crack
//! progression from the live bonded network.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::dirt_bond::{BondConfig, BondMetrics};
use dirt_core::dirt_fixes::FixesPlugin;
use dirt_core::prelude::*;
use dirt_core::soil_core::BondStore;
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufWriter, Write as IoWrite};

#[derive(Clone, Deserialize)]
struct BenchSpec {
    nx: usize,
    ny: usize,
    spacing_m: f64,
    radius_m: f64,
    height_m: f64,
    width_m: f64,
    thickness_m: f64,
    imperfection_m: f64,
    compression_rate_m_s: f64,
    record_every: usize,
}

impl BenchSpec {
    fn from_argv() -> Self {
        let args: Vec<String> = std::env::args().collect();
        let cfg_path = args
            .get(1)
            .map(String::as_str)
            .unwrap_or("examples/bench_potyondy_cundall_bpm/config.toml");
        let text = fs::read_to_string(cfg_path).expect("read config");
        let value: toml::Value = text.parse().expect("parse config");
        let table = value
            .get("potyondy_specimen")
            .unwrap_or_else(|| panic!("{cfg_path} must contain [potyondy_specimen]"));
        table.clone().try_into().expect("parse [potyondy_specimen]")
    }
}

struct Recorder {
    spec: BenchSpec,
    curve: Option<BufWriter<File>>,
    cracks: Option<BufWriter<File>>,
    initial_pos: BTreeMap<u32, [f64; 3]>,
    top_tags: BTreeSet<u32>,
    bottom_tags: BTreeSet<u32>,
    live_bonds_last: BTreeSet<(u32, u32)>,
    initialized: bool,
}

impl Recorder {
    fn new(spec: BenchSpec) -> Self {
        Self {
            spec,
            curve: None,
            cracks: None,
            initial_pos: BTreeMap::new(),
            top_tags: BTreeSet::new(),
            bottom_tags: BTreeSet::new(),
            live_bonds_last: BTreeSet::new(),
            initialized: false,
        }
    }
}

fn main() {
    let spec = BenchSpec::from_argv();
    write_specimen_particles(&spec);

    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin)
        .add_plugins(DemBondPlugin);

    app.add_resource(Recorder::new(spec));
    app.add_update_system(
        record_compression_state,
        ParticleSimScheduleSet::PostFinalIntegration,
    );
    app.start();
}

fn write_specimen_particles(spec: &BenchSpec) {
    assert!(
        spec.spacing_m >= 2.0 * spec.radius_m,
        "specimen spacing must not overlap particles"
    );
    let out = "examples/bench_potyondy_cundall_bpm/data/specimen_particles.csv";
    fs::create_dir_all("examples/bench_potyondy_cundall_bpm/data").expect("create data dir");
    let mut w = BufWriter::new(File::create(out).expect("create specimen particle csv"));
    for j in 0..spec.ny {
        let y = -0.5 * spec.height_m + (j as f64) * spec.spacing_m;
        let stagger = if j % 2 == 0 {
            0.0
        } else {
            0.5 * spec.spacing_m
        };
        for i in 0..spec.nx {
            let x = -0.5 * spec.width_m + (i as f64) * spec.spacing_m + stagger;
            let z = spec.imperfection_m * ((i * 17 + j * 31) as f64).sin();
            writeln!(w, "{x:.9},{y:.9},{z:.9}").expect("write particle");
        }
    }
}

fn init_recorder(
    atoms: &Atom,
    registry: &AtomDataRegistry,
    input: &Input,
    recorder: &mut Recorder,
) -> bool {
    let Some(bonds) = registry.get::<BondStore>() else {
        return false;
    };
    if bonds.bonds.is_empty() {
        return false;
    }

    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_potyondy_cundall_bpm".to_string());
    fs::create_dir_all(format!("{out_dir}/data")).expect("create output data dir");
    let mut curve = BufWriter::new(
        File::create(format!("{out_dir}/data/dirt_stress_strain.csv")).expect("curve csv"),
    );
    let mut cracks =
        BufWriter::new(File::create(format!("{out_dir}/data/dirt_cracks.csv")).expect("crack csv"));
    writeln!(
        curve,
        "step,strain,stress_pa,stress_norm,intact_bonds,broken_bonds,crossing_bonds"
    )
    .unwrap();
    writeln!(cracks, "step,strain,x_m,y_m,z_m,mode").unwrap();

    let nlocal = atoms.nlocal as usize;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for i in 0..nlocal {
        y_min = y_min.min(atoms.pos[i][1] as f64);
        y_max = y_max.max(atoms.pos[i][1] as f64);
        recorder.initial_pos.insert(
            atoms.tag[i],
            [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ],
        );
    }
    let row_tol = 0.51 * recorder.spec.spacing_m;
    for i in 0..nlocal {
        let y = atoms.pos[i][1] as f64;
        if (y - y_min).abs() <= row_tol {
            recorder.bottom_tags.insert(atoms.tag[i]);
        }
        if (y - y_max).abs() <= row_tol {
            recorder.top_tags.insert(atoms.tag[i]);
        }
    }
    recorder.live_bonds_last = live_bond_set(atoms, &bonds);
    recorder.curve = Some(curve);
    recorder.cracks = Some(cracks);
    recorder.initialized = true;
    true
}

fn live_bond_set(atoms: &Atom, bonds: &BondStore) -> BTreeSet<(u32, u32)> {
    let mut set = BTreeSet::new();
    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal.min(bonds.bonds.len()) {
        let tag_i = atoms.tag[i];
        for b in &bonds.bonds[i] {
            let a = tag_i.min(b.partner_tag);
            let c = tag_i.max(b.partner_tag);
            set.insert((a, c));
        }
    }
    set
}

fn avg_y_for_tags(atoms: &Atom, tags: &BTreeSet<u32>) -> Option<f64> {
    let nlocal = atoms.nlocal as usize;
    let mut sum = 0.0;
    let mut n = 0usize;
    for i in 0..nlocal {
        if tags.contains(&atoms.tag[i]) {
            sum += atoms.pos[i][1] as f64;
            n += 1;
        }
    }
    (n > 0).then_some(sum / n as f64)
}

fn tag_index(atoms: &Atom, tag: u32) -> Option<usize> {
    (0..atoms.nlocal as usize).find(|&i| atoms.tag[i] == tag)
}

fn record_compression_state(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    bond_config: Res<BondConfig>,
    bond_metrics: Res<BondMetrics>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut recorder: ResMut<Recorder>,
) {
    if !recorder.initialized && !init_recorder(&atoms, &registry, &input, &mut recorder) {
        return;
    }

    let step = run_state.total_cycle;
    if step % recorder.spec.record_every != 0 {
        return;
    }

    let Some(bonds) = registry.get::<BondStore>() else {
        return;
    };
    let Some(dem) = registry.get::<DemAtom>() else {
        return;
    };

    let y_bottom =
        avg_y_for_tags(&atoms, &recorder.bottom_tags).unwrap_or(-0.5 * recorder.spec.height_m);
    let y_top = avg_y_for_tags(&atoms, &recorder.top_tags).unwrap_or(0.5 * recorder.spec.height_m);
    let strain = (step as f64 * atoms.dt * recorder.spec.compression_rate_m_s
        / recorder.spec.height_m)
        .max(0.0);
    let section_area = recorder.spec.width_m * recorder.spec.thickness_m;
    let mid_y = 0.5 * (y_top + y_bottom);
    let mut axial_force = 0.0;
    let mut crossing_bonds = 0usize;

    let nlocal = atoms.nlocal as usize;
    for i in 0..nlocal.min(bonds.bonds.len()) {
        let yi = atoms.pos[i][1] as f64;
        for b in &bonds.bonds[i] {
            if atoms.tag[i] > b.partner_tag {
                continue;
            }
            let Some(j) = tag_index(&atoms, b.partner_tag) else {
                continue;
            };
            let yj = atoms.pos[j][1] as f64;
            if (yi - mid_y) * (yj - mid_y) > 0.0 {
                continue;
            }
            let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
            let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
            let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(f64::MIN_POSITIVE);
            let r_i = dem.radius[i];
            let r_j = dem.radius[j];
            let r_b = bond_config.bond_radius_ratio * r_i.min(r_j);
            let area = PI * r_b * r_b;
            let k_n = match bond_config.youngs_modulus {
                Some(e) => e * area / b.r0,
                None => bond_config.normal_stiffness,
            };
            let f_n = k_n * (dist - b.r0);
            axial_force += (f_n * dy / dist).abs();
            crossing_bonds += 1;
        }
    }

    let stress = if crossing_bonds > 0 {
        axial_force / section_area
    } else {
        0.0
    };
    let live_now = live_bond_set(&atoms, &bonds);
    let broken_bonds = bond_metrics.total_bonds_broken;

    let newly_broken: Vec<(u32, u32)> = recorder
        .live_bonds_last
        .difference(&live_now)
        .copied()
        .collect();
    for (a, b) in newly_broken {
        let pa = recorder.initial_pos.get(&a).copied().unwrap_or([0.0; 3]);
        let pb = recorder.initial_pos.get(&b).copied().unwrap_or([0.0; 3]);
        if let Some(ref mut w) = recorder.cracks {
            writeln!(
                w,
                "{step},{strain:.9},{:.9},{:.9},{:.9},combined_stress",
                0.5 * (pa[0] + pb[0]),
                0.5 * (pa[1] + pb[1]),
                0.5 * (pa[2] + pb[2])
            )
            .ok();
        }
    }
    recorder.live_bonds_last = live_now;

    let live_bond_count = recorder.live_bonds_last.len();
    if let Some(ref mut w) = recorder.curve {
        writeln!(
            w,
            "{step},{strain:.9},{stress:.6},{:.9},{},{},{}",
            stress / 199.1e6,
            live_bond_count,
            broken_bonds,
            crossing_bonds
        )
        .ok();
    }
}
