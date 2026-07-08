//! bench_mindlin_rescale_tangential — isolates the Mindlin unloading-rescale
//! tangential variants against the documented LAMMPS recurrence.
//!
//! Two identical spheres are held at prescribed overlaps. The contact is first
//! tangentially loaded at fixed peak overlap, then normally unloaded with zero
//! tangential velocity. With no damping (`restitution = 1`) and a high Coulomb
//! cap, the force law reduces to a direct history check:
//!
//! - `history`: displacement history stays constant on normal unloading, so
//!   `F_t ∝ a`.
//! - `mindlin_rescale`: displacement history is gated by `a/a_prev`, so
//!   `F_t ∝ a^2` during unloading.
//! - `mindlin_rescale/force`: elastic-force history is gated by `a/a_prev`.
//! - `linear_nohistory`: no history survives, so unloading force is zero.

use dirt_core::dirt_atom::{DemAtom, MaterialTable};
use dirt_core::dirt_granular::contact::{contact_force_core, ForcePass};
use dirt_core::dirt_granular::tangential::{ContactHistoryStore, CONTACT_HISTORY_LEN};
use dirt_core::soil_core::{Atom, AtomDataRegistry, Neighbor};
use serde::Deserialize;
use std::fs;
use std::io::Write as IoWrite;

#[derive(Deserialize)]
struct Config {
    material: MaterialCfg,
    scenario: ScenarioCfg,
}

#[derive(Deserialize)]
struct MaterialCfg {
    name: String,
    youngs_mod: f64,
    poisson_ratio: f64,
    restitution: f64,
    friction: f64,
}

#[derive(Deserialize)]
struct ScenarioCfg {
    radius: f64,
    peak_overlap: f64,
    final_overlap: f64,
    dt: f64,
    load_steps: usize,
    unload_steps: usize,
    tangential_velocity: f64,
}

fn build_system(
    mat: &MaterialCfg,
    sc: &ScenarioCfg,
    tangential_model: &str,
) -> (Atom, Neighbor, AtomDataRegistry, MaterialTable) {
    let mut mt = MaterialTable::new();
    mt.contact_model = "hertz".to_string();
    mt.tangential_model = tangential_model.to_string();
    mt.add_material(
        &mat.name,
        mat.youngs_mod,
        mat.poisson_ratio,
        mat.restitution,
        mat.friction,
        0.0,
        0.0,
    );
    mt.build_pair_tables();

    let r = sc.radius;
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);

    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = sc.dt;

    for tag in [0u32, 1u32] {
        atom.push_test_atom(tag, [0.0, 0.0, 0.0], r, mass);
        dem.radius.push(r);
        dem.density.push(density);
        dem.inv_inertia.push(1.0 / (0.4 * mass * r * r));
        dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
        dem.omega.push([0.0; 3]);
        dem.ang_mom.push([0.0; 3]);
        dem.torque.push([0.0; 3]);
        dem.body_id.push(0.0);
        hist.contacts.push(Vec::new());
    }
    atom.nlocal = 2;
    atom.natoms = 2;

    let mut nb = Neighbor::new();
    nb.newton = true;
    nb.neighbor_offsets = vec![0, 1, 1];
    nb.neighbor_indices = vec![1];

    let mut reg = AtomDataRegistry::new();
    reg.register(dem);
    reg.register(hist);
    (atom, nb, reg, mt)
}

fn run_model(mat: &MaterialCfg, sc: &ScenarioCfg, model: &str, out: &mut fs::File) {
    let (mut atom, nb, reg, mt) = build_system(mat, sc, model);
    let r = sc.radius;
    let total_steps = sc.load_steps + sc.unload_steps;

    for step in 0..total_steps {
        let (phase, overlap, vt) = if step < sc.load_steps {
            ("load", sc.peak_overlap, sc.tangential_velocity)
        } else {
            let k = step - sc.load_steps;
            let frac = if sc.unload_steps > 1 {
                k as f64 / (sc.unload_steps - 1) as f64
            } else {
                1.0
            };
            let overlap = sc.peak_overlap + frac * (sc.final_overlap - sc.peak_overlap);
            ("unload", overlap, 0.0)
        };

        atom.pos[0] = [0.0, 0.0, 0.0];
        atom.pos[1] = [(2.0 * r - overlap) as dirt_core::soil_core::Real, 0.0, 0.0];
        atom.vel[0] = [0.0, 0.0, 0.0];
        atom.vel[1] = [0.0, vt as dirt_core::soil_core::Real, 0.0];
        atom.force[0] = [0.0, 0.0, 0.0];
        atom.force[1] = [0.0, 0.0, 0.0];

        contact_force_core(&mut atom, &nb, &reg, &mt, None, ForcePass::All);

        let hist = reg.expect::<ContactHistoryStore>("read history");
        let h = hist.contacts[0]
            .iter()
            .find(|(t, _, _)| *t == 1)
            .map(|(_, s, _)| *s)
            .unwrap_or([0.0; CONTACT_HISTORY_LEN]);
        let hmag = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]).sqrt();
        let contact_radius = h[7];

        writeln!(
            out,
            "{model},{phase},{step},{overlap:.12e},{contact_radius:.12e},{vt:.12e},{:.12e},{:.12e},{hmag:.12e}",
            (atom.force[0][0] as f64).abs(),
            atom.force[0][1] as f64,
        )
        .unwrap();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let config_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "examples/bench_mindlin_rescale_tangential/config.toml".to_string());
    let text = fs::read_to_string(&config_path)
        .unwrap_or_else(|e| panic!("cannot read {config_path}: {e}"));
    let cfg: Config =
        toml::from_str(&text).unwrap_or_else(|e| panic!("cannot parse {config_path}: {e}"));

    let out_dir = "examples/bench_mindlin_rescale_tangential/data";
    fs::create_dir_all(out_dir).ok();
    let results = format!("{out_dir}/mindlin_rescale_tangential_results.csv");
    let mut f =
        fs::File::create(&results).unwrap_or_else(|e| panic!("cannot create {results}: {e}"));
    writeln!(f, "model,phase,step,overlap,a,vt,fn,ft,hmag").unwrap();

    for model in [
        "history",
        "mindlin_rescale",
        "mindlin_rescale/force",
        "linear_nohistory",
    ] {
        run_model(&cfg.material, &cfg.scenario, model, &mut f);
    }

    println!("=== Mindlin unloading-rescale tangential benchmark ===");
    println!("  results -> {results}");
}
