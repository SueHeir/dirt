//! bench_nohistory_tangential — validates the history-free `linear_nohistory`
//! tangential contact model against the documented LAMMPS `pair_granular`
//! velocity-Coulomb law, and demonstrates that it is genuinely distinct from the
//! history-based Mindlin path.
//!
//! # Experiment
//!
//! Two identical glass spheres are held at a FIXED normal overlap `δ` (so the
//! normal force `F_n = k_n δ` is constant) and driven with a prescribed, reversing
//! relative TANGENTIAL velocity — a triangle "load → unload → reverse" path
//! `0 → +V → −V → +V → −V`. This is a loading path along which contact history
//! matters: the Mindlin spring accumulates displacement `ξ`, so at the same
//! instantaneous `v_t` the force differs between loading and unloading.
//!
//! At every step the pair contact force ([`contact_force_core`]) is evaluated for
//! BOTH tangential models on independent state, and the recorder writes
//! `(model, step, v_t, F_n, F_t, |ξ|)` to a CSV. `sweep.py` then checks:
//!
//! - **linear_nohistory** reproduces `F_t = -min(μ|F_n|, γ_t|v_t|)` (velocity-
//!   Coulomb, LAMMPS `tangential linear_nohistory`) with `|ξ| ≡ 0` and is
//!   path-independent (equal `v_t` ⇒ equal `F_t`).
//! - **history** (Mindlin) accumulates `|ξ| > 0`, is path-dependent, and retains a
//!   nonzero elastic force at the `v_t = 0` crossings — the distinguishing behavior.
//!
//! Because positions are held fixed and velocities are prescribed each step, the
//! comparison is exact and deterministic — no time integration, no flakiness.
//!
//! ```bash
//! cargo run --release --example bench_nohistory_tangential --no-default-features \
//!     --features precision-double -- examples/bench_nohistory_tangential/config.toml
//! ```

use dirt_core::dirt_atom::{DemAtom, MaterialTable};
use dirt_core::dirt_granular::contact::{contact_force_core, ForcePass};
use dirt_core::dirt_granular::tangential::ContactHistoryStore;
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
    overlap: f64,
    dt: f64,
    v_amp: f64,
    legs: usize,
    steps_per_leg: usize,
}

/// Build a two-sphere system at fixed overlap with the given tangential model.
/// Returns `(atoms, neighbor, registry)`. Sphere 0 at the origin, sphere 1 along
/// +x at centre distance `2R − δ` so the overlap is exactly `δ`.
fn build_system(mat: &MaterialCfg, sc: &ScenarioCfg, tangential_model: &str) -> (Atom, Neighbor, AtomDataRegistry, MaterialTable) {
    let mut mt = MaterialTable::new();
    mt.contact_model = "hertz".to_string();
    mt.tangential_model = tangential_model.to_string();
    // rolling/twisting friction = 0 ⇒ tangential force is the only in-plane force.
    mt.add_material(
        &mat.name,
        mat.youngs_mod,
        mat.poisson_ratio,
        mat.restitution,
        mat.friction,
        0.0, // rolling friction
        0.0, // cohesion energy
    );
    mt.build_pair_tables();

    let r = sc.radius;
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);

    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = sc.dt;

    for (tag, x) in [(0u32, 0.0f64), (1u32, 2.0 * r - sc.overlap)] {
        atom.push_test_atom(tag, [x, 0.0, 0.0], r, mass);
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

    // Half neighbour list (newton on): atom 0 → {1}.
    let mut nb = Neighbor::new();
    nb.newton = true;
    nb.neighbor_offsets = vec![0, 1, 1];
    nb.neighbor_indices = vec![1];

    let mut reg = AtomDataRegistry::new();
    reg.register(dem);
    reg.register(hist);

    (atom, nb, reg, mt)
}

/// Triangle wave of tangential velocity: waypoints `0, +V, −V, +V, −V, …`.
fn tangential_velocity(step: usize, sc: &ScenarioCfg) -> f64 {
    let leg = (step / sc.steps_per_leg).min(sc.legs - 1);
    let frac = (step % sc.steps_per_leg) as f64 / sc.steps_per_leg as f64;
    let waypoint = |k: usize| -> f64 {
        if k == 0 {
            0.0
        } else if k % 2 == 1 {
            sc.v_amp
        } else {
            -sc.v_amp
        }
    };
    let a = waypoint(leg);
    let b = waypoint(leg + 1);
    a + (b - a) * frac
}

fn run_model(
    mat: &MaterialCfg,
    sc: &ScenarioCfg,
    tangential_model: &str,
    out: &mut fs::File,
) {
    let (mut atom, nb, reg, mt) = build_system(mat, sc, tangential_model);
    let total_steps = sc.legs * sc.steps_per_leg;

    for step in 0..total_steps {
        let vt = tangential_velocity(step, sc);

        // Prescribe the kinematic state: fixed positions (constant overlap), zero
        // spin, all relative tangential motion on sphere 1 along +y. Zero the force
        // accumulator (contact_force_core ADDS into it).
        atom.vel[0] = [0.0, 0.0, 0.0];
        atom.vel[1] = [0.0, vt, 0.0];
        atom.force[0] = [0.0, 0.0, 0.0];
        atom.force[1] = [0.0, 0.0, 0.0];

        contact_force_core(&mut atom, &nb, &reg, &mt, None, ForcePass::All);

        // Normal force is along +x (contact normal), tangential along +y.
        let f_n = (atom.force[0][0] as f64).abs();
        let f_t = atom.force[0][1] as f64;

        // Stored tangential spring displacement magnitude for the 0–1 contact.
        let xi = {
            let hist = reg.expect::<ContactHistoryStore>("read history");
            let s = hist.contacts[0]
                .iter()
                .find(|(t, _, _)| *t == 1)
                .map(|(_, s, _)| *s)
                .unwrap_or([0.0; 7]);
            (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt()
        };

        writeln!(
            out,
            "{},{},{:.12e},{:.12e},{:.12e},{:.12e}",
            tangential_model, step, vt, f_n, f_t, xi
        )
        .unwrap();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let config_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| "examples/bench_nohistory_tangential/config.toml".to_string());
    let text = fs::read_to_string(&config_path)
        .unwrap_or_else(|e| panic!("cannot read {config_path}: {e}"));
    let cfg: Config = toml::from_str(&text)
        .unwrap_or_else(|e| panic!("cannot parse {config_path}: {e}"));

    let out_dir = "examples/bench_nohistory_tangential/data";
    fs::create_dir_all(out_dir).ok();
    let results = format!("{out_dir}/nohistory_tangential_results.csv");
    let mut f = fs::File::create(&results)
        .unwrap_or_else(|e| panic!("cannot create {results}: {e}"));
    writeln!(f, "model,step,vt,fn,ft,xi").unwrap();

    println!("=== history-free tangential model benchmark ===");
    println!("  material: {} (mu = {})", cfg.material.name, cfg.material.friction);
    println!(
        "  fixed overlap = {:.3e} m, V = {} m/s, {} legs x {} steps",
        cfg.scenario.overlap, cfg.scenario.v_amp, cfg.scenario.legs, cfg.scenario.steps_per_leg
    );

    run_model(&cfg.material, &cfg.scenario, "linear_nohistory", &mut f);
    run_model(&cfg.material, &cfg.scenario, "history", &mut f);

    println!("  results -> {results}");
}
