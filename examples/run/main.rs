//! run — the generic, config-only DIRT driver.
//!
//! This is a *single* example binary that assembles the standard DIRT plugin
//! stack and then hands everything — geometry, materials, particle insertion,
//! walls, body forces, run stages, and output — to a TOML config supplied on
//! the command line. There is **no per-case Rust**: every case is a `config.toml`.
//!
//! ```bash
//! cargo run --release --example run -- examples/run/config.toml
//! ```
//!
//! Shipped configs (each a complete, self-contained scenario — pick one on argv):
//!
//! | Config | Scenario | Geometry | Insertion |
//! |--------|----------|----------|-----------|
//! | `config.toml`        | two-sphere settle       | floor plane           | 2 fixed spheres |
//! | `pour_settle.toml`   | granular pour & settle  | box (5 plane walls)   | polydisperse cloud |
//! | `pour_cylinder.toml` | bidisperse silo pour    | cylinder wall + floor | two materials |
//! | `shear_box.toml`     | Lees–Edwards simple shear | triperiodic cell    | dense glass pack |
//! | `compression_box.toml` | uniaxial compression    | triperiodic cell    | loose glass pack |
//!
//! The plugin set below is deliberately a superset of what any one simple case
//! needs. Each plugin is a no-op when its config section is absent (gravity
//! defaults to zero body force only if you set it; walls/fixes register nothing
//! when no `[[wall]]` / `[[*]]` fixes are declared), so the same driver runs a
//! two-particle rebound, a settle-into-a-box, and a granular pour — chosen
//! entirely by the config.
//!
//! Standard dump output is produced via the config's `[dump]` (per-atom CSV /
//! binary snapshots) and/or `[vtp]` sections, handled by the core `PrintPlugin`.
//!
//! # Declarative loading (`[loading]`) — deformation without a schedule
//!
//! A shear or compression box is a *loaded* simulation: the cell is deformed at
//! a controlled rate up to a target strain. We express that **declaratively** —
//! the config says *what* loading to apply, never *how* to sequence it:
//!
//! ```toml
//! [loading]
//! type         = "shear"     # "shear" (Lees–Edwards xy) | "compression" (uniaxial)
//! rate         = 50.0        # engineering strain rate (1/s), magnitude
//! target_strain = 2.0        # dimensionless total strain to reach
//! dt           = 2.0e-7      # timestep (s)
//! # axis       = "z"         # compression axis (compression only; default z)
//! # thermo     = 10000       # console/dump cadence (default: target reached / 20)
//! ```
//!
//! The **runner** — not the user — interprets this spec: it derives the step
//! count from `target_strain / (|rate| · dt)`, wires the underlying box-deform
//! driver (`DeformPlugin`), and owns the single step loop that carries the cell
//! to the target strain. There is **no user-authored stage schedule** and no
//! imperative sequence of operations: `[loading]` is a *description* of the
//! deformation, expanded by [`expand_loading`] into the driver's internal run.
//! Authoring both a `[loading]` block and an explicit `[[run]]`/`[deform]`
//! section is rejected — the point is that loading is declarative, not scripted.

use std::env;
use std::process;

use dirt_core::prelude::*;
use dirt_core::soil_core::toml;

fn main() {
    // Declarative-loading pre-pass: if the config carries a `[loading]` block,
    // expand it (in the runner, not the config) into the internal `[deform]` +
    // single `[[run]]` the driver executes, then seed the `Config`/`Input`
    // resources so the standard plugins pick them up unchanged. Configs without
    // `[loading]` fall through untouched and are loaded by the InputPlugin.
    let mut app = App::new();
    seed_from_loading(&mut app);

    app
        // Infrastructure: CLI parse + TOML load, comm, domain, neighbor, run, print.
        .add_plugins(CorePlugins)
        // DEM: atom data + insertion, Velocity Verlet, Hertz-Mindlin contact, rotation.
        .add_plugins(GranularDefaultPlugins)
        // Optional, config-driven physics — each inert unless its section is present.
        .add_plugins(GravityPlugin) // [gravity] body force
        .add_plugins(WallPlugin) // [[wall]] container / boundary faces
        .add_plugins(FixesPlugin) // [[addforce]] / [[freeze]] / [[viscous]] / ...
        .add_plugins(DeformPlugin); // [deform] / [loading] box deformation driver

    // Everything case-specific comes from the TOML config path on argv.
    app.start();
}

/// If argv points at a config with a `[loading]` section, load it, expand the
/// declarative loading spec into the driver's internal run, and seed the
/// `Config`/`Input` resources (which makes the `InputPlugin` a no-op). Configs
/// without `[loading]` are left entirely to the normal `InputPlugin` path.
fn seed_from_loading(app: &mut App) {
    let args: Vec<String> = env::args().collect();
    // Let InputPlugin handle --generate-config and the usage/error paths.
    let Some(path) = args.get(1) else { return };
    if path.starts_with("--") {
        return;
    }

    let mut table = load_toml(path);
    if !table.contains_key("loading") {
        return; // no declarative loading — normal config, hand off to InputPlugin.
    }

    expand_loading(&mut table);

    // Mirror InputPlugin's output-dir resolution: prefer [output].dir, else the
    // config file's parent directory.
    let output_dir = table
        .get("output")
        .and_then(|v| v.as_table())
        .and_then(|t| t.get("dir"))
        .and_then(|v| v.as_str())
        .map(|s| s.to_string())
        .or_else(|| {
            std::path::Path::new(path)
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .map(|p| p.to_string_lossy().into_owned())
        });

    app.add_resource(Input {
        filename: path.clone(),
        output_dir,
    });
    app.add_resource(Config { table });
}

/// Expand a declarative `[loading]` block into the driver's internal `[deform]`
/// + single `[[run]]` stage, in place.
///
/// This is where the runner *interprets* the loading description and takes
/// ownership of the step loop. The user supplies only *what* to load
/// (deformation type, strain rate, target strain, timestep); the runner derives
/// the number of steps and the box-deform wiring. It is deliberately not a
/// scripting hook: exactly one loading spec, no operation sequencing.
///
/// Supported `type`s:
/// - `"shear"`       — Lees–Edwards simple shear, `xy` at engineering
///   shear-strain rate `rate` (1/s). Requires periodic x and y.
/// - `"compression"` — uniaxial engineering-strain compression on `axis`
///   (default `z`); the box shrinks at magnitude `rate`.
///
/// Stopping is strain-controlled: `steps = ceil(target_strain / (|rate|·dt))`.
fn expand_loading(table: &mut toml::Table) {
    // Reject mixing declarative loading with a hand-authored schedule/deform.
    for forbidden in ["run", "deform"] {
        if table.contains_key(forbidden) {
            die(&format!(
                "[loading] is declarative and owns the run; do not also specify \
                 `[{forbidden}]` — remove one. Use `[loading]` alone (the runner \
                 derives the steps and deformation), or drop `[loading]` and drive \
                 the box yourself with `[[run]]`/`[deform]`."
            ));
        }
    }

    let loading = table
        .get("loading")
        .and_then(|v| v.as_table())
        .unwrap_or_else(|| die("`[loading]` must be a table"))
        .clone();

    let ltype = req_str(&loading, "type");
    let rate = req_f64(&loading, "rate");
    let target_strain = req_f64(&loading, "target_strain");
    let dt = req_f64(&loading, "dt");

    if rate == 0.0 {
        die("[loading] `rate` must be non-zero");
    }
    if dt <= 0.0 {
        die("[loading] `dt` must be > 0");
    }
    if target_strain <= 0.0 {
        die("[loading] `target_strain` must be > 0 (magnitude of total strain)");
    }

    // Runner owns the loop length: steps to reach the target engineering strain.
    // Round to the nearest integer when the exact quotient lands on one (within a
    // small relative tolerance — floating-point makes e.g. 1.0/(100·2e-7) evaluate
    // to 50000.0000000001); otherwise ceil so we never undershoot the target.
    let raw = target_strain / (rate.abs() * dt);
    let nearest = raw.round();
    let steps = if (raw - nearest).abs() <= 1e-6 * nearest.max(1.0) {
        nearest as i64
    } else {
        raw.ceil() as i64
    };
    if steps <= 0 {
        die("[loading] derived a non-positive step count — check rate/target_strain/dt");
    }

    let thermo = opt_i64(&loading, "thermo").unwrap_or((steps / 20).max(1));
    let name = opt_str(&loading, "type_name")
        .or_else(|| opt_str(&loading, "name"))
        .unwrap_or_else(|| ltype.clone());

    // Build the internal [deform] table for the requested loading type.
    let mut deform = toml::Table::new();
    match ltype.as_str() {
        "shear" => {
            // Lees–Edwards simple shear: xy tilt at engineering shear-strain rate.
            deform.insert("xy".to_string(), erate_axis(rate));
        }
        "compression" => {
            // Uniaxial compression on `axis`: box shrinks → negative engineering
            // strain rate regardless of the sign the user wrote for `rate`.
            let axis = opt_str(&loading, "axis").unwrap_or_else(|| "z".to_string());
            if !matches!(axis.as_str(), "x" | "y" | "z") {
                die("[loading] `axis` must be one of \"x\", \"y\", \"z\"");
            }
            deform.insert(axis, erate_axis(-rate.abs()));
        }
        other => die(&format!(
            "[loading] unknown `type` = {other:?}; expected \"shear\" or \"compression\""
        )),
    }
    table.insert("deform".to_string(), toml::Value::Table(deform));

    // Build the single internal run stage the runner drives to the target strain.
    let mut stage = toml::Table::new();
    stage.insert("name".to_string(), toml::Value::String(name));
    stage.insert("steps".to_string(), toml::Value::Integer(steps));
    stage.insert("dt".to_string(), toml::Value::Float(dt));
    stage.insert("thermo".to_string(), toml::Value::Integer(thermo));
    table.insert(
        "run".to_string(),
        toml::Value::Array(vec![toml::Value::Table(stage)]),
    );
}

/// Build an `{ style = "erate", rate = <rate> }` axis-deform inline table.
fn erate_axis(rate: f64) -> toml::Value {
    let mut t = toml::Table::new();
    t.insert(
        "style".to_string(),
        toml::Value::String("erate".to_string()),
    );
    t.insert("rate".to_string(), toml::Value::Float(rate));
    toml::Value::Table(t)
}

// ── small TOML accessors (accept ints where floats are expected) ─────────────

fn req_str(t: &toml::Table, key: &str) -> String {
    opt_str(t, key).unwrap_or_else(|| die(&format!("[loading] missing required string `{key}`")))
}

fn opt_str(t: &toml::Table, key: &str) -> Option<String> {
    t.get(key).and_then(|v| v.as_str()).map(|s| s.to_string())
}

fn req_f64(t: &toml::Table, key: &str) -> f64 {
    match t.get(key) {
        Some(toml::Value::Float(f)) => *f,
        Some(toml::Value::Integer(i)) => *i as f64,
        Some(_) => die(&format!("[loading] `{key}` must be a number")),
        None => die(&format!("[loading] missing required number `{key}`")),
    }
}

fn opt_i64(t: &toml::Table, key: &str) -> Option<i64> {
    match t.get(key) {
        Some(toml::Value::Integer(i)) => Some(*i),
        Some(toml::Value::Float(f)) => Some(*f as i64),
        Some(_) => die(&format!("[loading] `{key}` must be an integer")),
        None => None,
    }
}

/// Print a friendly error and exit (matches the driver's other user-facing
/// failure paths — a bad config is a user error, not a panic/backtrace).
fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    process::exit(1);
}
