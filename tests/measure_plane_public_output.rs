use std::fs;
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn run_measure_plane_case(name: &str, config_body: &str) -> Output {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock must be after Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("dirt-{name}-{}-{stamp}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temporary config directory");
    let config = dir.join("config.toml");
    fs::write(&config, config_body).expect("write temporary config");

    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
    let output = Command::new(cargo)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args([
            "run",
            "--quiet",
            "--no-default-features",
            "--features",
            "precision-double",
            "--example",
            "measure_plane_throughput",
            "--",
        ])
        .arg(&config)
        .output()
        .expect("run measure_plane_throughput example");

    fs::remove_dir_all(&dir).expect("remove temporary config directory");
    output
}

#[test]
fn staged_variable_dt_rates_are_current_in_thermo_output() {
    let output = run_measure_plane_case(
        "measure-plane-variable-dt-thermo",
        r#"
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.10
x_high = 0.10
y_low = -0.02
y_high = 0.02
z_low = -0.02
z_high = 0.02
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.2
bin_size = 0.012

[thermo]
columns = ["step", "crossings_gate", "flow_rate_gate", "cross_rate_gate"]

[[dem.materials]]
name = "glass"
youngs_mod = 8.7e9
poisson_ratio = 0.3
restitution = 0.9
friction = 0.0

[[particles.insert]]
material = "glass"
count = 1
radius = 0.005
density = 2500.0
velocity_x = 10.0
region = { type = "block", min = [-0.0301, -0.0001, -0.0001], max = [-0.0299, 0.0001, 0.0001] }
seed = 7

[[measure_plane]]
name = "gate"
point = [0.0, 0.0, 0.0]
normal = [1.0, 0.0, 0.0]
report_interval = 3

[[run]]
name = "slow"
dt = 0.001
steps = 2
thermo = 3

[[run]]
name = "fast"
dt = 0.002
steps = 2
thermo = 3
"#,
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "measure plane staged run failed\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );

    let row = stdout
        .lines()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.first().copied() == Some("3") && fields.len() >= 4 {
                Some(fields)
            } else {
                None
            }
        })
        .last()
        .unwrap_or_else(|| panic!("missing step-3 thermo row in stdout:\n{stdout}"));

    let crossings: f64 = row[1].parse().expect("parse crossings_gate");
    let flow_rate: f64 = row[2].parse().expect("parse flow_rate_gate");
    let cross_rate: f64 = row[3].parse().expect("parse cross_rate_gate");

    let radius = 0.005_f64;
    let density = 2500.0_f64;
    let mass = density * (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
    let elapsed = 0.001 + 2.0 * 0.002;

    assert!((crossings - 1.0).abs() < 1e-12);
    assert!((cross_rate - (1.0 / elapsed)).abs() < 1e-3);
    assert!((flow_rate - (mass / elapsed)).abs() < 5e-8);
}
