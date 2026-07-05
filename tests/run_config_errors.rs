use std::fs;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn empty_discrete_radius_without_region_errors_before_region_sampling() {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock must be after Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "dirt-empty-discrete-radius-{}-{stamp}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).expect("create temporary config directory");
    let config = dir.join("config.toml");

    fs::write(
        &config,
        r#"
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = 0.0
x_high = 0.04
y_low = 0.0
y_high = 0.04
z_low = 0.0
z_high = 0.06
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.2
bin_size = 0.012

[[dem.materials]]
name = "glass"
youngs_mod = 8.7e9
poisson_ratio = 0.3
restitution = 0.3
friction = 0.5

[[particles.insert]]
material = "glass"
count = 2
radius = { distribution = "discrete", values = [], weights = [] }
density = 2500.0
seed = 1

[[run]]
name = "settle"
dt = 1.0e-5
steps = 0
thermo = 1
"#,
    )
    .expect("write temporary config");

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
            "run",
            "--",
        ])
        .arg(&config)
        .output()
        .expect("run generic config example");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !output.status.success(),
        "malformed config should exit non-zero"
    );
    assert!(
        stderr.contains(
            "ERROR: invalid radius in [[particles.insert]]: discrete radius requires at least one value"
        ),
        "stderr should contain typed config error, got:\n{stderr}"
    );
    assert!(
        !stderr.contains("panicked at"),
        "malformed config must not panic, got:\n{stderr}"
    );

    fs::remove_dir_all(&dir).expect("remove temporary config directory");
}
