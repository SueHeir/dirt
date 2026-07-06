use std::fs;
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn run_config_case_with_extra(name: &str, particles_insert: &str, extra_config: &str) -> Output {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock must be after Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("dirt-{name}-{}-{stamp}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temporary config directory");
    let config = dir.join("config.toml");

    fs::write(
        &config,
        format!(
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
{particles_insert}
density = 2500.0
seed = 1

{extra_config}

[[run]]
name = "settle"
dt = 1.0e-5
steps = 0
thermo = 1
"#
        ),
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

    fs::remove_dir_all(&dir).expect("remove temporary config directory");
    output
}

fn run_config_case(name: &str, particles_insert: &str) -> Output {
    run_config_case_with_extra(name, particles_insert, "")
}

fn run_config_error_case(name: &str, particles_insert: &str) -> String {
    let output = run_config_case(name, particles_insert);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        !output.status.success(),
        "malformed config should exit non-zero"
    );
    assert!(
        !stderr.contains("panicked at"),
        "malformed config must not panic, got:\n{stderr}"
    );

    stderr
}

fn run_wall_config_error_case(name: &str, wall_config: &str) -> String {
    let output = run_config_case_with_extra(
        name,
        r#"
count = 1
radius = 0.001
"#,
        wall_config,
    );
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        !output.status.success(),
        "malformed wall config should exit non-zero"
    );
    assert!(
        !stderr.contains("panicked at"),
        "malformed wall config must not panic, got:\n{stderr}"
    );

    stderr
}

#[test]
fn empty_discrete_radius_without_region_errors_before_region_sampling() {
    let stderr = run_config_error_case(
        "empty-discrete-radius",
        r#"
count = 2
radius = { distribution = "discrete", values = [], weights = [] }
"#,
    );
    assert!(
        stderr.contains(
            "ERROR: invalid radius in [[particles.insert]]: discrete radius requires at least one value"
        ),
        "stderr should contain typed config error, got:\n{stderr}"
    );
}

#[test]
fn infinite_immediate_insert_velocity_errors_before_normal_distribution() {
    let stderr = run_config_error_case(
        "infinite-immediate-velocity",
        r#"
count = 1
radius = 0.001
velocity = inf
"#,
    );
    assert!(
        stderr.contains(
            "ERROR: velocity in [[particles.insert]] must be finite and non-negative, got inf"
        ),
        "stderr should contain typed config error, got:\n{stderr}"
    );
}

#[test]
fn infinite_rate_insert_velocity_errors_before_normal_distribution() {
    let stderr = run_config_error_case(
        "infinite-rate-velocity",
        r#"
rate = 1
rate_interval = 1
radius = 0.001
velocity = inf
"#,
    );
    assert!(
        stderr.contains(
            "ERROR: velocity in rate-based [[particles.insert]] must be finite and non-negative, got inf"
        ),
        "stderr should contain typed config error, got:\n{stderr}"
    );
}

#[test]
fn too_large_immediate_radius_errors_before_region_sampling() {
    let stderr = run_config_error_case(
        "too-large-immediate-radius",
        r#"
count = 1
radius = 1.0
"#,
    );
    assert!(
        stderr.contains(
            "ERROR: [[particles.insert]] default insertion region is smaller than particle: radius 1 exceeds domain extent"
        ),
        "stderr should contain typed config error, got:\n{stderr}"
    );
    assert!(
        !stderr.contains("cannot sample empty range"),
        "bad region should be rejected before sampling, got:\n{stderr}"
    );
}

#[test]
fn too_large_rate_radius_errors_before_region_sampling() {
    let stderr = run_config_error_case(
        "too-large-rate-radius",
        r#"
rate = 1
rate_interval = 1
radius = 1.0
"#,
    );
    assert!(
        stderr.contains(
            "ERROR: rate-based [[particles.insert]] default insertion region is smaller than particle: radius 1 exceeds domain extent"
        ),
        "stderr should contain typed config error, got:\n{stderr}"
    );
    assert!(
        !stderr.contains("cannot sample empty range"),
        "bad region should be rejected before sampling, got:\n{stderr}"
    );
}

#[test]
fn valid_default_region_still_runs_immediate_insert() {
    let output = run_config_case(
        "valid-default-region",
        r#"
count = 1
radius = 0.001
"#,
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success(),
        "valid insertion config should run, stderr:\n{stderr}"
    );
    assert!(
        stdout.contains("DemAtomInsert: inserting 1 particles"),
        "stdout should show the valid insert path ran, got:\n{stdout}"
    );
}

#[test]
fn cylinder_wall_missing_center_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "cylinder-wall-missing-center",
        r#"
[[wall]]
type = "cylinder"
axis = "z"
radius = 0.01
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: cylinder wall requires 'center' [c0, c1]"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn cylinder_wall_missing_radius_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "cylinder-wall-missing-radius",
        r#"
[[wall]]
type = "cylinder"
axis = "z"
center = [0.02, 0.02]
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: cylinder wall requires 'radius'"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn cylinder_wall_wrong_center_length_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "cylinder-wall-wrong-center-length",
        r#"
[[wall]]
type = "cylinder"
axis = "z"
center = [0.02, 0.02, 0.02]
radius = 0.01
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: cylinder wall 'center' must have 2 elements"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn cylinder_wall_bad_axis_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "cylinder-wall-bad-axis",
        r#"
[[wall]]
type = "cylinder"
axis = "q"
center = [0.02, 0.02]
radius = 0.01
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: cylinder wall axis must be x, y, or z, got 'q'"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn sphere_wall_missing_center_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "sphere-wall-missing-center",
        r#"
[[wall]]
type = "sphere"
radius = 0.01
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: sphere wall requires 'center' [x, y, z]"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn sphere_wall_missing_radius_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "sphere-wall-missing-radius",
        r#"
[[wall]]
type = "sphere"
center = [0.02, 0.02, 0.02]
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: sphere wall requires 'radius'"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn sphere_wall_wrong_center_length_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "sphere-wall-wrong-center-length",
        r#"
[[wall]]
type = "sphere"
center = [0.02, 0.02]
radius = 0.01
material = "glass"
"#,
    );
    assert!(
        stderr.contains("ERROR: sphere wall 'center' must have 3 elements"),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}

#[test]
fn unknown_wall_type_errors_without_panic() {
    let stderr = run_wall_config_error_case(
        "unknown-wall-type",
        r#"
[[wall]]
type = "capsule"
normal_z = 1.0
material = "glass"
"#,
    );
    assert!(
        stderr.contains(
            "ERROR: unknown wall type in [[wall]]: 'capsule'. Expected 'plane', 'cylinder', 'sphere', or 'region'"
        ),
        "stderr should contain typed wall config error, got:\n{stderr}"
    );
}
