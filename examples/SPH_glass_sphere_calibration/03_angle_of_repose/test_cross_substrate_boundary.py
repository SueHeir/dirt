#!/usr/bin/env python3
"""Regression guard for the standalone DEM scope boundary."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
CALIBRATION = Path(__file__).resolve().parents[1] / "calibration.yaml"
REPOSE_README = Path(__file__).with_name("README.md")
CARGO = REPO_ROOT / "Cargo.toml"


def main() -> None:
    readme = REPOSE_README.read_text(encoding="utf-8")
    normalized_readme = " ".join(readme.lower().split())
    cargo = CARGO.read_text(encoding="utf-8").lower()
    assert not CALIBRATION.exists(), (
        "a deleted SPH calibration hand-off must not be recreated by this DEM example"
    )
    assert "dev_soil_sph" not in cargo, "the DIRT manifest must not reintroduce the removed SPH crate"
    assert "standalone **dem** formation study" in normalized_readme
    assert "must not imply that `μ_r` is exported to an sph constitutive model" in normalized_readme
    assert "implementation diagnostic, not that validation" in normalized_readme
    print("PASS: rolling-friction study is explicitly standalone DEM with no SPH hand-off")


if __name__ == "__main__":
    main()
