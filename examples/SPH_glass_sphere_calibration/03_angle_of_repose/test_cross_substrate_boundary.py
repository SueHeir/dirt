#!/usr/bin/env python3
"""Regression guard for the DIRT-to-SPH rolling-friction hand-off boundary."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CALIBRATION = ROOT / "calibration.yaml"
REPOSE_README = Path(__file__).with_name("README.md")


def main() -> None:
    calibration = CALIBRATION.read_text(encoding="utf-8")
    readme = REPOSE_README.read_text(encoding="utf-8")
    assert "mu_r_pinned: null" in calibration, (
        "an unqualified numerical rolling-friction value must not cross the DIRT-to-SPH boundary"
    )
    assert "mu_r_status: \"withheld:" in calibration
    # A null field is not sufficient if a conflicting prose header still labels
    # a numerical rolling value as canonical.  That value can be copied into a
    # downstream SPH/DEM setup even though no transferable repose closure exists.
    assert "μ_r=0.10" not in calibration
    assert "Rolling friction is deliberately omitted" in calibration
    assert "does **not** currently produce an SPH constitutive input" in readme
    assert "cross-substrate validation" in readme
    print("PASS: unqualified DIRT rolling friction is withheld from the SPH hand-off")


if __name__ == "__main__":
    main()
