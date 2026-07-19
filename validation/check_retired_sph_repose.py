#!/usr/bin/env python3
"""Check the repository boundary left by DIRT's SPH-suite removal.

This is deliberately a source-tree audit, not a physics benchmark.  A PASS
only establishes that the retired DIRT command and its handoff artefacts are
not silently present in this revision.  It cannot establish an angle of
repose, a glass reference, or a transferable rolling-friction coefficient.

Run from any directory with:

    python3 validation/check_retired_sph_repose.py
"""

from __future__ import annotations

import subprocess
import sys
import re
from pathlib import Path


RETIRED_PREFIXES = (
    "examples/SPH_glass_sphere_calibration/",
    "crates/dirt_sph/",
)
RETIRED_TOKENS = (
    "dev_soil_sph",
    "sphcal_angle_of_repose",
    "SPH_glass_sphere_calibration",
)
REQUIRED_DISCLOSURES = (
    "no runnable surface",
    "not a dirt validation claim",
    "not an executable pass gate",
    "not a calibration pass",
)


def repo_root() -> Path:
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode:
        raise RuntimeError("must run inside a DIRT git worktree")
    return Path(result.stdout.strip())


def tracked_files(root: Path) -> list[str]:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip())
    return result.stdout.splitlines()


def main() -> int:
    root = repo_root()
    files = tracked_files(root)
    errors: list[str] = []

    retired_paths = [
        path for path in files if path.startswith(RETIRED_PREFIXES)
    ]
    if retired_paths:
        errors.append("retired SPH source paths are tracked: " + ", ".join(retired_paths))

    cargo_sources = ["Cargo.toml", "Cargo.lock"]
    for source in cargo_sources:
        path = root / source
        if path.exists():
            text = path.read_text(encoding="utf-8")
            for token in RETIRED_TOKENS:
                if token in text:
                    errors.append(f"{source} still refers to retired token {token!r}")

    note = root / "docs/src/retired/sph-glass-angle-of-repose.md"
    if not note.is_file():
        errors.append(f"missing retirement note: {note.relative_to(root)}")
    else:
        note_text = re.sub(r"\s+", " ", note.read_text(encoding="utf-8").lower())
        for disclosure in REQUIRED_DISCLOSURES:
            if disclosure not in note_text:
                errors.append(f"retirement note lacks disclosure: {disclosure!r}")

    if errors:
        print("RETIREMENT BOUNDARY: FAIL")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("RETIREMENT BOUNDARY: PASS")
    print("  DIRT contains no retired SPH repose executable or SPH dependency.")
    print("  This source audit is not a calibration result and does not satisfy its gate.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"RETIREMENT BOUNDARY: FAIL\n  - {error}")
        raise SystemExit(1)
