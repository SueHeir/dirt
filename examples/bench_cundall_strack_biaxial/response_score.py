#!/usr/bin/env python3
"""Score a qualified independent biaxial trajectory without choosing a window.

This is deliberately separate from candidate admission.  Admission establishes
that an external solver ran the same apparatus; scoring then requires a
predeclared, source-cited criterion for *each* response required by the goal.
In particular, it refuses interpolation, peak picking, and a tolerance that is
not bound to the candidate package.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


SERIES = ("stress_ratio", "volumetric_strain", "fabric_anisotropy")
REQUIRED_COLUMNS = ("state",) + SERIES


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not set(REQUIRED_COLUMNS) <= set(rows[0]):
        raise RuntimeError("response must contain registered state and every required series")
    parsed = []
    for row in rows:
        try:
            parsed.append({name: float(row[name]) for name in REQUIRED_COLUMNS})
        except (KeyError, ValueError) as exc:
            raise RuntimeError("response contains non-numeric required values") from exc
    if any(not all(math.isfinite(value) for value in row.values()) for row in parsed):
        raise RuntimeError("response contains non-finite required values")
    if any(b["state"] <= a["state"] for a, b in zip(parsed, parsed[1:])):
        raise RuntimeError("response states must be strictly increasing")
    return parsed


def score(package_path: str | Path):
    """Score two SHA-256-bound, state-identical response tables.

    The package is intentionally strict: the states must be identical rather
    than resampled, and every tolerance needs a non-empty external citation.
    This avoids silently changing the comparison window or fitting a response
    before a positive replication is claimed.
    """
    package_path = Path(package_path)
    package = json.loads(package_path.read_text())
    if set(package) != {"dirt_response", "reference_response", "criteria"}:
        raise RuntimeError("score package must bind dirt_response, reference_response, and criteria")
    root = package_path.parent
    artifacts = {}
    for name, item in package.items():
        if set(item) != {"path", "sha256"}:
            raise RuntimeError(f"score package {name} must contain path and sha256")
        path = root / item["path"]
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"score package {name} artifact hash mismatch")
        artifacts[name] = path

    criteria = json.loads(artifacts["criteria"].read_text())
    if set(criteria) != set(SERIES):
        raise RuntimeError("criteria must set every required response series exactly once")
    limits = {}
    for name in SERIES:
        entry = criteria[name]
        if set(entry) != {"max_abs_error", "reference"} or not str(entry["reference"]).strip():
            raise RuntimeError("each criterion needs a positive bound and external reference")
        try:
            limit = float(entry["max_abs_error"])
        except (TypeError, ValueError) as exc:
            raise RuntimeError("criterion bound must be numeric") from exc
        if not math.isfinite(limit) or limit <= 0.0:
            raise RuntimeError("criterion bound must be finite and positive")
        limits[name] = limit

    dirt, reference = _rows(artifacts["dirt_response"]), _rows(artifacts["reference_response"])
    if [row["state"] for row in dirt] != [row["state"] for row in reference]:
        raise RuntimeError("DIRT and reference states must be identical; interpolation is forbidden")
    errors = {name: max(abs(a[name] - b[name]) for a, b in zip(dirt, reference)) for name in SERIES}
    return {"passed": all(errors[name] <= limits[name] for name in SERIES), "errors": errors, "limits": limits}
