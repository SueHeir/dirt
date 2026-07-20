#!/usr/bin/env python3
"""Artifact-bound admission for an externally produced biaxial trajectory.

The protocol ledger answers whether two apparatuses are comparable.  That is
necessary but not sufficient: a future all-``yes`` ledger must not become a
reference merely by assertion.  This module binds every proposed positive
candidate to immutable input/output artifacts and verifies that its response
table actually contains each response class required by the replication goal.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

from candidate_admission import REQUIRED_PROTOCOL_FIELDS, decide

REQUIRED_COLUMNS = {
    "state_registration": "state",
    "stress_or_deviatoric_path": "stress_ratio",
    "volumetric_strain_or_dilatancy_path": "volumetric_strain",
    "contact_or_fabric_evolution": "fabric_anisotropy",
}


def sha256(path: Path) -> str:
    """Return the SHA-256 identity of one candidate artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def admit(package_path: str | Path):
    """Return the candidate decision only after artifacts and series are checked.

    The JSON package names a ledger, an immutable solver input, and a response
    CSV.  Relative paths are resolved beside the package so a manifest cannot
    silently depend on the reviewer's current directory.
    """
    package_path = Path(package_path)
    package = json.loads(package_path.read_text())
    required = {"candidate", "ledger", "input", "response"}
    if set(package) != required:
        raise RuntimeError("candidate package must contain candidate, ledger, input, and response")
    root = package_path.parent
    artifacts = {}
    for name in ("ledger", "input", "response"):
        item = package[name]
        if set(item) != {"path", "sha256"}:
            raise RuntimeError(f"candidate package {name} must contain path and sha256")
        path = root / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise RuntimeError(f"candidate package {name} artifact hash mismatch")
        artifacts[name] = path

    decision = decide(artifacts["ledger"])
    if not decision.eligible:
        return decision

    with artifacts["response"].open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not set(REQUIRED_COLUMNS.values()) <= set(rows[0]):
        raise RuntimeError("candidate response lacks required registered response columns")
    for row in rows:
        try:
            values = [float(row[column]) for column in REQUIRED_COLUMNS.values()]
        except (KeyError, ValueError) as exc:
            raise RuntimeError("candidate response has non-numeric required values") from exc
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("candidate response has non-finite required values")
    states = [float(row[REQUIRED_COLUMNS["state_registration"]]) for row in rows]
    if any(b <= a for a, b in zip(states, states[1:])):
        raise RuntimeError("candidate response state must be strictly increasing")
    if set(decision.failures) & set(REQUIRED_PROTOCOL_FIELDS):
        raise RuntimeError("inconsistent eligible candidate decision")
    return decision
