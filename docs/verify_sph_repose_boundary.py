#!/usr/bin/env python3
"""Audit whether the retired DIRT repose claim has a maintained SPH target.

This is deliberately a *withholding* audit.  A zero exit status establishes
only that the specified SPH revision has no compatible executable/parameter
surface; it is never a calibration pass and must not be used as one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import unittest
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

RETIRED = "examples/SPH_glass_sphere_calibration/03_angle_of_repose"
DOI = "10.1016/S0378-4371(99)00183-1"
# A generic "glass" example is not a repose experiment.  Keeping it out of
# this predicate prevents the audit from treating unrelated material examples
# as a reason to infer a calibration surface.  Conversely, an actual repose
# or calibration path must stop this withholding audit and trigger a new,
# model-specific review.
TARGET_RE = re.compile(r"(?:angle[_ -]?of[_ -]?repose|repose|calibrat)", re.IGNORECASE)
ROLLING_FIELD_RE = re.compile(
    r"^\s*pub\s+(?:rolling_friction|rolling_resistance|mu_r)\s*:", re.MULTILINE
)


def git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout


def material_block(source: str) -> str:
    match = re.search(r"pub struct MaterialParams\s*\{(?P<body>.*?)\n\}", source, re.S)
    if not match:
        raise ValueError("MaterialParams struct was not found")
    return match.group("body")


def inspect(paths: list[str], params_source: str) -> dict[str, object]:
    """Classify the solver interface without looking at any DIRT result."""
    candidates = [p for p in paths if p.startswith("examples/") and TARGET_RE.search(p)]
    params = material_block(params_source)
    return {
        "candidate_examples": candidates,
        # Only public fields of the maintained constitutive interface count.
        # Documentation mentions, comments, and unrelated names such as
        # `mu_ref` must not be promoted into a rolling-contact interface.
        "rolling_parameter_present": bool(ROLLING_FIELD_RE.search(params)),
        "material_fields": re.findall(r"pub\s+(\w+)\s*:", params),
        "material_source_sha256": hashlib.sha256(params_source.encode()).hexdigest(),
    }


def crossref_identity() -> dict[str, object]:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(DOI, safe="")
    with urllib.request.urlopen(url, timeout=20) as response:  # nosec B310: fixed HTTPS host
        message = json.load(response)["message"]
    title = message.get("title", [""])[0]
    venue = message.get("container-title", [""])[0]
    year = message.get("published", {}).get("date-parts", [[None]])[0][0]
    if (message.get("DOI", "").lower() != DOI.lower()
            or title != "Rolling friction in the dynamic simulation of sandpile formation"
            or not venue.startswith("Physica A") or year != 1999):
        raise ValueError("Crossref record does not match the declared negative-control citation")
    return {"doi": DOI, "title": title, "venue": venue, "year": year,
            "disposition": "identity_only_protocol_incompatible"}


def audit(dirt_repo: Path, sph_repo: Path, rev: str) -> dict[str, object]:
    dirt_paths = git(dirt_repo, "ls-tree", "-r", "--name-only", "HEAD").splitlines()
    sph_paths = git(sph_repo, "ls-tree", "-r", "--name-only", rev).splitlines()
    source = git(sph_repo, "show", f"{rev}:crates/sph_constitutive/src/lib.rs")
    surface = inspect(sph_paths, source)
    report: dict[str, object] = {
        "audit_kind": "cross_substrate_precondition",
        "dirt_revision": git(dirt_repo, "rev-parse", "HEAD").strip(),
        "dirt_retired_path_present": RETIRED in dirt_paths,
        "sph_revision": git(sph_repo, "rev-parse", rev).strip(),
        "sph_surface": surface,
        "external_reference": crossref_identity(),
    }
    if report["dirt_retired_path_present"]:
        raise ValueError("DIRT still contains the retired SPH executable; retirement claim is stale")
    if surface["candidate_examples"] or surface["rolling_parameter_present"]:
        raise ValueError("SPH has changed: re-evaluate the calibration with a maintained protocol")
    report["status"] = "WITHHELD_NO_MAINTAINED_TARGET"
    report["not_a_calibration_pass"] = True
    return report


class BoundaryAuditTests(unittest.TestCase):
    def test_continuum_surface_is_not_a_rolling_interface(self) -> None:
        result = inspect(["examples/simple_shear/main.rs"], "pub struct MaterialParams {\n pub mu_s: f64,\n pub d: f64,\n}")
        self.assertEqual(result["candidate_examples"], [])
        self.assertFalse(result["rolling_parameter_present"])

    def test_glass_name_alone_is_not_a_repose_surface(self) -> None:
        result = inspect(["examples/glass_pour/main.rs"], "pub struct MaterialParams {\n pub mu_s: f64,\n}")
        self.assertEqual(result["candidate_examples"], [])

    def test_rolling_parameter_is_detected(self) -> None:
        result = inspect([], "pub struct MaterialParams {\n pub rolling_friction: f64,\n}")
        self.assertTrue(result["rolling_parameter_present"])

    def test_repose_executable_is_detected(self) -> None:
        result = inspect(["examples/glass_angle_of_repose/main.rs"], "pub struct MaterialParams {\n pub mu_s: f64,\n}")
        self.assertEqual(result["candidate_examples"], ["examples/glass_angle_of_repose/main.rs"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dirt-repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--sph-repo", type=Path)
    parser.add_argument("--rev", default="origin/main")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return 0 if unittest.main(argv=[sys.argv[0]], exit=False).result.wasSuccessful() else 1
    if args.sph_repo is None:
        parser.error("--sph-repo is required unless --self-test is used")
    try:
        print(json.dumps(audit(args.dirt_repo, args.sph_repo, args.rev), indent=2, sort_keys=True))
    except (OSError, subprocess.CalledProcessError, ValueError, urllib.error.URLError, KeyError) as error:
        print(f"AUDIT_ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
