#!/usr/bin/env python3
"""Print reproducible provenance facts for DIRT's retired SPH repose request.

This is deliberately a forensic report, not a scientific validator.  It never
scores an angle, chooses a rolling-friction value, or returns a calibration
"pass".  With ``--online`` it compares bibliographic identity from two external
catalogues to the title stored in the pre-removal Git object.  That comparison
does not establish experimental suitability.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
REMOVAL_COMMIT = "f7fe1a4"
ARCHIVED_README = "examples/SPH_glass_sphere_calibration/03_angle_of_repose/README.md"
ARCHIVED_CASE = "examples/SPH_glass_sphere_calibration/03_angle_of_repose"
DOI = "10.1016/S0378-4371(99)00183-1"
ARCHIVED_TITLE = "Rolling friction in the dynamic simulation of sandpile formation"


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO_ROOT, text=True, stderr=subprocess.STDOUT
    )


def normalise(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "DIRT-retirement-audit/1"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def source_claim(readme: str) -> str:
    for line in readme.splitlines():
        if "[22, 26]" in line:
            return line.strip()
    raise ValueError("archived README has no [22, 26] claim")


def report_local(readme: str) -> bool:
    source_revision = git("rev-parse", f"{REMOVAL_COMMIT}^").strip()
    removed = git("diff", "--name-only", f"{REMOVAL_COMMIT}^", REMOVAL_COMMIT, "--", ARCHIVED_CASE)
    live = git("ls-tree", "-r", "--name-only", "HEAD", "--", ARCHIVED_CASE)
    print("RETIREMENT EVIDENCE REPORT — NOT A VALIDATION RESULT")
    print(f"archived source: {source_revision}:{ARCHIVED_README}")
    print(f"archived source SHA-256: {hashlib.sha256(readme.encode()).hexdigest()}")
    print(f"archived numerical-claim line: {source_claim(readme)}")
    print(f"files deleted by {REMOVAL_COMMIT}: {len(removed.splitlines())}")
    print(f"files now present at retired path: {len(live.splitlines())}")
    if not removed.splitlines() or live.splitlines():
        print("LOCAL SCOPE CHECK: INCONSISTENT", file=sys.stderr)
        return False
    print("LOCAL SCOPE CHECK: historical path remains absent")
    print("LIMIT: absence of this path says nothing about a renamed solver or any experiment.")
    return True


def report_online() -> bool:
    encoded = urllib.parse.quote(DOI, safe="")
    crossref = fetch_json(f"https://api.crossref.org/works/{encoded}")["message"]
    openalex = fetch_json(f"https://api.openalex.org/works/https://doi.org/{encoded}")
    titles = {"Crossref": crossref.get("title", [""])[0], "OpenAlex": openalex.get("title", "")}
    print(f"external DOI queried: {DOI}")
    ok = True
    for catalogue, title in titles.items():
        matches = normalise(title) == normalise(ARCHIVED_TITLE)
        print(f"{catalogue} title: {title!r} ({'matches archived title' if matches else 'TITLE MISMATCH'})")
        ok &= matches
    authors = ", ".join(a.get("family", "?") for a in crossref.get("author", []))
    print(f"Crossref author surnames: {authors}")
    print("CATALOGUE IDENTITY CHECK: " + ("consistent" if ok else "inconclusive/inconsistent"))
    print("LIMIT: catalogue metadata is not evidence of bead material, apparatus, contact laws, preparation, estimator, uncertainty, or an admissible target.")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--online", action="store_true", help="also query Crossref and OpenAlex")
    args = parser.parse_args()
    try:
        readme = git("show", f"{REMOVAL_COMMIT}^:{ARCHIVED_README}")
        local_ok = report_local(readme)
        online_ok = report_online() if args.online else True
    except (OSError, ValueError, subprocess.CalledProcessError, json.JSONDecodeError) as exc:
        print(f"EVIDENCE REPORT INCONCLUSIVE: {exc}", file=sys.stderr)
        return 2
    if not (local_ok and online_ok):
        return 1
    print("No scientific or calibration conclusion follows from this report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
