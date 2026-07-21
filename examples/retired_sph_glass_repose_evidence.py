#!/usr/bin/env python3
"""Audit the evidence boundary of DIRT's retired SPH repose request.

This is deliberately a forensic report, not a scientific validator. It never
scores an angle, chooses rolling friction, or returns a calibration "pass".
With ``--online`` it compares the *archived citation itself* with Crossref and
OpenAlex. Catalogue identity can expose bad provenance; it cannot establish an
experimental target or validate a solver.
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


def archived_authors(readme: str) -> list[str]:
    """Extract the first bibliography entry's surnames from the frozen README."""
    match = re.search(
        r'^1\.\s+([^\n]+?),\s+"Rolling friction in the dynamic\s+simulation of sandpile formation"',
        readme,
        flags=re.MULTILINE,
    )
    if not match:
        raise ValueError("archived README has no parseable first bibliography entry")
    names = [part.strip() for part in match.group(1).split(",")]
    surnames = [name.rsplit(" ", 1)[-1] for name in names]
    if not surnames or any(not surname for surname in surnames):
        raise ValueError("archived bibliography has no author surnames")
    return surnames


def report_local(readme: str) -> bool:
    source_revision = git("rev-parse", f"{REMOVAL_COMMIT}^").strip()
    removed = git("diff", "--name-only", f"{REMOVAL_COMMIT}^", REMOVAL_COMMIT, "--", ARCHIVED_CASE)
    live = git("ls-tree", "-r", "--name-only", "HEAD", "--", ARCHIVED_CASE)
    print("RETIREMENT EVIDENCE REPORT — NOT A VALIDATION RESULT")
    print(f"archived source: {source_revision}:{ARCHIVED_README}")
    print(f"archived source SHA-256: {hashlib.sha256(readme.encode()).hexdigest()}")
    print(f"archived numerical-claim line: {source_claim(readme)}")
    print(f"archived citation surnames: {', '.join(archived_authors(readme))}")
    print(f"files deleted by {REMOVAL_COMMIT}: {len(removed.splitlines())}")
    print(f"files now present at retired path: {len(live.splitlines())}")
    if not removed.splitlines() or live.splitlines():
        print("LOCAL SCOPE CHECK: INCONSISTENT", file=sys.stderr)
        return False
    print("LOCAL SCOPE CHECK: historical path remains absent")
    print("LIMIT: absence of this path says nothing about a renamed solver or any experiment.")
    return True


def catalogue_authors(openalex: dict) -> list[str]:
    names = []
    for authorship in openalex.get("authorships", []):
        name = authorship.get("author", {}).get("display_name", "").strip()
        if name:
            names.append(name.rsplit(" ", 1)[-1])
    return names


def report_online(readme: str) -> bool:
    encoded = urllib.parse.quote(DOI, safe="")
    crossref = fetch_json(f"https://api.crossref.org/works/{encoded}")["message"]
    openalex = fetch_json(f"https://api.openalex.org/works/https://doi.org/{encoded}")
    titles = {"Crossref": crossref.get("title", [""])[0], "OpenAlex": openalex.get("title", "")}
    archived = archived_authors(readme)
    print(f"external DOI queried: {DOI}")
    titles_match = True
    for catalogue, title in titles.items():
        matches = normalise(title) == normalise(ARCHIVED_TITLE)
        print(f"{catalogue} title: {title!r} ({'matches archived title' if matches else 'TITLE MISMATCH'})")
        titles_match &= matches
    external_authors = {
        "Crossref": [author.get("family", "?") for author in crossref.get("author", [])],
        "OpenAlex": catalogue_authors(openalex),
    }
    author_agreement = {name: authors == archived for name, authors in external_authors.items()}
    for catalogue, authors in external_authors.items():
        verdict = "matches archived citation" if author_agreement[catalogue] else "DIFFERS FROM ARCHIVED CITATION"
        print(f"{catalogue} author surnames: {', '.join(authors)} ({verdict})")
    if titles_match and not any(author_agreement.values()):
        print("PROVENANCE OUTCOME: title identity is confirmed, but both catalogues contradict the archived author list.")
    elif titles_match:
        print("PROVENANCE OUTCOME: title identity confirmed; author discrepancy not reproduced by both catalogues.")
    else:
        print("PROVENANCE OUTCOME: inconclusive/inconsistent catalogue identity.")
    print("LIMIT: this is bibliographic counterevidence only. It does not identify bead material, apparatus, contact laws, preparation, estimator, uncertainty, or an admissible angle target.")
    return titles_match


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--online", action="store_true", help="also query Crossref and OpenAlex")
    args = parser.parse_args()
    try:
        readme = git("show", f"{REMOVAL_COMMIT}^:{ARCHIVED_README}")
        local_ok = report_local(readme)
        online_ok = report_online(readme) if args.online else True
    except (OSError, ValueError, subprocess.CalledProcessError, json.JSONDecodeError) as exc:
        print(f"EVIDENCE REPORT INCONCLUSIVE: {exc}", file=sys.stderr)
        return 2
    if not (local_ok and online_ok):
        return 1
    print("No scientific or calibration conclusion follows from this report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
