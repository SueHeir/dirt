#!/usr/bin/env python3
"""Adversarially identify the source cited for DIRT's retired SPH claim.

This is an evidence-admission audit, not a calibration test. It derives the
cited work's title from the immutable README immediately before retirement;
then Crossref discovers a DOI from that title and Crossref and OpenAlex check
the resulting record independently. No DOI, title, author list, angle, or
material property is embedded as an expected passing value in this program.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path


RETIRE_COMMIT = "f7fe1a44b3d744eef3b7e068c42b97c8c10ad2dc"
HISTORICAL_README = "examples/SPH_glass_sphere_calibration/03_angle_of_repose/README.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True, stderr=subprocess.STDOUT)


def normalize(text: str) -> str:
    return " ".join(text.lower().split())


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "dirt-retired-claim-audit/2"})
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)


def historical_readme(repo: Path) -> str:
    return git(repo, "show", f"{RETIRE_COMMIT}^:{HISTORICAL_README}")


def archived_claim_and_title(readme: str) -> tuple[str, str]:
    """Extract the retired claim and first supplied reference from source."""
    lower = normalize(readme)
    require("measured glass repose band" in lower, "historical repose claim not found")
    band_match = re.search(r"\[\s*(\d+(?:\.\d+)?)\s*,\s*(\d+(?:\.\d+)?)\s*\]", readme)
    require(band_match is not None, "historical repose claim has no numeric band")
    band = f"[{band_match.group(1)}, {band_match.group(2)}]"
    references = readme.split("## References", 1)
    require(len(references) == 2, "historical README has no reference section")
    match = re.search(r'^1\.\s+[\s\S]*?"([^"]+)"', references[1], flags=re.MULTILINE)
    require(match is not None, "historical first reference has no quoted title")
    title = normalize(match.group(1))
    require(title, "historical first reference has an empty title")
    return band, title


def crossref_discover(title: str) -> dict:
    """Discover exactly one Crossref work by the source-derived exact title."""
    query = urllib.parse.urlencode({"query.bibliographic": title, "rows": 20})
    candidates = fetch_json(f"https://api.crossref.org/works?{query}")["message"]["items"]
    matches = [item for item in candidates if normalize(" ".join(item.get("title", []))) == title]
    require(len(matches) == 1, f"Crossref did not uniquely resolve archived title ({len(matches)} exact matches)")
    doi = matches[0].get("DOI")
    require(isinstance(doi, str) and doi, "Crossref match has no DOI")
    return matches[0]


def audit_catalogues(title: str) -> None:
    """Use Crossref discovery plus two direct catalogue records as witnesses."""
    discovered = crossref_discover(title)
    doi = discovered["DOI"]
    encoded_doi = urllib.parse.quote(doi, safe="")
    crossref_work = fetch_json(f"https://api.crossref.org/works/{encoded_doi}")["message"]
    openalex_work = fetch_json(f"https://api.openalex.org/works/https://doi.org/{encoded_doi}")
    crossref_title = normalize(" ".join(crossref_work["title"]))
    openalex_title = normalize(openalex_work["title"])
    require(crossref_title == openalex_title, "catalogues disagree on cited-paper title")
    require(crossref_title == title, "catalogue title does not match archived citation")
    require("simulation" in crossref_title, "archived citation is not self-described as a simulation")
    print(f"SOURCE-DERIVED CATALOGUE CONSENSUS: {crossref_title} ({doi}).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    repo = args.repo.resolve()
    require(git(repo, "rev-parse", "--is-inside-work-tree").strip() == "true", "--repo is not a Git worktree")
    readme = historical_readme(repo)
    band, title = archived_claim_and_title(readme)
    audit_catalogues(title)
    print(f"INADMISSIBLE BY DESIGN: archived band {band} is paired with a self-described simulation, not a primary dry-glass-bead measurement; no angle, mu_r, solver, or calibration is validated.")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
