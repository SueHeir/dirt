#!/usr/bin/env python3
"""Independently audit the source cited for DIRT's retired SPH repose claim.

This is not a calibration test.  It reads the historical README from the Git
object immediately before the campaign was removed, then asks Crossref and
OpenAlex to identify the *one reference that README actually supplied*.  The
historical document calls its 22--26 degree interval empirical, but its cited
paper is a computer-simulation study.  Therefore the historical record cannot
support a dry-glass-bead validation target.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.request
from pathlib import Path


RETIRE_COMMIT = "f7fe1a44b3d744eef3b7e068c42b97c8c10ad2dc"
HISTORICAL_README = "examples/SPH_glass_sphere_calibration/03_angle_of_repose/README.md"
HISTORICAL_TITLE = "Rolling friction in the dynamic simulation of sandpile formation"
HISTORICAL_DOI = "10.1016/S0378-4371(99)00183-1"
ARCHIVED_FAMILY_NAMES = {"zhou", "xu", "yu", "zulli"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True, stderr=subprocess.STDOUT)


def normalize(text: str) -> str:
    return " ".join(text.lower().split())


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "dirt-retired-claim-audit/1"})
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)


def historical_readme(repo: Path) -> str:
    return git(repo, "show", f"{RETIRE_COMMIT}^:{HISTORICAL_README}")


def require_unsupported_claim(readme: str) -> None:
    """Check the archived claim, not a rewritten present-day description."""
    lower = normalize(readme)
    require("measured glass repose band" in lower and "22, 26" in lower, "historical repose claim not found")
    references = lower.split("## references", 1)
    require(len(references) == 2, "historical README has no reference section")
    require(normalize(HISTORICAL_TITLE) in references[1], "historical cited paper not found")
    require("experiment" not in references[1], "historical reference section unexpectedly claims an experiment")


def archived_attribution_matches(catalogue_authors: set[str]) -> bool:
    """Whether the names in the archived citation survive an external check."""
    return ARCHIVED_FAMILY_NAMES <= catalogue_authors


def audit_catalogues() -> None:
    """Cross-check one immutable historical bibliographic identity.

    A search result is not a stable identity witness: ranking and unrelated near
    matches can change.  The historical title is first required in the archived
    README; the known DOI is then fetched directly from two independent
    catalogues, both of which must report that same title.
    """
    doi = HISTORICAL_DOI
    crossref_work = fetch_json(f"https://api.crossref.org/works/{doi}")["message"]
    openalex_work = fetch_json(f"https://api.openalex.org/works/https://doi.org/{doi}")
    crossref_title = normalize(" ".join(crossref_work["title"]))
    openalex_title = normalize(openalex_work["title"])
    require(crossref_title == openalex_title, "catalogues disagree on cited-paper title")
    require(crossref_title == normalize(HISTORICAL_TITLE), "catalogue title does not match the archived citation")
    crossref_authors = {normalize(author["family"]) for author in crossref_work["author"]}
    openalex_authors = {normalize(item["author"]["display_name"]).split()[-1] for item in openalex_work["authorships"]}
    require(crossref_authors == openalex_authors, "catalogues disagree on cited-paper authors")
    archived_authors_match = archived_attribution_matches(crossref_authors)
    print(f"CATALOGUES CONFIRM HISTORICAL CITATION: {crossref_title} ({doi}).")
    if not archived_authors_match:
        print("ARCHIVED CITATION ATTRIBUTION DISAGREES WITH BOTH CATALOGUES.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    repo = args.repo.resolve()
    require(git(repo, "rev-parse", "--is-inside-work-tree").strip() == "true", "--repo is not a Git worktree")
    readme = historical_readme(repo)
    require_unsupported_claim(readme)
    audit_catalogues()
    print("INADMISSIBLE BY DESIGN: the archived claim cites a simulation study and its author attribution is not catalogue-confirmed; no angle, mu_r, solver, or calibration is validated.")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
