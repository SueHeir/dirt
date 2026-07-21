#!/usr/bin/env python3
"""Audit the boundary of DIRT's retired SPH glass angle-of-repose campaign.

This is deliberately an admission audit, not a numerical validation.  It checks
one repository-history fact and asks two independent live catalogues to identify
the two citations that had been near the retired campaign.  Neither citation is
admitted as a dry-glass-bead target, so this program must never print a physical
PASS result.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path


RETIRE_COMMIT = "f7fe1a44b3d744eef3b7e068c42b97c8c10ad2dc"
RETIRED_CASE = "examples/SPH_glass_sphere_calibration/03_angle_of_repose"
REQUIRED_ENTRY_POINTS = ("README.md", "main.rs", "sweep.py")

# These are intentionally rejected candidates.  The check is bibliographic:
# metadata cannot establish apparatus, bead distribution, wall friction, or
# uncertainty, so it cannot promote either work into a validation target.
REJECTED_REFERENCES = (
    {
        "doi": "10.1073/pnas.2107965118",
        "title_terms": ("expression", "angle", "repose", "cohesive", "granular"),
        "authors": {"elekes", "parteli"},
        "reason": "cohesive-material theory, not a matched dry glass-bead experiment",
    },
    {
        "doi": "10.1016/S0378-4371(99)00183-1",
        "title_terms": ("rolling", "friction", "dynamic", "simulation", "sandpile"),
        "authors": {"zhou", "wright", "yang", "xu", "yu"},
        "reason": "simulation study, not an independent primary measurement",
    },
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True, stderr=subprocess.STDOUT)


def normalize(text: str) -> str:
    return " ".join(text.lower().split())


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "dirt-retirement-audit/1"})
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)


def audit_retirement(repo: Path) -> None:
    require(git(repo, "rev-parse", "--is-inside-work-tree").strip() == "true", "--repo is not a Git worktree")
    require(
        subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", RETIRE_COMMIT, "HEAD"], check=False).returncode == 0,
        "the SPH-removal commit is not an ancestor of HEAD",
    )
    deleted = set(
        git(repo, "diff", "--name-only", "--diff-filter=D", f"{RETIRE_COMMIT}^", RETIRE_COMMIT, "--", RETIRED_CASE).splitlines()
    )
    expected = {f"{RETIRED_CASE}/{name}" for name in REQUIRED_ENTRY_POINTS}
    require(expected <= deleted, "retirement commit does not remove all historical entry points")
    present = git(repo, "ls-tree", "-r", "--name-only", "HEAD", "--", RETIRED_CASE).splitlines()
    require(not present, "retired SPH angle-of-repose case is present at HEAD")
    print(f"HISTORY CONFIRMED: {RETIRE_COMMIT[:7]} removed the SPH repose entry points; HEAD has none.")


def audit_reference(reference: dict[str, object]) -> None:
    doi = str(reference["doi"])
    encoded = urllib.parse.quote(doi, safe="")
    crossref = fetch_json(f"https://api.crossref.org/works/{encoded}")["message"]
    openalex = fetch_json(f"https://api.openalex.org/works/https://doi.org/{encoded}")
    titles = (normalize(crossref["title"][0]), normalize(openalex["title"]))
    expected_terms = tuple(reference["title_terms"])
    require(all(all(term in title for term in expected_terms) for title in titles), f"catalogue title mismatch for {doi}")
    crossref_authors = {normalize(author["family"]) for author in crossref["author"]}
    openalex_authors = {normalize(item["author"]["display_name"]).split()[-1] for item in openalex["authorships"]}
    expected_authors = set(reference["authors"])
    require(expected_authors <= crossref_authors and expected_authors <= openalex_authors, f"catalogue author mismatch for {doi}")
    print(f"CATALOGUES CONFIRM {doi}; REJECTED AS TARGET: {reference['reason']}.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--online", action="store_true", help="require Crossref and OpenAlex; fail closed otherwise")
    args = parser.parse_args()
    require(args.online, "--online is required: local citation text is not independent evidence")
    audit_retirement(args.repo.resolve())
    for reference in REJECTED_REFERENCES:
        audit_reference(reference)
    print("INCONCLUSIVE BY DESIGN: this audit validates no angle, mu_r, solver, or calibration.")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
