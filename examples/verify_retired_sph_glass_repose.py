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

# A DIRT retirement audit must not silently turn a similarly named file in a
# neighbouring solver into a "replacement".  These tokens intentionally look
# only for a runnable repose/calibration surface, not for generic friction code:
# granular SPH may quite properly have frictional constitutive parameters
# without implementing this experiment.
REPLACEMENT_TERMS = ("angle_of_repose", "angle-of-repose", "angle of repose", "repose_calibration")

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


def retired_tree_paths(repo: Path) -> set[str]:
    """Return the complete historical campaign manifest, not just its runner.

    Checking only a README, binary, and sweep leaves a misleadingly easy path
    for a partial restoration: an old config, result table, or helper could be
    brought back under the retired directory and later wired into a new claim.
    The removal commit itself is the independent witness of what the campaign
    comprised, so use that complete tree as the negative-scope manifest.
    """
    return set(
        git(repo, "ls-tree", "-r", "--name-only", f"{RETIRE_COMMIT}^", "--", RETIRED_CASE).splitlines()
    )


def tree_blobs(repo: Path, revision: str, pathspec: str | None = None) -> dict[str, str]:
    """Return regular-file paths and Git blob IDs from an immutable tree."""
    args = ["ls-tree", "-r", revision]
    if pathspec is not None:
        args.extend(("--", pathspec))
    blobs: dict[str, str] = {}
    for line in git(repo, *args).splitlines():
        metadata, path = line.split("\t", 1)
        mode, kind, blob = metadata.split()
        if mode in {"100644", "100755"} and kind == "blob":
            blobs[path] = blob
    return blobs


def require_retired_surface_absent(
    historical: set[str], deleted_by_retirement: set[str], present_at_head: set[str]
) -> None:
    """Fail if removal was partial or if any historical campaign file returns."""
    require(historical, "retirement predecessor contains no SPH repose campaign files")
    missing_from_removal = sorted(historical - deleted_by_retirement)
    require(
        not missing_from_removal,
        "retirement commit leaves historical SPH repose files behind: " + ", ".join(missing_from_removal),
    )
    restored = sorted(historical & present_at_head)
    require(
        not restored,
        "retired SPH repose files have been restored at HEAD: " + ", ".join(restored),
    )


def require_no_relocated_historical_blobs(historical: dict[str, str], present_at_head: dict[str, str]) -> None:
    """Reject a byte-for-byte relocation of the retired campaign.

    A path-only check cannot distinguish deletion from moving an old runner or
    generated input to a new directory.  This check is deliberately narrower
    than a semantic source-code search: it rejects only a retained historical
    Git blob, and does not assert that no independently implemented experiment
    exists elsewhere.
    """
    historical_blob_ids = set(historical.values())
    relocated = sorted(path for path, blob in present_at_head.items() if blob in historical_blob_ids)
    require(
        not relocated,
        "historical SPH repose campaign blobs were relocated at HEAD: " + ", ".join(relocated),
    )


def audit_retirement(repo: Path) -> None:
    require(git(repo, "rev-parse", "--is-inside-work-tree").strip() == "true", "--repo is not a Git worktree")
    require(
        subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", RETIRE_COMMIT, "HEAD"], check=False).returncode == 0,
        "the SPH-removal commit is not an ancestor of HEAD",
    )
    deleted = set(
        git(repo, "diff", "--name-only", "--diff-filter=D", f"{RETIRE_COMMIT}^", RETIRE_COMMIT, "--", RETIRED_CASE).splitlines()
    )
    present = set(git(repo, "ls-tree", "-r", "--name-only", "HEAD", "--", RETIRED_CASE).splitlines())
    historical = retired_tree_paths(repo)
    require_retired_surface_absent(historical, deleted, present)
    historical_blobs = tree_blobs(repo, f"{RETIRE_COMMIT}^", RETIRED_CASE)
    require(set(historical_blobs) == historical, "historical SPH repose manifest/blob listing disagrees")
    require_no_relocated_historical_blobs(historical_blobs, tree_blobs(repo, "HEAD"))
    print(
        f"HISTORY CONFIRMED: {RETIRE_COMMIT[:7]} removed all {len(historical)} SPH repose campaign files; "
        "HEAD has neither those paths nor byte-identical relocated blobs."
    )


def replacement_candidates(paths: set[str], texts: dict[str, str]) -> list[str]:
    """Return candidate runnable surfaces, conservatively, without deciding physics.

    This is a negative-scope guard, not a search for a physical reference.  A
    hit requires human scientific review; no hit says only that the checked
    tree has no obvious replacement under the retired experiment's vocabulary.
    """
    runnable_suffixes = (".rs", ".py", ".toml")
    candidates = {
        path
        for path in paths
        if path.endswith(runnable_suffixes) and any(term in path.lower() for term in REPLACEMENT_TERMS)
    }
    for path, text in texts.items():
        lowered = text.lower()
        if path.endswith(runnable_suffixes) and any(term in lowered for term in REPLACEMENT_TERMS):
            candidates.add(path)
    return sorted(candidates)


def audit_maintained_sph_boundary(soil_sph: Path) -> None:
    """Check the named maintained SPH checkout for an obvious runnable replacement.

    We deliberately inspect Git's committed HEAD, not a caller's working tree,
    so an uncommitted local experiment cannot be mistaken for maintained work.
    """
    require(
        git(soil_sph, "rev-parse", "--is-inside-work-tree").strip() == "true",
        "--soil-sph is not a Git worktree",
    )
    revision = git(soil_sph, "rev-parse", "HEAD").strip()
    paths = set(git(soil_sph, "ls-tree", "-r", "--name-only", revision).splitlines())
    tracked_text = {
        path: git(soil_sph, "show", f"{revision}:{path}")
        for path in paths
        if path.endswith((".toml", ".rs", ".py"))
    }
    candidates = replacement_candidates(paths, tracked_text)
    require(
        not candidates,
        "maintained SPH checkout has a candidate repose/calibration surface; "
        "a human must assess it before claiming retirement: " + ", ".join(candidates),
    )
    print(
        f"MAINTAINED-SPH BOUNDARY: {soil_sph} at {revision[:12]} has no obvious committed "
        "repose/calibration replacement; this is not a solver validation."
    )


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
    parser.add_argument(
        "--soil-sph",
        type=Path,
        help="optional maintained dev_soil_sph checkout to audit at its committed HEAD",
    )
    args = parser.parse_args()
    require(args.online, "--online is required: local citation text is not independent evidence")
    audit_retirement(args.repo.resolve())
    for reference in REJECTED_REFERENCES:
        audit_reference(reference)
    if args.soil_sph is not None:
        audit_maintained_sph_boundary(args.soil_sph.resolve())
    print("INCONCLUSIVE BY DESIGN: this audit validates no angle, mu_r, solver, or calibration.")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
