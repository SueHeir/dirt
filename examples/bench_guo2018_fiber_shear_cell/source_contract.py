#!/usr/bin/env python3
"""Guard a source-receipted but geometrically non-equivalent Guo draft.

The primary article is now hash-bound by ``evidence_contract.py``.  It states
that its periodic control-cell walls and blades are built from rigidly connected
spheres (pp. 5-6).  The retained DIRT input uses finite plane walls, so it is
not source-equivalent and must not run as a Guo replication.  This is neither
a physics validation nor a substitute for the two published-observable gates.
"""
import argparse
import json
import unittest
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOI = "10.1002/aic.16397"
CROSSREF = f"https://api.crossref.org/works/{DOI}"


def published_record() -> dict:
    """Fetch bibliographic identity only; geometry comes from the primary paper."""
    with urllib.request.urlopen(CROSSREF, timeout=30) as response:
        record = json.load(response)["message"]
    title = " ".join(record.get("title", []))
    if record.get("DOI", "").lower() != DOI or "Flexible Fiber Flows" not in title:
        raise ValueError("Crossref record does not identify Guo et al. doi:10.1002/aic.16397")
    return {"title": title, "doi": record["DOI"]}


def require_published_control_cell(config: Path) -> None:
    """Reject the retained plane-wall draft against the paper's wall topology."""
    text = config.read_text()
    if 'type = "plane"' in text:
        raise RuntimeError(
            "BLOCKED: Guo pp. 5-6 specifies walls/blades built from rigidly connected spheres, "
            "but config.toml uses plane walls; this is not a source-equivalent Guo control cell"
        )
    raise RuntimeError(
        "BLOCKED: source-equivalence is not established for this control-cell configuration; "
        "do not run it as a Guo replication"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.toml")
    parser.add_argument("--verify-doi", action="store_true", help="query Crossref for bibliographic identity")
    parser.add_argument("--require-runnable", action="store_true", help="raise unless the configured numerical cell matches the source protocol")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(SourceContractTests)
        raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    if not args.verify_doi and not args.require_runnable:
        parser.error("choose --verify-doi and/or --require-runnable")
    if args.verify_doi:
        record = published_record()
        print(f"EXTERNAL BIBLIOGRAPHY: {record['doi']} — {record['title']}")
        print("Bibliography is independently checked; method evidence is hash-bound in reference_provenance.json.")
    if args.require_runnable:
        require_published_control_cell(args.config)


class SourceContractTests(unittest.TestCase):
    def test_rejects_plane_walls_against_primary_wall_topology(self):
        with self.assertRaisesRegex(RuntimeError, "rigidly connected spheres"):
            require_published_control_cell(HERE / "config.toml")


if __name__ == "__main__":
    main()
