#!/usr/bin/env python3
"""Guard a source-receipted but not-yet-reconstructible Guo replication.

The primary article is hash-bound by ``evidence_contract.py``. It states that
its periodic control-cell walls and blades are built from rigidly connected
spheres (pp. 5-6), but it omits their geometry. The assigned gravity load is
derivable from reported normal stress and periodic area, whereas wall geometry
is not. This audit therefore has no DIRT solver input. It rejects a
hypothetical plane/servo substitute and never certifies an unimplemented
sphere-built cell.
"""
import argparse
import json
import tempfile
import tomllib
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


def protocol_mismatches(config: Path) -> list[str]:
    """Return source-visible mismatches; never infer an unimplemented equivalence.

    The primary article specifies both contact geometry and normal loading.
    Replacing the former by analytic planes or the latter by a feedback servo
    changes the boundary-value problem, even if a nominal normal force happens
    to be the same.  This check deliberately identifies only disproven
    equivalence; a clean list is *not* a runnable certification.
    """
    parsed = tomllib.loads(config.read_text())
    walls = parsed.get("wall", [])
    mismatches = []
    if any(wall.get("type") == "plane" for wall in walls):
        mismatches.append(
            "Guo pp. 5-6 specifies walls/blades built from rigidly connected spheres, "
            "but the draft uses analytic plane walls"
        )
    if any("servo" in wall for wall in walls):
        mismatches.append(
            "Guo p. 5 prescribes normal stress by assigning a weight to the upper wall, "
            "but the draft uses a force-feedback servo"
        )
    return mismatches


def require_published_control_cell(config: Path) -> None:
    """Reject a known mismatch; do not certify an unimplemented source cell."""
    mismatches = protocol_mismatches(config)
    if mismatches:
        raise RuntimeError("BLOCKED: source-equivalence is false:\n- " + "\n- ".join(mismatches))
    from source_geometry_audit import require_reproducible_wall_body
    require_reproducible_wall_body()


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
    def test_rejects_plane_walls_and_servo_against_primary_protocol(self):
        with tempfile.NamedTemporaryFile("w", suffix=".toml") as config:
            config.write(
                '[[wall]]\n'
                'type = "plane"\n'
                '[wall.servo]\n'
                'target_force = 1.0\n'
            )
            config.flush()
            with self.assertRaisesRegex(
                RuntimeError, "(?s)rigidly connected spheres.*force-feedback servo"
            ):
                require_published_control_cell(Path(config.name))

    def test_audit_declaration_cannot_become_a_runnable_cell(self):
        with self.assertRaisesRegex(RuntimeError, "cannot be uniquely transcribed"):
            require_published_control_cell(HERE / "config.toml")

    def test_sphere_spelling_cannot_certify_an_unaudited_implementation(self):
        with tempfile.NamedTemporaryFile("w", suffix=".toml") as config:
            config.write('[[wall]]\ntype = "sphere"\n')
            config.flush()
            self.assertEqual(protocol_mismatches(Path(config.name)), [])
            with self.assertRaisesRegex(RuntimeError, "cannot be uniquely transcribed"):
                require_published_control_cell(Path(config.name))


if __name__ == "__main__":
    main()
