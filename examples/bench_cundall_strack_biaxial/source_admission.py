#!/usr/bin/env python3
"""Fail-closed source admission for the Cundall--Strack biaxial case.

The result concerns the cited primary source, not a DIRT run. A different PDF,
or a paper without each affirmative scope statement, is rejected before any
claim about replication eligibility can be made.
"""
from __future__ import annotations

import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

PRIMARY_SOURCE_SHA256 = "9fe966e2e470e66848d74a1f0c45ae157931e3b9c3566ba3a0521b3fbf0f1b84"

# Exact affirmative statements from the primary PDF. These document why the
# source is not a numerical response oracle; absence of a keyword is not used.
REQUIRED_SCOPE_STATEMENTS = {
    "qualitative verification": "the purpose of the verification is to compare force vector plots",
    "qualitative scope": "the verification is primarily qualitative",
    "digitised geometry": "centre locations and radii were obtained from fig.",
    "unknown pre-load geometry": "geometry of the assembly and the wall locations as they existed before loading are unknown",
    "control mismatch": "numerical test is strain-controlled whereas the original is stress-controlled",
}

MISSING_REPLICATION_OBSERVABLES = (
    "registered state map (disc centres/radii and beam or wall positions)",
    "boundary-control and wall-contact-law history",
    "stress/deviatoric trajectory at common states",
    "volumetric strain or dilatancy trajectory at common states",
    "contact or fabric evolution at common states",
)

# Tables 1--3 and Figs. 10--13 report these *individual wall-force snapshots*.
# They are retained solely to stop a future benchmark from silently treating a
# handful of unregistered states as a strain trajectory.  The values are not a
# curve, do not have common strain/volume coordinates, and are therefore never
# returned as an admissible validation oracle.
WALL_FORCE_SNAPSHOTS = (
    ("experiment-A", 0.39),
    ("experiment-B", 0.33),
    ("BALL-1", 0.43),
    ("BALL-2", 0.40),
    ("BALL-3", 0.46),
    ("BALL-4", 0.47),
)


class SourceAdmissionError(RuntimeError):
    """The proposed source cannot be admitted for this replication claim."""


@dataclass(frozen=True)
class AdmissionResult:
    sha256: str
    pages: int
    statement_pages: dict[str, tuple[int, ...]]
    eligible: bool
    missing_observables: tuple[str, ...]
    wall_force_snapshots: tuple[tuple[str, float], ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalise(text: str) -> str:
    return " ".join(text.lower().split())


def audit_primary_source(path: str | Path, reader_factory: Callable | None = None) -> AdmissionResult:
    """Authenticate the primary PDF and return its fail-closed admission result."""
    source = Path(path)
    actual = sha256(source)
    if actual != PRIMARY_SOURCE_SHA256:
        raise SourceAdmissionError(
            f"primary-source checksum mismatch: expected {PRIMARY_SOURCE_SHA256}, got {actual}"
        )
    if reader_factory is None:
        try:
            from pypdf import PdfReader
        except ImportError as exc:
            raise SourceAdmissionError("source audit requires pypdf; install requirements.txt") from exc
        reader_factory = PdfReader

    pages = [_normalise(page.extract_text() or "") for page in reader_factory(str(source)).pages]
    statement_pages: dict[str, tuple[int, ...]] = {}
    for label, statement in REQUIRED_SCOPE_STATEMENTS.items():
        hits = tuple(index + 1 for index, page in enumerate(pages) if statement in page)
        if not hits:
            raise SourceAdmissionError(f"primary source lacks required scope statement: {label}")
        statement_pages[label] = hits

    return AdmissionResult(
        sha256=actual,
        pages=len(pages),
        statement_pages=statement_pages,
        eligible=False,
        missing_observables=MISSING_REPLICATION_OBSERVABLES,
        wall_force_snapshots=WALL_FORCE_SNAPSHOTS,
    )


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {argv[0]} PRIMARY_SOURCE.pdf", file=sys.stderr)
        return 64
    try:
        result = audit_primary_source(argv[1])
    except (OSError, SourceAdmissionError) as exc:
        print(f"SOURCE_REJECTED: {exc}", file=sys.stderr)
        return 1
    print(f"PRIMARY_SOURCE_SHA256={result.sha256}")
    print("SCOPE_STATEMENTS=" + ", ".join(
        f"{label}@{','.join(map(str, pages))}" for label, pages in result.statement_pages.items()
    ))
    print("REPLICATION_UNAVAILABLE: cited source is qualitative and lacks registered response data")
    print("WALL_FORCE_SNAPSHOTS_ONLY=" + ", ".join(
        f"{label}:{ratio:.2f}" for label, ratio in result.wall_force_snapshots
    ))
    for missing in result.missing_observables:
        print(f"MISSING: {missing}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
