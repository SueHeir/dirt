#!/usr/bin/env python3
"""Admission ledger for an independent biaxial response candidate.

An external solver output is not automatically a reference.  This small,
data-driven boundary records the apparatus features that must be common before
the response-series contract is even considered.  In particular it prevents a
periodic/virial trace from being relabelled as a finite moving-wall resultant
comparison merely because it has the same particle count or friction.
"""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path


REQUIRED_PROTOCOL_FIELDS = (
    "dimension_and_out_of_plane_confinement",
    "lateral_boundary_and_resultant",
    "axial_boundary_and_resultant",
    "particle_and_wall_contact_law",
    "preparation",
    "loading_path",
    "state_coordinate",
    "stress_or_deviatoric_path",
    "volumetric_strain_or_dilatancy_path",
    "contact_or_fabric_evolution",
)

# The archived deck is a negative control, not an admissible comparison.  Bind
# the ledger to its exact bytes so the protocol findings can be independently
# reproduced and cannot silently drift with an untracked local input.
ARCHIVED_LAMMPS_SHA256 = "e69d0b4e99b102b4e949d4172d9862606277a44d5f7b1e2452fcfb2f24928777"
# The ledger is evidence too: without a content binding, a later edit could
# silently sever its claimed protocol facts from the archived deck.  This is a
# negative-control receipt, not a checksum for a positive reference.
ARCHIVED_LAMMPS_LEDGER_SHA256 = "5a6cb3b6b5b5e08b15e20a774aa09db6b0ef381f45779936974148fad4f5cf94"


def _sha256(path: str | Path) -> str:
    """Return a content address without interpreting an evidence artifact."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit_archived_lammps_deck(path: str | Path) -> None:
    """Authenticate and structurally verify the rejected periodic input.

    This audit intentionally verifies only the facts used to reject the
    candidate.  It cannot make the candidate eligible: the response-series
    fields remain absent from the ledger and are separately required.
    """
    raw = Path(path).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != ARCHIVED_LAMMPS_SHA256:
        raise RuntimeError(
            "archived LAMMPS negative-control checksum mismatch: "
            f"expected {ARCHIVED_LAMMPS_SHA256}, got {actual}"
        )
    text = " ".join(raw.decode("utf-8").lower().split())
    required_fragments = (
        "dimension 2",
        "boundary p p p",
        "fix planar all enforce2d",
        "fix seat all press/berendsen",
        "fix drive all deform 1 y erate -25.0",
        "thermo_style custom step temp press v_e v_sxx v_syy lx ly atoms",
    )
    absent = [fragment for fragment in required_fragments if fragment not in text]
    if absent:
        raise RuntimeError("archived LAMMPS deck lacks expected protocol facts: " + ", ".join(absent))


def audit_archived_lammps_ledger(path: str | Path) -> None:
    """Authenticate the complete rejection ledger before interpreting it.

    Authentication is intentionally performed before CSV parsing.  The ledger
    cannot be selectively amended to make a periodic candidate appear more
    comparable while retaining the same archived-deck receipt.
    """
    actual = _sha256(path)
    if actual != ARCHIVED_LAMMPS_LEDGER_SHA256:
        raise RuntimeError(
            "archived LAMMPS negative-control ledger checksum mismatch: "
            f"expected {ARCHIVED_LAMMPS_LEDGER_SHA256}, got {actual}"
        )


@dataclass(frozen=True)
class CandidateDecision:
    candidate: str
    eligible: bool
    failures: tuple[str, ...]

    @property
    def reason(self) -> str:
        return "eligible" if not self.failures else "protocol/evidence mismatch: " + ", ".join(self.failures)


def decide(path: str | Path) -> CandidateDecision:
    """Admit only a complete, apparatus-equivalent independently produced trace."""
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"candidate", "field", "common", "evidence", "limitation"}
    if not rows or any(set(row) != required for row in rows):
        raise RuntimeError("invalid external-candidate admission ledger schema")
    candidates = {row["candidate"] for row in rows}
    if len(candidates) != 1:
        raise RuntimeError("each admission ledger must describe exactly one candidate")
    fields = [row["field"] for row in rows]
    if set(fields) != set(REQUIRED_PROTOCOL_FIELDS) or len(fields) != len(set(fields)):
        raise RuntimeError("admission ledger must cover each required protocol/evidence field once")
    failures = tuple(row["field"] for row in rows if row["common"] != "yes")
    if any(not row["evidence"].strip() or not row["limitation"].strip() for row in rows):
        raise RuntimeError("admission ledger requires evidence and limitation for every field")
    return CandidateDecision(next(iter(candidates)), not failures, failures)
