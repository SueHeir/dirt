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
