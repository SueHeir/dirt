#!/usr/bin/env python3
"""Evidence eligibility for a Cundall--Strack response replication.

This module is intentionally independent of the DIRT recorder.  It decides
whether an external data set can be used to score a replication *before* a
solver response is read or a numerical tolerance is considered.
"""
from __future__ import annotations

from dataclasses import dataclass


REQUIRED_SERIES = (
    "state_registration",
    "stress_or_deviatoric_path",
    "volumetric_strain_or_dilatancy_path",
    "contact_or_fabric_evolution",
)


@dataclass(frozen=True)
class EvidenceDecision:
    candidate: str
    eligible: bool
    missing: tuple[str, ...]
    protocol_failures: tuple[str, ...]

    @property
    def reason(self) -> str:
        parts = []
        if self.missing:
            parts.append("missing external series: " + ", ".join(self.missing))
        if self.protocol_failures:
            parts.append("protocol mismatch: " + ", ".join(self.protocol_failures))
        return "; ".join(parts) if parts else "eligible"


def decide(candidate: str, support: dict[str, bool], protocol_failures=()) -> EvidenceDecision:
    """Admit only a complete, protocol-comparable external response data set."""
    unknown = set(support) - set(REQUIRED_SERIES)
    if unknown:
        raise ValueError("unknown evidence classes: " + ", ".join(sorted(unknown)))
    missing = tuple(series for series in REQUIRED_SERIES if not support.get(series, False))
    failures = tuple(protocol_failures)
    return EvidenceDecision(candidate, not missing and not failures, missing, failures)
