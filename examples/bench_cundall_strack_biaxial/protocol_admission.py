#!/usr/bin/env python3
"""Fail closed before a cross-code trace is presented as a replication score."""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "reference", "lammps", "protocol_admission.csv")


def admission_failures(path=MANIFEST):
    """Return required protocol fields that are not common between solvers."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"field", "dirt", "lammps", "required_for_scored_replication"}
    if not rows or any(set(row) != required for row in rows):
        raise RuntimeError("invalid cross-code protocol admission manifest")
    if len({row["field"] for row in rows}) != len(rows):
        raise RuntimeError("protocol admission fields must be unique")
    return [row["field"] for row in rows
            if row["required_for_scored_replication"] == "yes"
            and row["dirt"].strip() != row["lammps"].strip()]


if __name__ == "__main__":
    failed = admission_failures()
    if failed:
        raise SystemExit("INELIGIBLE cross-code replication protocol: " + ", ".join(failed))
