#!/usr/bin/env python3
"""Run and independently check DIRT's Cundall--Strack-inspired wall cell.

The paper gives two *configurations* (Fig. 10), but does not report a loading
history or a strain/state marker that identifies either configuration.  A
DIRT time/strain window therefore cannot be called source stage A or B.  This
driver checks the live, directly measured wall reactions and deliberately has
no external PASS claim.
"""
import csv
import math
import os
import subprocess
import sys
from replication_contract import REQUIRED_SERIES, decide

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(HERE, "data", "biaxial_results.csv")
REFERENCE = os.path.join(HERE, "data", "cundall_strack_stages.csv")
REGISTRATION = os.path.join(HERE, "data", "source_state_registration.csv")
PROTOCOL = os.path.join(HERE, "data", "cundall_strack_protocol.csv")
EVIDENCE_INVENTORY = os.path.join(HERE, "data", "external_evidence_inventory.csv")
PLOTS = os.path.join(HERE, "plots")
REQUIRED = ("axial_strain", "f_h_mean", "f_v_mean", "wall_force_ratio", "contacts")
# These are executable specimen-integrity floors, not external response
# tolerances.  A loose random-insertion transient is not a dense assembly.
MIN_DENSE_PHI = 0.50
MIN_DENSE_COORDINATION = 4.0
# Every recorder value is serialized with ``{:.8e}``, so this is an output
# round-trip bound (not a physical or source-comparison tolerance).
CSV_RELATIVE_PRECISION = 5.0e-8


def build_run():
    subprocess.run(["cargo", "build", "--release", "--no-default-features",
                    "--features", "precision-double", "--example",
                    "bench_cundall_strack_biaxial"], cwd=ROOT, check=True)
    if os.path.exists(RESULTS):
        os.remove(RESULTS)
    subprocess.run([os.path.join(ROOT, "target", "release", "examples",
                                 "bench_cundall_strack_biaxial"),
                    os.path.join(HERE, "config.toml")], cwd=ROOT, check=True)


def read(path):
    with open(path, newline="") as f:
        return [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(f)]



def read_reference(path=REFERENCE):
    """Read the primary-source transcription without converting it into an oracle.

    Fig. 10 supplies two labelled configurations, not a trajectory coordinate.
    Keeping this parser separate from ``evaluate`` makes that boundary testable:
    a future registered protocol may consume the rows, but this runner cannot
    silently turn its own time windows into Fig. 10 stages.
    """
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    selected = [r for r in rows if r["figure"].startswith("Fig. 10")]
    values = {r["stage"]: float(r["fh_over_fv"]) for r in selected}
    if values != {"A": 0.39, "B": 0.33}:
        raise RuntimeError("source transcription must contain Fig. 10 A=0.39 and B=0.33")
    if any("photoelastic" not in r["method"].lower() for r in selected):
        raise RuntimeError("Fig. 10 rows must identify the photoelastic experiment")
    return values


def read_protocol(path=PROTOCOL):
    """Audit the source facts that are sufficient to define the A-to-B load event.

    This deliberately does *not* promote the event into a DIRT state map: the
    paper says that it digitized the Fig. 10(a) disc locations, but does not
    publish those coordinates or the four initial wall positions.
    """
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    by_stage = {r["stage"]: r for r in rows}
    if set(by_stage) != {"A", "B"}:
        raise RuntimeError("source protocol must contain A and B")
    a, b = by_stage["A"], by_stage["B"]
    if (int(a["disc_count"]), int(b["disc_count"])) != (197, 197):
        raise RuntimeError("source protocol must retain the 197-disc census")
    if a["radius_census"] != "40:6;35:7;30:16;25:33;20:33;18:36;15:28;10:38":
        raise RuntimeError("unexpected source radius census")
    if abs(float(b["vertical_load_factor"]) - 1.04348) > 1e-8 or \
       abs(float(b["horizontal_load_factor"]) - 0.95652) > 1e-8:
        raise RuntimeError("source A-to-B load factors must be transcribed exactly")
    return by_stage


def audit_external_evidence(path=EVIDENCE_INVENTORY):
    """Audit whether the cited source can support this goal's required observables.

    This is deliberately an *evidence* audit, not a DIRT response criterion.
    A source comparison is only eligible when its state map, stress path,
    dilatancy path, and fabric/contact evolution are independently available.
    The primary paper's Fig. 10 has none of those trajectory artifacts, so the
    negative result here is expected and prevents a snapshot label from being
    promoted into a synthetic stress--strain validation.
    """
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    required = {
        "state_registration",
        "stress_or_deviatoric_path",
        "volumetric_strain_or_dilatancy_path",
        "contact_or_fabric_evolution",
    }
    if {row.get("observable") for row in rows} != required:
        raise RuntimeError("external evidence inventory must cover every replication observable")
    if any(row.get("required_for_replication") != "yes" for row in rows):
        raise RuntimeError("external evidence inventory must not downgrade a goal observable")
    if any(not row.get("source_location", "").strip() or not row.get("limitation", "").strip()
           for row in rows):
        raise RuntimeError("external evidence inventory requires source locations and limitations")
    missing = [row["observable"] for row in rows if row.get("source_support") == "no"]
    unsupported = [row["observable"] for row in rows if row.get("source_support") not in {"yes", "no"}]
    if unsupported:
        raise RuntimeError("external evidence support must be yes or no")
    return missing


def replication_evidence_decision():
    """Classify the cited source before any response score is attempted.

    The former bundled LAMMPS curve has been deliberately removed: its
    periodic, 2-D virial measurement was not the finite-wall, 3-D resultant
    protocol measured here.  Keeping it in a response plot invited a visual
    comparison that the protocol contract already rejected.
    """
    primary_missing = set(audit_external_evidence())
    primary = decide("Cundall--Strack 1979 primary source", {
        series: series not in primary_missing for series in REQUIRED_SERIES
    })
    return primary


def source_registration(path=REGISTRATION):
    """Load an independently justified state map or fail closed.

    Figure 10 supplies neither a loading history nor a coordinate that can
    identify A/B in a DIRT trace.  In particular, force ratio is forbidden as
    a selector: selecting the closest ratio would make the comparison circular.
    """
    if not os.path.exists(path):
        raise RuntimeError(
            "external validation unavailable: no traceable Fig. 10 state "
            "registration; refusing to select DIRT states from force ratio"
        )
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    required = {"stage", "source_state_coordinate", "dirt_state_coordinate", "basis"}
    if not rows or any(not required.issubset(row) for row in rows):
        raise RuntimeError("invalid source-state registration schema")
    if {row["stage"] for row in rows} != {"A", "B"}:
        raise RuntimeError("source-state registration must map both Fig. 10 stages")
    if any(not row["basis"].strip() for row in rows):
        raise RuntimeError("source-state registration requires a non-empty basis")
    return rows


def evaluate(rows):
    """Validate measurement integrity, not agreement with source values."""
    if len(rows) < 2:
        raise RuntimeError("expected at least two recorded states")
    if any(any(key not in row or not math.isfinite(row[key]) for key in REQUIRED)
           for row in rows):
        raise RuntimeError("non-finite or incomplete recorder output")

    loaded = [row for row in rows if row["f_v_mean"] > 0.0]
    if not loaded:
        raise RuntimeError("the loading platen recorded no positive reaction")

    # This recomputation intentionally uses the separately emitted mean
    # reactions, not the recorder's precomputed ratio.  CSV formatting bounds
    # the comparison rather than furnishing a material-response tolerance.
    residuals = [abs(row["wall_force_ratio"] - row["f_h_mean"] / row["f_v_mean"])
                 for row in loaded]
    max_ratio_residual = max(residuals)
    ratio_recomputed = all(
        residual <= CSV_RELATIVE_PRECISION * max(1.0, abs(row["wall_force_ratio"]))
        for residual, row in zip(residuals, loaded)
    )
    checks = {
        "rows": len(rows),
        "forward_compression": rows[-1]["axial_strain"] > rows[0]["axial_strain"],
        "lateral_compression": rows[-1]["lateral_strain"] < rows[0]["lateral_strain"],
        "positive_contacts": max(row["contacts"] for row in rows) > 0.0,
        "dense_solid_fraction": min(row["phi"] for row in rows) >= MIN_DENSE_PHI,
        "dense_coordination": min(row["coordination"] for row in rows) >= MIN_DENSE_COORDINATION,
        "positive_platen_reaction": True,
        "ratio_recomputed": ratio_recomputed,
        "max_ratio_residual": max_ratio_residual,
        "ratio_min": min(row["wall_force_ratio"] for row in loaded),
        "ratio_max": max(row["wall_force_ratio"] for row in loaded),
    }
    return all(checks[key] for key in ("forward_compression", "lateral_compression", "positive_contacts",
                                        "dense_solid_fraction", "dense_coordination",
                                        "positive_platen_reaction", "ratio_recomputed")), checks


def plot(rows, checks, passed):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f"(matplotlib unavailable; plot not regenerated: {exc})")
        return

    os.makedirs(PLOTS, exist_ok=True)
    x = [row["axial_strain"] for row in rows]
    fig, (state, ratio) = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
    state.plot(x, [row["volumetric_strain"] for row in rows],
               label="DIRT volumetric strain")
    state.plot(x, [row["fabric_anisotropy"] for row in rows],
               label="DIRT contact-fabric anisotropy")
    state.plot(x, [row["coordination"] for row in rows],
               label="DIRT coordination")
    state.set_ylabel("DIRT state observables")
    state.legend(fontsize=7, loc="best")
    state.grid(alpha=.25)

    ratio.plot(x, [row["wall_force_ratio"] for row in rows], "o-", ms=2.5,
               label="DIRT direct x-wall/platen resultant")
    ratio.set_ylabel(r"$F_H/F_V$")
    ratio.set_title("Wall measurement integrity: " + ("PASS" if passed else "FAIL") +
                    "; recorder integrity only")
    ratio.legend(fontsize=7)
    ratio.grid(alpha=.25)
    ratio.set_xlabel("axial platen strain from first recorder state")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "stress_volume_response.png"), dpi=180)
    plt.close(fig)


def main():
    command = sys.argv[1] if len(sys.argv) > 1 else "all"
    if command == "audit":
        source = read_reference()
        protocol = read_protocol()
        missing = audit_external_evidence()
        decision = replication_evidence_decision()
        print("SOURCE-TRAJECTORY REPLICATION INELIGIBLE: " + ", ".join(missing))
        print("audited source facts only: Fig. 10 A={A:.2f}, B={B:.2f}; "
              "A-to-B V={v}, H={h}".format(
                  A=source["A"], B=source["B"],
                  v=protocol["B"]["vertical_load_factor"],
                  h=protocol["B"]["horizontal_load_factor"]))
        print(f"evidence decision: {decision.candidate}: INELIGIBLE; {decision.reason}")
        return
    if command == "external":
        decision = replication_evidence_decision()
        print(f"{decision.candidate}: INELIGIBLE; {decision.reason}")
        raise SystemExit("external replication unavailable: no traceable complete trajectory")
    if command in ("all", "run"):
        build_run()
    if command in ("all", "graph"):
        rows = read(RESULTS)
        source = read_reference()
        protocol = read_protocol()
        missing = audit_external_evidence()
        decision = replication_evidence_decision()
        passed, checks = evaluate(rows)
        plot(rows, checks, passed)
        print("wall-reaction measurement: " + "; ".join(
            f"{key}={value}" for key, value in checks.items()))
        if not passed:
            raise SystemExit("WALL-REACTION MEASUREMENT FAILED")
        print("MEASUREMENT CHECK PASSED; source A-to-B protocol audited "
              f"(V={protocol['B']['vertical_load_factor']}, H={protocol['B']['horizontal_load_factor']}); "
              "Fig. 10 values audited but not used as targets "
              "because the paper provides no state-registration rule")
        print("external replication unavailable; missing source evidence: "
              + ", ".join(missing))
        print(f"evidence decision: {decision.candidate}: INELIGIBLE; {decision.reason}")
    if command not in ("all", "run", "graph", "audit", "external"):
        raise SystemExit(f"unknown command {command!r}; use all, run, graph, audit, or external")


if __name__ == "__main__":
    main()
