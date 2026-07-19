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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(HERE, "data", "biaxial_results.csv")
REFERENCE = os.path.join(HERE, "data", "cundall_strack_stages.csv")
REGISTRATION = os.path.join(HERE, "data", "source_state_registration.csv")
PROTOCOL = os.path.join(HERE, "data", "cundall_strack_protocol.csv")
EVIDENCE_INVENTORY = os.path.join(HERE, "data", "external_evidence_inventory.csv")
LAMMPS_RESULTS = os.path.join(HERE, "reference", "lammps", "lammps_results.csv")
PLOTS = os.path.join(HERE, "plots")
REQUIRED = ("axial_strain", "f_h_mean", "f_v_mean", "wall_force_ratio", "contacts")
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
        "positive_contacts": max(row["contacts"] for row in rows) > 0.0,
        "positive_platen_reaction": True,
        "ratio_recomputed": ratio_recomputed,
        "max_ratio_residual": max_ratio_residual,
        "ratio_min": min(row["wall_force_ratio"] for row in loaded),
        "ratio_max": max(row["wall_force_ratio"] for row in loaded),
    }
    return all(checks[key] for key in ("forward_compression", "positive_contacts",
                                        "positive_platen_reaction", "ratio_recomputed")), checks


def _interpolate(points, x):
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    raise RuntimeError("comparison strain lies outside a trajectory")


def compare_independent_lammps(dirt_rows, path=LAMMPS_RESULTS):
    """Report a solver-to-solver response comparison without a fitted gate.

    The source paper does not provide a trajectory.  LAMMPS is therefore an
    independently executable *equivalent-protocol* reference, not a substitute
    primary-source curve.  There is intentionally no tolerance here: this
    function reports measured disagreement, which reviewers can inspect along
    with the raw external trajectory.
    """
    lammps = read(path)
    dirt = [(r["axial_strain"], r["f_v_mean"]) for r in dirt_rows
            if r["f_v_mean"] > 0.0]
    ext = [(r["axial_strain"], r["syy"]) for r in lammps if r["syy"] > 0.0]
    grid = [0.01 + 0.005 * i for i in range(12)]
    if dirt[0][0] > grid[0] or dirt[-1][0] < grid[-1] or ext[0][0] > grid[0] or ext[-1][0] < grid[-1]:
        raise RuntimeError("DIRT and LAMMPS trajectories must span 0.01 <= strain <= 0.065")
    d = [_interpolate(dirt, x) for x in grid]
    l = [_interpolate(ext, x) for x in grid]
    d = [v / d[0] for v in d]
    l = [v / l[0] for v in l]
    d_mean, l_mean = sum(d) / len(d), sum(l) / len(l)
    denominator = math.sqrt(sum((v - d_mean) ** 2 for v in d) * sum((v - l_mean) ** 2 for v in l))
    correlation = sum((a - d_mean) * (b - l_mean) for a, b in zip(d, l)) / denominator
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(d, l)) / len(d))
    return {"samples": len(grid), "normalized_axial_correlation": correlation,
            "normalized_axial_rmse": rmse, "dirt_normalized": d, "lammps_normalized": l}


def plot(rows, checks, passed, source, comparison):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f"(matplotlib unavailable; plot not regenerated: {exc})")
        return

    os.makedirs(PLOTS, exist_ok=True)
    x = [row["axial_strain"] for row in rows]
    lammps = read(LAMMPS_RESULTS)
    lammps_x = [row["axial_strain"] for row in lammps if row["syy"] > 0.0]
    lammps_y_raw = [row["syy"] for row in lammps if row["syy"] > 0.0]
    dirt_vertical = [row["f_v_mean"] for row in rows if row["f_v_mean"] > 0.0]
    dirt_vertical_x = [row["axial_strain"] for row in rows if row["f_v_mean"] > 0.0]

    fig, (cross_code, state, ratio) = plt.subplots(3, 1, figsize=(7.2, 9.4), sharex=True)
    # This is the decisive falsification panel.  The normalization is fixed at
    # the first common strain sample in compare_independent_lammps(), not fit
    # to improve agreement.  LAMMPS is an analogue, so disagreement is shown
    # rather than turned into a source-replication verdict.
    cross_code.plot(dirt_vertical_x, [value / dirt_vertical[0] for value in dirt_vertical],
                    "o-", ms=2.5, label="DIRT platen reaction / first positive reaction")
    cross_code.plot(lammps_x, [value / lammps_y_raw[0] for value in lammps_y_raw],
                    "-", lw=1.4, label="independent LAMMPS $\\sigma_{yy}$ / initial")
    cross_code.axvspan(0.01, 0.065, color="0.85", alpha=.55,
                       label="fixed comparison interval")
    cross_code.set_ylabel("normalized axial response")
    cross_code.set_title(
        "Independent analogue: NOT a source replication "
        f"($r={comparison['normalized_axial_correlation']:.3f}$, "
        f"NRMSE={comparison['normalized_axial_rmse']:.3f}$)"
    )
    cross_code.legend(fontsize=7, loc="best")
    cross_code.grid(alpha=.25)

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
    for stage, value in sorted(source.items()):
        ratio.axhline(value, color="tab:orange", ls="--", lw=1.2,
                     label="Cundall--Strack Fig. 10 (unregistered)" if stage == "A" else None)
        ratio.annotate(f"Fig. 10 {stage}: {value:.2f}", (x[-1], value),
                       xytext=(-5, 4), textcoords="offset points", ha="right", fontsize=7)
    ratio.set_ylabel(r"$F_H/F_V$")
    ratio.set_title("Wall measurement integrity: " + ("PASS" if passed else "FAIL") +
                    "; source stages are not state-registered")
    ratio.legend(fontsize=7)
    ratio.grid(alpha=.25)
    ratio.set_xlabel("axial platen strain from first recorder state")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "stress_volume_response.png"), dpi=180)
    plt.close(fig)


def main():
    command = sys.argv[1] if len(sys.argv) > 1 else "all"
    if command == "external":
        rows = read(RESULTS)
        comparison = compare_independent_lammps(rows)
        print("independent LAMMPS comparison: " + "; ".join(
            f"{key}={value}" for key, value in comparison.items()
            if key not in {"dirt_normalized", "lammps_normalized"}))
        # These observed values are deliberately not transformed into a PASS
        # band.  The two solvers disagree, and the primary source lacks the
        # trajectory needed to adjudicate a claimed Cundall--Strack replication.
        raise SystemExit("external comparison reports non-reproduction; no source-derived acceptance gate exists")
    if command in ("all", "run"):
        build_run()
    if command in ("all", "graph"):
        rows = read(RESULTS)
        source = read_reference()
        protocol = read_protocol()
        missing = audit_external_evidence()
        passed, checks = evaluate(rows)
        comparison = compare_independent_lammps(rows)
        plot(rows, checks, passed, source, comparison)
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
        print("independent LAMMPS analogue: " + "; ".join(
            f"{key}={value}" for key, value in comparison.items()
            if key not in {"dirt_normalized", "lammps_normalized"}))
    if command not in ("all", "run", "graph", "external"):
        raise SystemExit(f"unknown command {command!r}; use all, run, graph, or external")


if __name__ == "__main__":
    main()
