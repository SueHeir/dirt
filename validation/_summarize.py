#!/usr/bin/env python3
"""Regenerate the precision baseline summary from archived run artifacts.

This is intentionally archive-first: completed runs are read from
validation/results/*.csv and non-OK runs from validation/results/*.status.json.
That keeps validation/cpu_precision_baseline.{csv,md},
validation/cpu_precision_final_states.json, and the precision plot reproducible
without re-running the long bulk simulations.
"""
import csv
import glob
import json
import os
import sys

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    bench_python = os.path.expanduser("~/.venvs/bench/bin/python3")
    if os.path.exists(bench_python) and os.path.abspath(sys.executable) != bench_python:
        os.execv(bench_python, [bench_python, *sys.argv])
    raise

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "validation", "results")
PLOTS = os.path.join(REPO, "validation", "plots")
PRECS = ["precision-double", "precision-mixed", "precision-single"]
CONTACT_EXAMPLES = [
    "bench_hertz_rebound",
    "bench_oblique_impact",
    "bench_rolling_decay",
    "bench_sliding_friction",
    "bench_sphere_haff_cooling",
    "bench_clump_haff_cooling",
    "bench_rod_haff_cooling",
    "bench_jkr_adhesion",
]
BULK_EXAMPLES = [
    "bench_angle_of_repose",
    "bench_column_collapse",
    "bench_hopper_beverloo",
    "bench_granular_conductivity",
    "bench_fiber_crossover",
    "bench_lebc_shear",
    "bench_plate_sinkage",
]


def fingerprint(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        return None
    hdr, body = rows[0], rows[1:]
    sig = 0.0
    for r in body:
        for c in r:
            try:
                sig += abs(float(c))
            except ValueError:
                pass
    return {
        "status": "ok",
        "rows": len(body),
        "sig": sig,
        "last": dict(zip(hdr, body[-1])),
    }


def precision_from_name(path, suffix):
    stem = os.path.basename(path)[: -len(suffix)]
    ex, p = stem.rsplit("__", 1)
    return ex, p


def load_archive():
    data = {}
    for f in sorted(glob.glob(os.path.join(RESULTS, "*.csv"))):
        ex, p = precision_from_name(f, ".csv")
        r = fingerprint(f)
        if r:
            data.setdefault(ex, {})[p] = r
    for f in sorted(glob.glob(os.path.join(RESULTS, "*.status.json"))):
        ex, p = precision_from_name(f, ".status.json")
        with open(f) as src:
            status = json.load(src)
        data.setdefault(ex, {})[p] = status
    return data


def ordered_examples(data):
    known = CONTACT_EXAMPLES + BULK_EXAMPLES
    extra = sorted(ex for ex in data if ex not in known)
    return [ex for ex in known if ex in data] + extra


def rel_delta(record, double_sig):
    if not record or record.get("status") != "ok" or not double_sig:
        return None
    return abs(record["sig"] - double_sig) / abs(double_sig)


def write_final_states(data, examples):
    final_states = {}
    for ex in examples:
        final_states[ex] = {}
        for p in PRECS:
            r = data.get(ex, {}).get(p)
            if not r:
                continue
            if r.get("status") == "ok":
                final_states[ex][p] = r.get("last", {})
            else:
                final_states[ex][p] = {
                    "status": r.get("status", "unknown"),
                    "secs": r.get("secs"),
                }
    with open(os.path.join(REPO, "validation", "cpu_precision_final_states.json"), "w") as f:
        json.dump(final_states, f, indent=2, sort_keys=True)
        f.write("\n")


def write_csv(data, examples):
    with open(os.path.join(REPO, "validation", "cpu_precision_baseline.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["example", "precision", "status", "rows", "signature_sum_abs"])
        for ex in examples:
            for p in PRECS:
                r = data.get(ex, {}).get(p, {})
                w.writerow(
                    [
                        ex,
                        p,
                        r.get("status", "missing"),
                        r.get("rows", ""),
                        repr(r.get("sig", "")),
                    ]
                )


def section_rows(data, examples):
    rows = []
    for ex in examples:
        d = data.get(ex, {}).get("precision-double", {})
        m = data.get(ex, {}).get("precision-mixed", {})
        s = data.get(ex, {}).get("precision-single", {})
        if d.get("status") != "ok":
            rows.append(f"| {ex} | {d.get('status', 'missing')} | {m.get('status', '-')} | {s.get('status', '-')} | - |")
            continue
        dsig = d["sig"]

        def fmt(o):
            delta = rel_delta(o, dsig)
            return o.get("status", "-") if delta is None else f"{delta:.2e}"

        rows.append(f"| {ex} | {dsig:.10g} | {fmt(m)} | {fmt(s)} | {d['rows']} |")
    return rows


def write_markdown(data):
    md = [
        "# CPU precision-validation baseline",
        "",
        "Deterministic fingerprint of each example output under each host-storage precision.",
        "`signature` = sum of |numeric cells| in the output CSV; `Delta vs double` is",
        "the relative difference of that signature from the double-precision run.",
        "Mixed/single store positions as f32, so they bound what the f32 GPU should reproduce.",
        "Raw completed outputs and non-OK status JSON files are archived under `validation/results/`.",
        "",
        "Default contact run: `python3 validation/precision_baseline.py`.",
        "Bulk/steady-state long run: `python3 validation/precision_baseline.py --set bulk --timeout 1200`.",
        "Combined archive regeneration: `python3 validation/_summarize.py`.",
        "",
        "![CPU precision deltas](plots/cpu_precision_deltas.png)",
        "",
        "*Mixed/single precision signature deltas relative to double. Completed runs are plotted;",
        "the dashed 10% line is a visible large-drift reference, and non-OK runs are listed below",
        "instead of being silently dropped.*",
        "",
        "## Contact physics",
        "",
        "| example | double signature | mixed Delta vs double | single Delta vs double | rows |",
        "|---|---|---|---|---|",
        *section_rows(data, [ex for ex in CONTACT_EXAMPLES if ex in data]),
        "",
        "## Bulk and steady state",
        "",
        "| example | double signature | mixed Delta vs double | single Delta vs double | rows |",
        "|---|---|---|---|---|",
        *section_rows(data, [ex for ex in BULK_EXAMPLES if ex in data]),
    ]
    non_ok = []
    for ex in ordered_examples(data):
        for p in PRECS:
            r = data.get(ex, {}).get(p, {})
            if r.get("status") and r.get("status") != "ok":
                secs = r.get("secs")
                suffix = f" after {secs:g} s" if isinstance(secs, (int, float)) else ""
                non_ok.append(f"- `{ex}` `{p}`: {r['status']}{suffix}.")
    if non_ok:
        md += ["", "## Non-OK statuses", "", *non_ok]
    with open(os.path.join(REPO, "validation", "cpu_precision_baseline.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def write_plot(data, examples):
    completed = []
    mixed = []
    single = []
    non_ok = []
    for ex in examples:
        d = data.get(ex, {}).get("precision-double", {})
        if d.get("status") != "ok":
            non_ok.append(ex)
            continue
        dm = rel_delta(data.get(ex, {}).get("precision-mixed", {}), d["sig"])
        ds = rel_delta(data.get(ex, {}).get("precision-single", {}), d["sig"])
        if dm is None or ds is None:
            non_ok.append(ex)
            continue
        completed.append(ex)
        mixed.append(max(dm, 1e-13))
        single.append(max(ds, 1e-13))

    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 5.6))
    x = list(range(len(completed)))
    width = 0.38
    ax.bar([i - width / 2 for i in x], mixed, width, label="mixed vs double", color="#31688e")
    ax.bar([i + width / 2 for i in x], single, width, label="single vs double", color="#35b779")
    ax.axhline(1e-1, color="#b23a48", linestyle="--", linewidth=1.2, label="10% large-drift reference")
    ax.set_yscale("log")
    ax.set_ylim(1e-13, max(2e-1, max(mixed + single) * 2 if mixed else 2e-1))
    ax.set_ylabel("relative signature delta")
    ax.set_title("CPU precision baseline: mixed/single deltas from double")
    ax.set_xticks(x)
    ax.set_xticklabels([ex.replace("bench_", "") for ex in completed], rotation=45, ha="right")
    ax.grid(axis="y", which="both", alpha=0.25)
    ax.legend(loc="upper left")
    if non_ok:
        ax.text(
            0.99,
            0.95,
            "Non-OK, no fingerprint: " + ", ".join(non_ok),
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#777777", "alpha": 0.9},
        )
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "cpu_precision_deltas.png"), dpi=180)
    plt.close(fig)


def main():
    data = load_archive()
    examples = ordered_examples(data)
    write_final_states(data, examples)
    write_csv(data, examples)
    write_plot(data, examples)
    write_markdown(data)
    print(f"merged {len(examples)} examples into the precision baseline summary")


if __name__ == "__main__":
    main()
