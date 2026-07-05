#!/usr/bin/env python3
"""CPU precision-validation baseline.

Runs the example benchmarks on the CPU under each host-storage precision
(double / mixed / single) and records a deterministic fingerprint of each run's
output. This is the reference baseline the GPU (always f32 ≈ mixed/single) is
validated against later: a GPU run of the same config should reproduce the
CPU-single / CPU-mixed fingerprint within f32 tolerance.

Usage:
    python3 validation/precision_baseline.py               # run the contact set
    python3 validation/precision_baseline.py --set bulk    # run long bulk/steady-state set
    python3 validation/precision_baseline.py --set all     # run contact + bulk sets
    python3 validation/precision_baseline.py ex1 ex2       # run specific examples

Output:
    validation/results/<example>__<precision>.csv   raw output (archived)
    validation/cpu_precision_baseline.csv            machine-readable table
    validation/cpu_precision_baseline.md             human-readable summary
"""
import csv
import glob
import json
import os
import argparse
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "validation", "results")
PRECISIONS = ["precision-double", "precision-mixed", "precision-single"]

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

EXAMPLE_SETS = {
    "contact": CONTACT_EXAMPLES,
    "bulk": BULK_EXAMPLES,
    "all": CONTACT_EXAMPLES + BULK_EXAMPLES,
}

DEFAULT_TIMEOUT = 1200  # seconds per run; bulk examples can take several minutes


def build(example, precision):
    subprocess.run(
        ["cargo", "build", "-q", "--release", "--example", example,
         "--no-default-features", "--features", precision],
        cwd=REPO, check=True,
    )


def newest_csv(example, min_mtime=0.0):
    """Newest CSV the example wrote under its data/ dir after min_mtime."""
    pat = os.path.join(REPO, "examples", example, "data", "**", "*.csv")
    files = [f for f in glob.glob(pat, recursive=True)
             if os.path.getmtime(f) >= min_mtime]
    return max(files, key=os.path.getmtime) if files else None


def fingerprint(csv_path):
    """Deterministic reduction: row count, sum-of-abs of all numeric cells, and
    the final data row. Sum-of-abs is one comparable scalar; the final row is
    the physical end-state. Both shift predictably with storage precision."""
    rows = []
    with open(csv_path) as f:
        reader = csv.reader(f)
        header = next(reader, [])
        for r in reader:
            rows.append(r)
    sig = 0.0
    for r in rows:
        for cell in r:
            try:
                sig += abs(float(cell))
            except ValueError:
                pass
    last = dict(zip(header, rows[-1])) if rows else {}
    return {"rows": len(rows), "sig": sig, "last": last, "header": header}


def status_path(example, precision):
    return os.path.join(RESULTS, f"{example}__{precision}.status.json")


def archive_status(example, precision, result):
    os.makedirs(RESULTS, exist_ok=True)
    with open(status_path(example, precision), "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")


def clear_status(example, precision):
    try:
        os.remove(status_path(example, precision))
    except FileNotFoundError:
        pass


def run_one(example, precision, timeout):
    config = os.path.join(REPO, "examples", example, "config.toml")
    # Run from REPO so the example writes into examples/<name>/data/.
    t0 = time.time()
    try:
        proc = subprocess.run(
            [os.path.join(REPO, "target", "release", "examples", example), config],
            cwd=REPO, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "secs": timeout}
    dt = time.time() - t0
    if proc.returncode != 0:
        return {"status": f"exit {proc.returncode}", "secs": dt,
                "stderr": proc.stderr[-400:]}
    out = newest_csv(example, t0)
    if not out:
        return {"status": "no-output", "secs": dt}
    fp = fingerprint(out)
    # Archive the raw output for later GPU diffing.
    os.makedirs(RESULTS, exist_ok=True)
    dst = os.path.join(RESULTS, f"{example}__{precision}.csv")
    with open(out) as s, open(dst, "w") as d:
        d.write(s.read())
    fp["status"] = "ok"
    fp["secs"] = dt
    clear_status(example, precision)
    return fp


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--set",
        choices=sorted(EXAMPLE_SETS),
        default="contact",
        help="named benchmark set to run; use 'bulk' for the long steady-state set",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help="per-example timeout in seconds",
    )
    parser.add_argument(
        "examples",
        nargs="*",
        help="explicit example names; overrides --set",
    )
    return parser.parse_args(argv)


def main():
    args = parse_args(sys.argv[1:])
    examples = args.examples or EXAMPLE_SETS[args.set]
    os.makedirs(RESULTS, exist_ok=True)
    # results[example][precision] = fingerprint dict
    results = {ex: {} for ex in examples}

    for precision in PRECISIONS:
        print(f"\n===== {precision} =====", flush=True)
        for ex in examples:
            try:
                build(ex, precision)
            except subprocess.CalledProcessError:
                results[ex][precision] = {"status": "build-fail"}
                archive_status(ex, precision, results[ex][precision])
                print(f"  {ex}: BUILD FAILED")
                continue
            r = run_one(ex, precision, args.timeout)
            results[ex][precision] = r
            if r["status"] == "ok":
                print(f"  {ex}: ok  rows={r['rows']} sig={r['sig']:.10g} "
                      f"({r['secs']:.1f}s)")
            else:
                archive_status(ex, precision, r)
                print(f"  {ex}: {r['status']} ({r.get('secs', 0):.1f}s)")

    # ---- machine-readable table ----
    csv_path = os.path.join(REPO, "validation", "cpu_precision_baseline.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["example", "precision", "status", "rows", "signature_sum_abs"])
        for ex in examples:
            for p in PRECISIONS:
                r = results[ex].get(p, {})
                w.writerow([ex, p, r.get("status", "—"),
                            r.get("rows", ""), repr(r.get("sig", ""))])

    # ---- human-readable summary with cross-precision agreement ----
    md = [
        "# CPU precision-validation baseline",
        "",
        "Deterministic fingerprint of each example's output under each host-storage",
        "precision. `signature` = sum of |numeric cells| in the output CSV; `Δ vs double`",
        "is the relative difference of that signature from the double-precision run.",
        "Mixed/single store positions as f32, so they bound what the f32 GPU should",
        "reproduce. Raw outputs archived under `validation/results/`.",
        "",
        "Benchmark sets: default/contact (`python3 validation/precision_baseline.py`),",
        "bulk/steady-state (`python3 validation/precision_baseline.py --set bulk`),",
        "or combined (`python3 validation/precision_baseline.py --set all`).",
        "",
        "| example | double signature | mixed Δ vs double | single Δ vs double | rows |",
        "|---|---|---|---|---|",
    ]
    for ex in examples:
        d = results[ex].get("precision-double", {})
        m = results[ex].get("precision-mixed", {})
        s = results[ex].get("precision-single", {})
        if d.get("status") != "ok":
            md.append(f"| {ex} | {d.get('status', '-')} | - | - | - |")
            continue
        dsig = d["sig"]

        def rel(o):
            if o.get("status") != "ok" or dsig == 0:
                return o.get("status", "-")
            return f"{abs(o['sig'] - dsig) / abs(dsig):.2e}"
        md.append(f"| {ex} | {dsig:.10g} | {rel(m)} | {rel(s)} | {d['rows']} |")
    md_path = os.path.join(REPO, "validation", "cpu_precision_baseline.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md) + "\n")

    print(f"\nWrote {csv_path}\n      {md_path}\n      {RESULTS}/")


if __name__ == "__main__":
    main()
