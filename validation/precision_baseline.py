#!/usr/bin/env python3
"""CPU precision-validation baseline.

Runs the example benchmarks on the CPU under each host-storage precision
(double / mixed / single) and records a deterministic fingerprint of each run's
output. This is the reference baseline the GPU (always f32 ≈ mixed/single) is
validated against later: a GPU run of the same config should reproduce the
CPU-single / CPU-mixed fingerprint within f32 tolerance.

Usage:
    python3 validation/precision_baseline.py            # run the default set
    python3 validation/precision_baseline.py ex1 ex2    # run specific examples

Output:
    validation/results/<example>__<precision>.csv   raw output (archived)
    validation/cpu_precision_baseline.csv            machine-readable table
    validation/cpu_precision_baseline.md             human-readable summary
"""
import csv
import glob
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "validation", "results")
PRECISIONS = ["precision-double", "precision-mixed", "precision-single"]

# Contact-physics benchmarks first: fast, deterministic, and they exercise the
# exact normal/tangential contact force the GPU kernels compute.
DEFAULT_EXAMPLES = [
    "bench_hertz_rebound",
    "bench_oblique_impact",
    "bench_rolling_decay",
    "bench_sliding_friction",
    "bench_sphere_haff_cooling",
    "bench_clump_haff_cooling",
    "bench_rod_haff_cooling",
    "bench_jkr_adhesion",
]

RUN_TIMEOUT = 600  # seconds per run


def build(example, precision):
    subprocess.run(
        ["cargo", "build", "-q", "--release", "--example", example,
         "--no-default-features", "--features", precision],
        cwd=REPO, check=True,
    )


def newest_csv(example):
    """Newest CSV the example wrote under its data/ dir."""
    pat = os.path.join(REPO, "examples", example, "data", "**", "*.csv")
    files = glob.glob(pat, recursive=True)
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


def run_one(example, precision):
    config = os.path.join(REPO, "examples", example, "config.toml")
    # Run from REPO so the example writes into examples/<name>/data/.
    t0 = time.time()
    try:
        proc = subprocess.run(
            [os.path.join(REPO, "target", "release", "examples", example), config],
            cwd=REPO, capture_output=True, text=True, timeout=RUN_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "secs": RUN_TIMEOUT}
    dt = time.time() - t0
    if proc.returncode != 0:
        return {"status": f"exit {proc.returncode}", "secs": dt,
                "stderr": proc.stderr[-400:]}
    out = newest_csv(example)
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
    return fp


def main():
    examples = sys.argv[1:] or DEFAULT_EXAMPLES
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
                print(f"  {ex}: BUILD FAILED")
                continue
            r = run_one(ex, precision)
            results[ex][precision] = r
            if r["status"] == "ok":
                print(f"  {ex}: ok  rows={r['rows']} sig={r['sig']:.10g} "
                      f"({r['secs']:.1f}s)")
            else:
                print(f"  {ex}: {r['status']} ({r.get('secs', 0):.1f}s)")

    # ---- machine-readable table ----
    csv_path = os.path.join(REPO, "validation", "cpu_precision_baseline.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
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
        "| example | double signature | mixed Δ vs double | single Δ vs double | rows |",
        "|---|---|---|---|---|",
    ]
    for ex in examples:
        d = results[ex].get("precision-double", {})
        m = results[ex].get("precision-mixed", {})
        s = results[ex].get("precision-single", {})
        if d.get("status") != "ok":
            md.append(f"| {ex} | {d.get('status', '—')} | — | — | — |")
            continue
        dsig = d["sig"]

        def rel(o):
            if o.get("status") != "ok" or dsig == 0:
                return o.get("status", "—")
            return f"{abs(o['sig'] - dsig) / abs(dsig):.2e}"
        md.append(f"| {ex} | {dsig:.10g} | {rel(m)} | {rel(s)} | {d['rows']} |")
    md_path = os.path.join(REPO, "validation", "cpu_precision_baseline.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md) + "\n")

    print(f"\nWrote {csv_path}\n      {md_path}\n      {RESULTS}/")


if __name__ == "__main__":
    main()
