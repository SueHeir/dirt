#!/usr/bin/env python3
"""Auto-generate the DIRT validation summary dashboard.

Produces a single machine-generated table — `docs/src/reference/validation-summary.md`
— listing every `bench_*` example with its reference, evidence type, checked
quantity, tolerance, latest MEASURED value, and PASS/FAIL, so a scientist
evaluating DIRT can see on one page *what is validated, against what, how
closely, right now*.

The design keeps it honest:

  * The set of benchmarks is discovered from the filesystem (`examples/bench_*`),
    so a benchmark cannot be silently dropped from the dashboard.
  * The static descriptors (reference source, evidence type, checked quantity,
    tolerance) come from the DECLARATIVE manifest `validation/validation_summary.toml`.
    A benchmark with no manifest entry (or a manifest entry with no benchmark) is
    a hard error — the two must stay in lockstep.
  * The dynamic, drift-prone facts — the LATEST MEASURED value and the PASS/FAIL
    verdict — are read from the actual benchmark run outputs (the hourly harness
    at `~/projects/automation`, or any results dir passed via --results-dir /
    $DIRT_BENCH_RESULTS). Nothing about the current numbers is hand-typed, so the
    page cannot drift from what the code actually does.

Usage:
    python3 validation/gen_validation_summary.py            # regenerate the page
    python3 validation/gen_validation_summary.py --check    # fail if out of date
    python3 validation/gen_validation_summary.py --results-dir DIR

Exit status is non-zero if the manifest and the `examples/bench_*` set disagree
(so it works as a CI gate), or, under --check, if the committed page is stale.
"""
import argparse
import glob
import json
import os
import re
import sys
import tomllib
from datetime import datetime, timezone

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(REPO, "validation", "validation_summary.toml")
OUT = os.path.join(REPO, "docs", "src", "reference", "validation-summary.md")

# Where the benchmark run outputs live. The hourly harness (run-bench.sh) writes
# one directory per benchmark id, each with a `last.json` describing the newest
# run: {status, verdict, ts, log, ...}. Overridable so the dashboard is not
# hard-wired to one machine.
DEFAULT_RESULTS = os.environ.get(
    "DIRT_BENCH_RESULTS",
    os.path.expanduser("~/projects/automation/results"),
)

# Verdict strings a run log may print, most-specific first, used only as the
# "latest measured value" fallback when a benchmark declares no `measure` regex.
VERDICT_RE = re.compile(
    r"(ALL CHECKS PASSED|CHECKS? FAILED(?:: \d+ of \d+)?|\d+/\d+ checks passed"
    r"|RESULT: (?:PASS|FAIL))"
)

STATUS_BADGE = {
    "PASS": "✅ PASS",
    "FAIL": "❌ FAIL",
    "TIMEOUT": "⏱ TIMEOUT",
    "NO_DRIVER": "— no driver",
    "MISSING": "— no run",
}


def discover_benches():
    """Every examples/bench_* directory, by id (basename), sorted."""
    dirs = glob.glob(os.path.join(REPO, "examples", "bench_*"))
    return sorted(os.path.basename(d) for d in dirs if os.path.isdir(d))


def load_manifest():
    with open(MANIFEST, "rb") as f:
        data = tomllib.load(f)
    return data.get("bench", {})


def load_latest_run(results_dir, bench_id):
    """Read the newest run record for a benchmark from the harness results dir.

    Returns a dict with status/verdict/ts/log (log text loaded), or None if the
    benchmark has never been run in this results dir.
    """
    last = os.path.join(results_dir, bench_id, "last.json")
    if not os.path.isfile(last):
        return None
    with open(last) as f:
        rec = json.load(f)
    log_text = ""
    log_path = rec.get("log", "")
    if log_path and os.path.isfile(log_path):
        with open(log_path, errors="replace") as f:
            log_text = f.read()
    rec["log_text"] = log_text
    return rec


def measured_value(entry, run):
    """The 'latest measured value' cell — always from the run, never the manifest.

    Prefer the benchmark's declared `measure` regex applied to the newest log;
    otherwise fall back to the run's verdict string. Returns (value, label).
    """
    if run is None:
        return "—", ""
    log_text = run.get("log_text", "")
    pat = entry.get("measure")
    if pat and log_text:
        m = re.search(pat, log_text)
        if m:
            val = m.group(1) if m.groups() else m.group(0)
            return val.strip(), entry.get("measure_label", "")
    # Fallback: the run's own verdict (harness field, else parsed from the log).
    verdict = (run.get("verdict") or "").strip()
    if not verdict and log_text:
        hits = VERDICT_RE.findall(log_text)
        if hits:
            verdict = hits[-1]
    return (verdict or "—"), "verdict"


def fmt_ts(ts):
    """Harness timestamp '20260703T220107Z' -> '2026-07-03 22:01 UTC'."""
    try:
        dt = datetime.strptime(ts, "%Y%m%dT%H%M%SZ").replace(tzinfo=timezone.utc)
        return dt.strftime("%Y-%m-%d %H:%M UTC")
    except (ValueError, TypeError):
        return ts or "—"


def build(results_dir):
    benches = discover_benches()
    manifest = load_manifest()

    # Lockstep check: manifest and filesystem must name the same benchmarks.
    missing_entry = [b for b in benches if b not in manifest]
    orphan_entry = [b for b in manifest if b not in benches]
    if missing_entry or orphan_entry:
        msg = []
        if missing_entry:
            msg.append(
                "benchmarks with no manifest entry (add them to "
                f"validation/validation_summary.toml): {', '.join(missing_entry)}"
            )
        if orphan_entry:
            msg.append(
                "manifest entries with no examples/bench_* directory (remove or "
                f"rename): {', '.join(orphan_entry)}"
            )
        raise SystemExit("ERROR: manifest/benchmark set out of sync:\n  - " +
                         "\n  - ".join(msg))

    rows = []
    newest_ts = ""
    n_pass = n_fail = n_other = n_norun = 0
    for bid in benches:
        entry = manifest[bid]
        run = load_latest_run(results_dir, bid)
        if run is None:
            status = "MISSING"
            n_norun += 1
            ts_disp = "—"
        else:
            status = run.get("status", "?")
            ts = run.get("ts", "")
            newest_ts = max(newest_ts, ts)
            ts_disp = fmt_ts(ts)
            if status == "PASS":
                n_pass += 1
            elif status == "FAIL":
                n_fail += 1
            else:
                n_other += 1

        value, vlabel = measured_value(entry, run)
        value_cell = f"{value}"
        if vlabel and vlabel != "verdict":
            value_cell = f"{value} ({vlabel})"

        rows.append({
            "id": bid,
            "reference": entry.get("reference", "—"),
            "reference_type": entry.get("reference_type", "—"),
            "quantity": entry.get("quantity", "—"),
            "tolerance": entry.get("tolerance", "—"),
            "value": value_cell,
            "status": STATUS_BADGE.get(status, status),
            "ts": ts_disp,
        })

    return rows, {
        "results_dir": results_dir,
        "newest_ts": fmt_ts(newest_ts) if newest_ts else "—",
        "n_pass": n_pass, "n_fail": n_fail, "n_other": n_other,
        "n_norun": n_norun, "n_total": len(benches),
    }


def render(rows, meta):
    L = []
    L.append("# Validation Summary Dashboard")
    L.append("")
    L.append(
        "<!-- AUTO-GENERATED by validation/gen_validation_summary.py — DO NOT EDIT "
        "BY HAND. Regenerate with `python3 validation/gen_validation_summary.py`. -->"
    )
    L.append("")
    L.append(
        "One page: **what every `bench_*` example validates, against what, how "
        "closely, and whether it passes right now.** The *reference*, *evidence "
        "type*, *checked quantity*, and *tolerance* columns are declared in "
        "[`validation/validation_summary.toml`]"
        "(https://github.com/SueHeir/dirt/blob/main/validation/validation_summary.toml); "
        "the **latest measured value** and **status** columns are read straight "
        "from the newest real benchmark run — never hand-typed — so this page "
        "cannot silently drift from what the code does. For the honest, per-bench "
        "narrative (where each test is weak) see "
        "[`examples/VALIDATION.md`]"
        "(https://github.com/SueHeir/dirt/blob/main/examples/VALIDATION.md)."
    )
    L.append("")
    L.append(
        f"**{meta['n_pass']} PASS · {meta['n_fail']} FAIL · "
        f"{meta['n_other']} other · {meta['n_norun']} not-yet-run** "
        f"of {meta['n_total']} benchmarks. "
        f"Newest run: {meta['newest_ts']}. "
        f"Results source: `{meta['results_dir']}`."
    )
    L.append("")
    L.append(
        "| Benchmark | Reference | Evidence type | Checked quantity | Tolerance | "
        "Latest measured value | Status | Last run |"
    )
    L.append("|---|---|---|---|---|---|---|---|")

    def esc(s):
        # Escape pipes (would break the table) and flatten any newlines.
        return str(s).replace("|", "\\|").replace("\n", " ")

    for r in rows:
        L.append(
            f"| `{r['id']}` | {esc(r['reference'])} | {esc(r['reference_type'])} | "
            f"{esc(r['quantity'])} | {esc(r['tolerance'])} | {esc(r['value'])} | "
            f"{esc(r['status'])} | {esc(r['ts'])} |"
        )
    L.append("")
    L.append("## How this page is generated")
    L.append("")
    L.append(
        "`validation/gen_validation_summary.py` discovers every "
        "`examples/bench_*` directory, joins it to its declarative descriptor in "
        "`validation/validation_summary.toml`, and reads the newest run record "
        "(`<results>/<bench>/last.json` + its log) written by the benchmark "
        "harness (`automation/bin/run-bench.sh`). The measured value is pulled "
        "from the run log via the benchmark's declared `measure` regex, falling "
        "back to the run's own PASS/FAIL verdict string. The script errors out if "
        "a benchmark is added or removed without updating the manifest, so a "
        "benchmark can never silently disappear from this dashboard."
    )
    L.append("")
    L.append(
        "Regenerate after a fresh benchmark run with:\n\n"
        "```sh\npython3 validation/gen_validation_summary.py\n```"
    )
    L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", default=DEFAULT_RESULTS,
                    help="benchmark results directory (default: $DIRT_BENCH_RESULTS "
                         "or ~/projects/automation/results)")
    ap.add_argument("--out", default=OUT, help="output markdown path")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit non-zero if the committed page differs")
    args = ap.parse_args()

    rows, meta = build(args.results_dir)
    text = render(rows, meta)

    if args.check:
        current = ""
        if os.path.isfile(args.out):
            with open(args.out) as f:
                current = f.read()
        # Ignore the volatile 'newest run / results source' line under --check:
        # it changes every harness run and would make the gate flap. Compare only
        # the structural table + descriptors.
        def strip_volatile(s):
            return "\n".join(
                ln for ln in s.splitlines()
                if not ln.startswith("**") or "PASS ·" not in ln
            )
        if strip_volatile(current) == strip_volatile(text):
            print("validation-summary.md up to date")
            return 0
        print("validation-summary.md is STALE — run "
              "`python3 validation/gen_validation_summary.py`", file=sys.stderr)
        return 1

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(text)
    print(f"wrote {args.out}")
    print(f"  {meta['n_pass']} PASS, {meta['n_fail']} FAIL, {meta['n_other']} "
          f"other, {meta['n_norun']} not-yet-run of {meta['n_total']} benchmarks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
