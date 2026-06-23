#!/usr/bin/env python3
"""Regenerate the combined precision baseline summary from ALL archived runs in
validation/results/ (so multiple harness invocations merge into one record)."""
import csv
import glob
import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "validation", "results")
PRECS = ["precision-double", "precision-mixed", "precision-single"]


def fp(path):
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
    return {"rows": len(body), "sig": sig, "last": dict(zip(hdr, body[-1]))}


data = {}
for f in sorted(glob.glob(os.path.join(RESULTS, "*.csv"))):
    ex, p = os.path.basename(f)[:-4].split("__")
    r = fp(f)
    if r:
        data.setdefault(ex, {})[p] = r

# final-states json
fs = {ex: {p: data[ex][p]["last"] for p in data[ex]} for ex in data}
json.dump(fs, open(os.path.join(REPO, "validation", "cpu_precision_final_states.json"), "w"), indent=2)

# csv table
with open(os.path.join(REPO, "validation", "cpu_precision_baseline.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["example", "precision", "rows", "signature_sum_abs"])
    for ex in sorted(data):
        for p in PRECS:
            if p in data[ex]:
                w.writerow([ex, p, data[ex][p]["rows"], repr(data[ex][p]["sig"])])

# md
md = ["# CPU precision-validation baseline", "",
      "Deterministic fingerprint of each example's output under each host-storage",
      "precision. `signature` = sum of |numeric cells|; `Δ vs double` is the relative",
      "difference of that signature from the double run. mixed/single store positions",
      "as f32, bounding what the f32 GPU should reproduce. Raw outputs (gitignored)",
      "regenerate via `python3 validation/precision_baseline.py`.", "",
      "| example | double signature | mixed Δ | single Δ | rows |",
      "|---|---|---|---|---|"]
for ex in sorted(data):
    d = data[ex].get("precision-double")
    if not d:
        md.append(f"| {ex} | (no double run) | — | — | — |")
        continue

    def rel(p):
        o = data[ex].get(p)
        if not o or d["sig"] == 0:
            return "—"
        return f"{abs(o['sig'] - d['sig']) / abs(d['sig']):.2e}"
    md.append(f"| {ex} | {d['sig']:.10g} | {rel('precision-mixed')} | {rel('precision-single')} | {d['rows']} |")
open(os.path.join(REPO, "validation", "cpu_precision_baseline.md"), "w").write("\n".join(md) + "\n")
print(f"merged {len(data)} examples into the baseline summary")
