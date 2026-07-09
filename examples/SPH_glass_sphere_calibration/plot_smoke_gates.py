#!/usr/bin/env python3
"""Regenerate the SPH glass calibration smoke-gate summary figure.

The measured values are the latest smoke-gate evidence recorded in the
top-level calibration README. This is intentionally only a figure generator:
the PASS/FAIL decisions remain in each subfolder's sweep.py stdout.
"""

from __future__ import annotations

from pathlib import Path


HERE = Path(__file__).resolve().parent

SHEAR = [
    {"phi": 0.19, "mu": 0.444, "lo": 0.15, "hi": 0.90},
    {"phi": 0.29, "mu": 0.490, "lo": 0.15, "hi": 0.90},
    {"phi": 0.38, "mu": 0.483, "lo": 0.15, "hi": 0.90},
]

ENDURING = [
    {"phi": 0.29, "frac": 0.05, "reference": 0.0},
    {"phi": 0.43, "frac": 0.15, "reference": 0.05},
    {"phi": 0.52, "frac": 0.39, "reference": 0.05},
]


def main() -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150})
    fig, (ax_mu, ax_contact) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    phis = [p["phi"] for p in SHEAR]
    ax_mu.plot(phis, [p["mu"] for p in SHEAR], "o-", color="#1f77b4", label="measured")
    ax_mu.plot(phis, [p["lo"] for p in SHEAR], "v", color="#555555", label="criterion low")
    ax_mu.plot(phis, [p["hi"] for p in SHEAR], "^", color="#555555", label="criterion high")
    ax_mu.set_xlabel("measured solid fraction Phi")
    ax_mu.set_ylabel("mu = |sigma_xy| / P")
    ax_mu.set_title("01 shear rheology")
    ax_mu.set_ylim(0.0, 1.02)
    ax_mu.grid(True, alpha=0.25)
    ax_mu.legend(fontsize=9, loc="best")

    ephis = [p["phi"] for p in ENDURING]
    ax_contact.plot(
        ephis,
        [p["frac"] for p in ENDURING],
        "o-",
        color="#d62728",
        label="measured",
    )
    ax_contact.plot(
        ephis,
        [p["reference"] for p in ENDURING],
        "s--",
        color="#555555",
        label="dilute-floor reference",
    )
    ax_contact.set_xlabel("measured solid fraction Phi")
    ax_contact.set_ylabel("sigma_contact / p_DEM")
    ax_contact.set_title("04 enduring contact")
    ax_contact.set_ylim(-0.02, 0.46)
    ax_contact.grid(True, alpha=0.25)
    ax_contact.legend(fontsize=9, loc="best")

    fig.suptitle("SPH glass calibration smoke-gate measurements vs criteria")
    fig.tight_layout()
    out = HERE / "plots" / "smoke_gates.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
