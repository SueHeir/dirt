#!/usr/bin/env python3
"""Independent single-contact Coulomb diagnostic for the repose study.

This is deliberately separate from ``sweep.py``: it reads the declared TOML
material, rather than sharing the sweep constants or its result parser. It is
an analytical diagnostic, not a fitted angle estimator and not an experimental
reference.

For an exposed non-cohesive grain on a planar free surface, force balance gives
``m g sin(theta) <= mu_p m g cos(theta)``, hence ``mu_p >= tan(theta)``.
The LAMMPS/DIRT SDS rolling term is a *couple* (``r_eff n x F_roll``); it can
affect rotation, but supplies no net tangential force to this one-contact
balance. A real heap has a contact network and geometric interlocking, so this
calculation is explicitly **not** a global upper bound or a campaign gate.
"""

import argparse
import math
from pathlib import Path
import tomllib


DEFAULT_CONFIG = Path(__file__).with_name("config.toml")
DEFAULT_BAND_DEG = (22.0, 26.0)


def declared_sliding_friction(config_path: Path) -> float:
    """Read the one mobile material's Coulomb sliding coefficient from TOML."""
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    materials = config.get("dem", {}).get("materials", [])
    if len(materials) != 1:
        raise ValueError("expected exactly one declared DEM material")
    model = config.get("dem", {}).get("rolling_model")
    if model != "sds":
        raise ValueError(f"expected SDS rolling model, found {model!r}")
    mu_p = materials[0].get("friction")
    if not isinstance(mu_p, (int, float)) or not math.isfinite(mu_p) or mu_p < 0:
        raise ValueError("material friction must be a finite non-negative number")
    return float(mu_p)


def single_contact_sliding_coefficient(angle_deg: float) -> float:
    """Coulomb coefficient required by a one-contact inclined-plane balance."""
    if not 0.0 <= angle_deg < 90.0:
        raise ValueError("angle must be in [0, 90) degrees")
    return math.tan(math.radians(angle_deg))


def single_contact_supports_band(mu_p: float, band_deg: tuple[float, float]) -> bool:
    """Whether a one-contact Coulomb control would support the band lower edge."""
    lo, hi = band_deg
    if not 0.0 <= lo <= hi < 90.0:
        raise ValueError("invalid target band")
    return mu_p >= single_contact_sliding_coefficient(lo)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--band", type=float, nargs=2, metavar=("LOW_DEG", "HIGH_DEG"),
                        default=DEFAULT_BAND_DEG)
    args = parser.parse_args()
    try:
        mu_p = declared_sliding_friction(args.config)
        band = tuple(args.band)
        required = single_contact_sliding_coefficient(band[0])
        possible = single_contact_supports_band(mu_p, band)
    except (OSError, ValueError, tomllib.TOMLDecodeError) as error:
        print(f"RESULT: ERROR ({error})")
        return 2
    print("=== Single-contact Coulomb diagnostic ===")
    print(f"  config: {args.config}")
    print(f"  declared Coulomb sliding coefficient: mu_p = {mu_p:.6g}")
    print(f"  target band: [{band[0]:.3g}, {band[1]:.3g}] deg")
    print(f"  one-contact mu_p at {band[0]:.3g} deg: tan(theta) = {required:.6g}")
    print("  SDS is a couple, not a tangential-force allowance in this control")
    if possible:
        print("RESULT: CONTROL SUPPORTS BAND (not a heap prediction or calibration)")
        return 0
    print("RESULT: CONTROL DOES NOT SUPPORT BAND (contact-network effects remain "
          "possible; campaign and external checks are still required)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
