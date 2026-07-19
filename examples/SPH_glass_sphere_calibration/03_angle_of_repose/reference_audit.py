#!/usr/bin/env python3
"""Audit bibliographic provenance separately from experimental comparability."""

import json
import math
import re
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

REQUIRED = ("doi", "title", "first_author_family", "published_year",
            "source_url", "comparability")
MATCH_FIELDS = ("sphere_radius_m", "density_kg_m3", "restitution",
                "sliding_friction", "container", "deposition", "release",
                "angle_definition", "band_deg")
MATCH_EVIDENCE_FIELDS = ("source_locator", "extraction_method", "observations")
OBSERVATION_FIELDS = ("value_deg", "source_locator")
REPORTED_BAND_FIELDS = ("reported_band_deg", "reported_band_locator")


def _normal(value):
    return re.sub(r"\s+", " ", str(value).casefold()).strip()


def _doi_url(doi):
    """Return the one durable bibliographic locator admitted by this audit.

    A record is evidence *about* a cited source, rather than a convenient URL
    chosen by the campaign.  Requiring the canonical DOI resolver prevents a
    locally-authored JSON record from attaching otherwise correct Crossref
    metadata to an unrelated page.
    """
    return "https://doi.org/" + doi


def _primary_source_locator(locator, doi, field):
    """Require a locator within the Crossref-checked primary DOI source.

    A local campaign artifact or an unanchored table label cannot establish an
    external measurement. This is a provenance condition, not a numerical gate.
    """
    if not isinstance(locator, str) or not locator.strip():
        raise ValueError(f"{field} must be text")
    if not locator.strip().startswith(_doi_url(doi) + "#"):
        raise ValueError(f"{field} must locate content in the primary DOI source")


def load_record(path):
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    missing = [field for field in REQUIRED if not record.get(field)]
    if missing:
        raise ValueError("missing required fields: " + ", ".join(missing))
    doi = str(record["doi"]).lower().strip().removeprefix("https://doi.org/")
    if not doi.startswith("10.") or any(char.isspace() for char in doi):
        raise ValueError("DOI must be a bare DOI beginning with 10.")
    if _normal(record["source_url"]) != _normal(_doi_url(doi)):
        raise ValueError("source_url must be the canonical HTTPS DOI resolver")
    if record["comparability"] not in ("incompatible", "not_assessed", "matched"):
        raise ValueError("invalid comparability")
    if record["comparability"] != "matched" and record.get("observations"):
        raise ValueError("a non-matched record must not contain numerical observations")
    if record["comparability"] == "matched":
        missing_match = [field for field in MATCH_FIELDS if field not in record]
        if missing_match:
            raise ValueError("matched record lacks: " + ", ".join(missing_match))
        missing_evidence = [field for field in MATCH_EVIDENCE_FIELDS
                            if not record.get(field)]
        if missing_evidence:
            raise ValueError("matched record lacks human-auditable evidence: "
                             + ", ".join(missing_evidence))
        if not isinstance(record["observations"], list):
            raise ValueError("matched record observations must be a list")
        if not record["observations"]:
            raise ValueError("matched record must contain numerical observations")
        values = []
        for index, observation in enumerate(record["observations"]):
            if not isinstance(observation, dict):
                raise ValueError(f"observation {index} must be an object")
            missing_observation = [field for field in OBSERVATION_FIELDS
                                   if field not in observation]
            if missing_observation:
                raise ValueError(f"observation {index} lacks: "
                                 + ", ".join(missing_observation))
            try:
                value = float(observation["value_deg"])
            except (TypeError, ValueError) as error:
                raise ValueError(f"observation {index} value_deg must be numeric") from error
            if not math.isfinite(value):
                raise ValueError(f"observation {index} value_deg must be finite")
            _primary_source_locator(observation["source_locator"], doi,
                                    f"observation {index} source_locator")
            values.append(value)
        # The reported comparison band must be recoverable from the cited raw
        # values.  This prevents a record from attaching arbitrary hand-entered
        # bounds to otherwise genuine bibliographic metadata.
        band = record["band_deg"]
        if (not isinstance(band, list) or len(band) != 2
                or not all(isinstance(value, (int, float)) and math.isfinite(value)
                           for value in band) or band[0] >= band[1]):
            raise ValueError("matched record band_deg must be two ordered finite numbers")
        if min(values) < band[0] or max(values) > band[1]:
            raise ValueError("matched record band_deg does not bound its observations")
        # `band_deg` is the value the executable compares against.  It must not
        # be an interval constructed locally around a set of source values: a
        # matched record must identify the interval as reported by the primary
        # source, with its own locator.  This does not prove that a transcription
        # is right (that remains a human review task), but it prevents the
        # campaign target from being silently manufactured by this JSON file.
        missing_reported_band = [field for field in REPORTED_BAND_FIELDS
                                 if field not in record]
        if missing_reported_band:
            raise ValueError("matched record lacks source-reported band: "
                             + ", ".join(missing_reported_band))
        reported_band = record["reported_band_deg"]
        if (not isinstance(reported_band, list) or len(reported_band) != 2
                or not all(isinstance(value, (int, float)) and math.isfinite(value)
                           for value in reported_band)
                or reported_band[0] >= reported_band[1]):
            raise ValueError("reported_band_deg must be two ordered finite numbers")
        if reported_band != band:
            raise ValueError("band_deg must equal the primary source reported_band_deg")
        _primary_source_locator(record["reported_band_locator"], doi,
                                "reported_band_locator")
        # A page/figure/table locator and a declared extraction method make the
        # source claim reviewable.  They are provenance fields, not a numerical
        # criterion inferred from the DIRT campaign.
        _primary_source_locator(record["source_locator"], doi, "source_locator")
        if not isinstance(record["extraction_method"], str):
            raise ValueError("matched record extraction_method must be text")
    record["doi"] = doi
    return record


def audit_metadata(record, opener=urlopen):
    request = Request("https://api.crossref.org/works/" + quote(record["doi"], safe="/"),
                      headers={"Accept": "application/json", "User-Agent": "dirt-reference-audit/1.0"})
    with opener(request, timeout=15) as response:
        message = json.load(response).get("message", {})
    returned_doi = str(message.get("DOI", "")).lower().strip()
    if returned_doi != record["doi"]:
        raise ValueError("Crossref returned a different DOI")
    title = (message.get("title") or [""])[0]
    author = (message.get("author") or [{}])[0].get("family", "")
    dates = (message.get("published-print") or message.get("published-online")
             or message.get("issued") or {}).get("date-parts", [[]])
    year = dates[0][0] if dates and dates[0] else None
    checks = (("title", record["title"], title),
              ("first_author_family", record["first_author_family"], author),
              ("published_year", record["published_year"], year))
    bad = [field for field, expected, actual in checks if _normal(expected) != _normal(actual)]
    if bad:
        raise ValueError("Crossref metadata disagrees: " + ", ".join(bad))


def protocol_matches(record, expected):
    if record["comparability"] != "matched":
        return False
    if any(record.get(field) != value for field, value in expected.items()):
        return False
    band = record["band_deg"]
    return isinstance(band, list) and len(band) == 2 and band[0] < band[1]


def audit_record(path, opener=urlopen):
    """Load a local record then verify its bibliography against Crossref.

    Crossref validates citation identity only.  It does not attest that an
    extracted angle, protocol match, or material transfer is scientifically
    valid; those remain an expert-review responsibility.
    """
    record = load_record(path)
    audit_metadata(record, opener=opener)
    return record


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("record")
    try:
        record = audit_record(parser.parse_args().record)
        print("EXTERNAL RECORD METADATA VERIFIED (Crossref only)")
        print("NOT VERIFIED: source extraction, protocol equivalence, material transfer, or calibration.")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit("REFERENCE AUDIT FAILED: " + str(error))
