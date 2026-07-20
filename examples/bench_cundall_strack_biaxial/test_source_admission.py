"""Regression tests for the primary-source admission boundary.

These tests use a synthetic reader only to exercise the parser.  The checksum
test remains a separate integration step against the archived PDF documented in
the receipt; synthetic text can never constitute a scientific reference.
"""
from __future__ import annotations

import hashlib
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE = Path(__file__).with_name("source_admission.py")
SPEC = importlib.util.spec_from_file_location("source_admission", MODULE)
assert SPEC and SPEC.loader
admission = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = admission
SPEC.loader.exec_module(admission)


class _Page:
    def __init__(self, text: str):
        self.text = text

    def extract_text(self) -> str:
        return self.text


class _Reader:
    def __init__(self, _path: str, pages: list[str]):
        self.pages = [_Page(page) for page in pages]


class SourceAdmissionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "primary.pdf"
        self.path.write_bytes(b"checksum-bound fixture")
        self.original_sha = admission.PRIMARY_SOURCE_SHA256
        admission.PRIMARY_SOURCE_SHA256 = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def tearDown(self) -> None:
        admission.PRIMARY_SOURCE_SHA256 = self.original_sha
        self.temp.cleanup()

    def _pages(self, omit: str | None = None) -> list[str]:
        statements = admission.REQUIRED_SCOPE_STATEMENTS
        return [text for label, text in statements.items() if label != omit]

    def test_qualified_scope_is_explicitly_ineligible(self) -> None:
        result = admission.audit_primary_source(
            self.path, reader_factory=lambda path: _Reader(path, self._pages())
        )
        self.assertFalse(result.eligible)
        self.assertEqual(set(result.statement_pages), set(admission.REQUIRED_SCOPE_STATEMENTS))
        self.assertEqual(result.statement_pages["qualitative verification"], (1,))
        self.assertEqual(len(result.missing_observables), 5)
        self.assertEqual(
            result.wall_force_snapshots,
            (("experiment-A", 0.39), ("experiment-B", 0.33),
             ("BALL-1", 0.43), ("BALL-2", 0.40),
             ("BALL-3", 0.46), ("BALL-4", 0.47)),
        )

    def test_missing_affirmative_scope_statement_is_rejected(self) -> None:
        with self.assertRaisesRegex(admission.SourceAdmissionError, "control mismatch"):
            admission.audit_primary_source(
                self.path,
                reader_factory=lambda path: _Reader(path, self._pages("control mismatch")),
            )

    def test_checksum_mismatch_is_rejected_before_pdf_reading(self) -> None:
        admission.PRIMARY_SOURCE_SHA256 = "0" * 64
        with self.assertRaisesRegex(admission.SourceAdmissionError, "checksum mismatch"):
            admission.audit_primary_source(self.path, reader_factory=lambda path: (_ for _ in ()).throw(AssertionError()))


if __name__ == "__main__":
    unittest.main()
