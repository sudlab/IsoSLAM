"""Tests for the isoslam module."""

from pathlib import Path
from typing import Any

import pytest  # type: ignore[import-not-found]

from isoslam import isoslam

BASE_DIR = Path.cwd()
RESOURCES = BASE_DIR / "tests" / "resources"


@pytest.mark.parametrize(
    ("bed_file", "expected_transcript"),
    [
        pytest.param(  # type: ignore[misc]
            RESOURCES / "bed" / "test_coding_introns.bed",
            {
                "ENST00000442898": [
                    ("9", 14940, 15080, "ENST00000442898", "-"),
                    ("9", 15149, 15908, "ENST00000442898", "-"),
                    ("9", 16061, 16717, "ENST00000442898", "-"),
                    ("9", 16876, 16964, "ENST00000442898", "-"),
                    ("9", 17166, 17343, "ENST00000442898", "-"),
                    ("9", 17479, 17718, "ENST00000442898", "-"),
                    ("9", 17855, 18027, "ENST00000442898", "-"),
                    ("9", 18174, 18380, "ENST00000442898", "-"),
                    ("9", 18492, 24850, "ENST00000442898", "-"),
                    ("9", 25004, 29601, "ENST00000442898", "-"),
                ]
            },
            id="bed coding introons",
        ),
    ],
)
def test_isoslam_extract_transcripts(
    bed_file: str | Path,
    expected_transcript: dict[Any, list[tuple[Any, int, int, Any, Any]]],
) -> None:
    """Test extraction of tanscript data from bed file using extract_transcripts()."""
    assert isoslam.extract_transcripts(bed_file) == expected_transcript
