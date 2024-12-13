"""Tests of the summary module."""

from pathlib import Path

import pandas as pd
import pytest

from isoslam import summary

BASE_DIR = Path.cwd()
RESOURCES = BASE_DIR / "tests" / "resources"

# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments


@pytest.mark.parametrize(
    ("pattern", "shape", "columns", "unique_filenames", "start_min", "start_max", "end_min", "end_max"),
    [
        pytest.param(
            "tests/**/*.tsv",
            (1696, 11),
            {
                "Read_UID",
                "Transcript_id",
                "Start",
                "End",
                "Chr",
                "Strand",
                "Assignment",
                "Conversions",
                "Convertible",
                "Coverage",
                "filename",
            },
            {"d0_0hr1_filtered_remapped_sorted", "d0_no4sU_filtered_remapped_sorted"},
            14940,
            25004,
            15080,
            29601,
            id="tsv files",
        ),
    ],
)
def test_append_files(
    pattern: str,
    shape: tuple[int, int],
    columns: list[str],
    unique_filenames: list[str],
    start_min: int,
    start_max: int,
    end_min: int,
    end_max: int,
) -> None:
    """Test appending of files."""
    all_files = summary.append_files(pattern)
    assert isinstance(all_files, pd.DataFrame)
    assert all_files.shape == shape
    assert set(all_files.columns) == columns
    assert set(all_files["filename"].unique()) == unique_filenames
    assert all_files["Start"].min() == start_min
    assert all_files["Start"].max() == start_max
    assert all_files["End"].min() == end_min
    assert all_files["End"].max() == end_max


@pytest.mark.parametrize(
    ("pattern", "expected", "count_max", "count_min"),
    [
        pytest.param(
            "tests/**/*.tsv",
            {"d0_0hr1_filtered_remapped_sorted", "d0_no4sU_filtered_remapped_sorted"},
            180,
            1,
            id="tsv files",
        )
    ],
)
def test_summary_counts(pattern: str, expected: set, count_max: int, count_min: int) -> None:
    """Test summary counts are correctly calculated."""
    summary_counts = summary.summary_counts(pattern)
    assert isinstance(summary_counts, pd.DataFrame)
    assert set(summary_counts["filename"].unique()) == expected
    assert summary_counts["count"].max() == count_max
    assert summary_counts["count"].min() == count_min
