"""Tests of the summary module."""

from pathlib import Path

import polars as pl
import pytest

from isoslam import summary

BASE_DIR = Path.cwd()
RESOURCES = BASE_DIR / "tests" / "resources"

# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments


@pytest.mark.parametrize(
    ("file_ext", "directory", "shape", "columns", "unique_filenames", "start_min", "start_max", "end_min", "end_max"),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "results",
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
            },
            {"d0_0hr1_filtered_remapped_sorted", "d0_no4sU_filtered_remapped_sorted"},
            14940,
            25004,
            15080,
            29601,
            id="tsv",
        ),
        pytest.param(
            ".parquet",
            RESOURCES / "parquet",
            (7000, 10),
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
            },
            {
                "d0_0hr1",
                "d0_0hr2",
                "d0_0hr3",
                "d0_0hr4",
                "d0_12hr1",
                "d0_12hr2",
                "d0_12hr3",
                "d0_12hr4",
                "d0_3hr1",
                "d0_3hr2",
                "d0_3hr3",
                "d0_3hr4",
                "d0_no4sU",
                "d16_0hr1",
                "d16_0hr2",
                "d16_0hr3",
                "d16_12hr1",
                "d16_12hr3",
                "d16_3hr1",
                "d16_3hr2",
                "d16_3hr3",
                "d16_no4sU",
                "d2_0hr1",
                "d2_0hr2",
                "d2_0hr3",
                "d2_0hr4",
                "d2_12hr1",
                "d2_12hr2",
                "d2_12hr3",
                "d2_12hr4",
                "d2_3hr1",
                "d2_3hr2",
                "d2_3hr3",
                "d2_3hr4",
                "d2_no4sU",
            },
            189339,
            244409112,
            189967,
            244409243,
            id="parquet",
        ),
    ],
)
def test_append_files(
    file_ext: str,
    directory: Path,
    shape: tuple[int, int],
    columns: list[str],
    unique_filenames: list[str],
    start_min: int,
    start_max: int,
    end_min: int,
    end_max: int,
) -> None:
    """Test appending of files."""
    all_files = summary.append_files(file_ext, directory, list(columns))
    assert isinstance(all_files, pl.DataFrame)
    assert all_files.shape == shape
    columns.add("filename")
    assert set(all_files.columns) == columns
    assert set(all_files["filename"].unique()) == unique_filenames
    assert all_files["Start"].min() == start_min
    assert all_files["Start"].max() == start_max
    assert all_files["End"].min() == end_min
    assert all_files["End"].max() == end_max


@pytest.mark.parametrize(
    ("file_ext", "directory", "expected", "count_max", "count_min"),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "results",
            {"d0_0hr1_filtered_remapped_sorted", "d0_no4sU_filtered_remapped_sorted"},
            180,
            1,
            id="tsv",
        ),
        pytest.param(
            ".parquet",
            RESOURCES / "parquet",
            {
                "d0_0hr1",
                "d0_0hr2",
                "d0_0hr3",
                "d0_0hr4",
                "d0_12hr1",
                "d0_12hr2",
                "d0_12hr3",
                "d0_12hr4",
                "d0_3hr1",
                "d0_3hr2",
                "d0_3hr3",
                "d0_3hr4",
                "d0_no4sU",
                "d16_0hr1",
                "d16_0hr2",
                "d16_0hr3",
                "d16_12hr1",
                "d16_12hr3",
                "d16_3hr1",
                "d16_3hr2",
                "d16_3hr3",
                "d16_no4sU",
                "d2_0hr1",
                "d2_0hr2",
                "d2_0hr3",
                "d2_0hr4",
                "d2_12hr1",
                "d2_12hr2",
                "d2_12hr3",
                "d2_12hr4",
                "d2_3hr1",
                "d2_3hr2",
                "d2_3hr3",
                "d2_3hr4",
                "d2_no4sU",
            },
            5,
            1,
            id="parquet",
        ),
    ],
)
def test_summary_counts(file_ext: str, directory: Path, expected: set, count_max: int, count_min: int) -> None:
    """Test summary counts are correctly calculated."""
    summary_counts = summary.summary_counts(file_ext, directory)
    assert isinstance(summary_counts, pl.DataFrame)
    assert set(summary_counts["filename"].unique()) == expected
    assert summary_counts["count"].max() == count_max
    assert summary_counts["count"].min() == count_min
