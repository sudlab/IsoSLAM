"""Tests of the summary module."""

from pathlib import Path
from re import Pattern
from typing import Any

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
    ("file_ext", "directory", "expected_files", "expected_numbers"),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "results",
            {"d0_0hr1_filtered_remapped_sorted", "d0_no4sU_filtered_remapped_sorted"},
            {
                "count_max": 180,
                "count_min": 1,
                "total_max": 267,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.026490066225165563,
            },
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
            {
                "count_max": 5,
                "count_min": 1,
                "total_max": 6,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.16666666666666666,
            },
            id="parquet",
        ),
    ],
)
def test_summary_counts(file_ext: str, directory: Path, expected_files: set, expected_numbers) -> None:
    """Test summary counts are correctly calculated."""
    summary_counts = summary.summary_counts(file_ext, directory)
    assert isinstance(summary_counts, pl.DataFrame)
    assert set(summary_counts["filename"].unique()) == expected_files
    assert summary_counts["conversion_count"].max() == expected_numbers["count_max"]
    assert summary_counts["conversion_count"].min() == expected_numbers["count_min"]
    assert summary_counts["conversion_total"].max() == expected_numbers["total_max"]
    assert summary_counts["conversion_total"].min() == expected_numbers["total_min"]
    assert summary_counts["conversion_percent"].max() == expected_numbers["percent_max"]
    assert summary_counts["conversion_percent"].min() == expected_numbers["percent_min"]


@pytest.mark.parametrize(
    ("df", "column", "regex", "expected"),
    [
        pytest.param(
            pl.DataFrame(
                {
                    "filename": [
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
                    ]
                }
            ),
            "filename",
            r"^d(\w+)_(\w+)hr(\w+)",
            pl.DataFrame(
                {
                    "filename": [
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
                    ],
                    "day": [
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        None,
                        "16",
                        "16",
                        "16",
                        "16",
                        "16",
                    ],
                    "hour": [
                        "0",
                        "0",
                        "0",
                        "0",
                        "12",
                        "12",
                        "12",
                        "12",
                        "3",
                        "3",
                        "3",
                        "3",
                        None,
                        "0",
                        "0",
                        "0",
                        "12",
                        "12",
                    ],
                    "replicate": [
                        "1",
                        "2",
                        "3",
                        "4",
                        "1",
                        "2",
                        "3",
                        "4",
                        "1",
                        "2",
                        "3",
                        "4",
                        None,
                        "1",
                        "2",
                        "3",
                        "1",
                        "3",
                    ],
                }
            ),
            id="simple",
        ),
        pytest.param(
            pl.DataFrame(
                {
                    "filename": [
                        "d0_0hr1_EKRN230046546-1A_HFWGNDSX7_L2",
                        "d0_0hr2_EKRN230046547-1A_HFWGNDSX7_L2",
                        "d0_0hr3_EKRN230046548-1A_HFWGNDSX7_L2",
                        "d0_0hr4_EKRN230046549-1A_HFWGNDSX7_L2",
                        "d0_12hr1_EKRN230046554-1A_HFWGNDSX7_L2",
                        "d0_12hr2_EKRN230046555-1A_HFWGNDSX7_L2",
                        "d0_12hr3_EKRN230046556-1A_HFWGNDSX7_L2",
                        "d0_12hr4_EKRN230046557-1A_HFWGNDSX7_L2",
                        "d0_no4sU_EKRN230046554-1A_HFWGNDSX7_L2",
                    ]
                }
            ),
            "filename",
            r"^d(\w+)_(\w+)hr(\w+)_",
            pl.DataFrame(
                {
                    "filename": [
                        "d0_0hr1_EKRN230046546-1A_HFWGNDSX7_L2",
                        "d0_0hr2_EKRN230046547-1A_HFWGNDSX7_L2",
                        "d0_0hr3_EKRN230046548-1A_HFWGNDSX7_L2",
                        "d0_0hr4_EKRN230046549-1A_HFWGNDSX7_L2",
                        "d0_12hr1_EKRN230046554-1A_HFWGNDSX7_L2",
                        "d0_12hr2_EKRN230046555-1A_HFWGNDSX7_L2",
                        "d0_12hr3_EKRN230046556-1A_HFWGNDSX7_L2",
                        "d0_12hr4_EKRN230046557-1A_HFWGNDSX7_L2",
                        "d0_no4sU_EKRN230046554-1A_HFWGNDSX7_L2",
                    ],
                    "day": [
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        None,
                    ],
                    "hour": [
                        "0",
                        "0",
                        "0",
                        "0",
                        "12",
                        "12",
                        "12",
                        "12",
                        None,
                    ],
                    "replicate": [
                        "1",
                        "2",
                        "3",
                        "4",
                        "1",
                        "2",
                        "3",
                        "4",
                        None,
                    ],
                }
            ),
            id="realistic",
        ),
    ],
)
def test_extract_hour_and_replicate(df: pl.DataFrame, column: str, regex: Pattern, expected: pl.DataFrame) -> None:
    """Test extraction of hour and replicate from filename."""
    pl.testing.assert_frame_equal(summary.extract_day_hour_and_replicate(df, column, regex), expected)


@pytest.mark.parametrize(
    (
        "file_ext",
        "directory",
        "columns",
        "groupby",
        "conversions_var",
        "conversions_threshold",
        "test_file",
        "expected_numbers",
    ),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "tsv",
            [
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
            ],
            [
                "Transcript_id",
                "Start",
                "End",
                "Chr",
                "Strand",
                "Assignment",
                "filename",
            ],
            "Conversions",
            1,
            "no4sU",
            {
                "shape": (6647, 14),
                "count_max": 5,
                "count_min": 1,
                "total_max": 6,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.16666666666666666,
                "unique_files": 35,
            },
            id="tsv",
        ),
        pytest.param(
            ".parquet",
            RESOURCES / "parquet",
            [
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
            ],
            [
                "Transcript_id",
                "Start",
                "End",
                "Chr",
                "Strand",
                "Assignment",
                "filename",
            ],
            "Conversions",
            1,
            "no4sU",
            {
                "shape": (6647, 14),
                "count_max": 5,
                "count_min": 1,
                "total_max": 6,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.16666666666666666,
                "unique_files": 35,
            },
            id="parquet",
        ),
    ],
)
def test_statistics_class(
    file_ext: str,
    directory: str,
    columns: list[str],
    groupby: list[str],
    conversions_var: str,
    conversions_threshold: int,
    test_file: str,
    expected_numbers: dict[str, Any],
) -> None:
    """Test instantiation of statistics class."""
    statistics = summary.Statistics(
        file_ext, directory, columns, groupby, conversions_var, conversions_threshold, test_file
    )
    assert statistics.file_ext == file_ext
    assert statistics.directory == directory
    assert statistics.columns == columns
    assert statistics.groupby == groupby
    assert statistics.conversions_var == conversions_var
    assert statistics.conversions_threshold == conversions_threshold
    assert statistics.test_file == test_file
    # Check the data frame
    assert isinstance(statistics.data, pl.DataFrame)
    assert statistics.shape == expected_numbers["shape"]
    assert statistics.data["conversion_count"].max() == expected_numbers["count_max"]
    assert statistics.data["conversion_count"].min() == expected_numbers["count_min"]
    assert statistics.data["conversion_total"].max() == expected_numbers["total_max"]
    assert statistics.data["conversion_total"].min() == expected_numbers["total_min"]
    assert statistics.data["conversion_percent"].max() == expected_numbers["percent_max"]
    assert statistics.data["conversion_percent"].min() == expected_numbers["percent_min"]
    assert statistics.unique == expected_numbers["unique_files"]
