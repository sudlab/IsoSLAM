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
# pylint: disable=protected-access


@pytest.mark.parametrize(
    ("file_ext", "directory", "shape", "unique_filenames", "start_min", "start_max", "end_min", "end_max"),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "results",
            (1696, 11),
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
            (7000, 11),
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
    unique_filenames: list[str],
    start_min: int,
    start_max: int,
    end_min: int,
    end_max: int,
) -> None:
    """Test appending of files."""
    all_files = summary.append_files(file_ext, directory)
    assert isinstance(all_files, pl.DataFrame)
    assert all_files.shape == shape
    assert set(all_files["filename"].unique()) == unique_filenames
    assert all_files["Start"].min() == start_min
    assert all_files["Start"].max() == start_max
    assert all_files["End"].min() == end_min
    assert all_files["End"].max() == end_max


@pytest.mark.parametrize(
    ("file_ext", "directory", "regex", "expected_files", "expected_numbers"),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "results",
            r"^d(\w+)_(\w+)hr(\w+)_",
            {"d0_0hr1_filtered_remapped_sorted"},
            {
                "count_max": 180,
                "count_min": 1,
                "total_max": 267,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.07142857142857142,
            },
            id="tsv",
        ),
        pytest.param(
            ".parquet",
            RESOURCES / "parquet",
            r"^d(\w+)_(\w+)hr(\w+)",
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
                "d16_0hr1",
                "d16_0hr2",
                "d16_0hr3",
                "d16_12hr1",
                "d16_12hr3",
                "d16_3hr1",
                "d16_3hr2",
                "d16_3hr3",
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
def test_summary_counts(
    file_ext: str, directory: Path, regex: Pattern, expected_files: set, expected_numbers, regtest
) -> None:
    """Test summary counts are correctly calculated."""
    summary_counts = summary.summary_counts(file_ext, directory, regex=regex)
    assert isinstance(summary_counts, pl.DataFrame)
    assert set(summary_counts["filename"].unique()) == expected_files
    assert summary_counts["conversion_count"].max() == expected_numbers["count_max"]
    assert summary_counts["conversion_count"].min() == expected_numbers["count_min"]
    assert summary_counts["conversion_total"].max() == expected_numbers["total_max"]
    assert summary_counts["conversion_total"].min() == expected_numbers["total_min"]
    assert summary_counts["conversion_percent"].max() == expected_numbers["percent_max"]
    assert summary_counts["conversion_percent"].min() == expected_numbers["percent_min"]
    print(summary_counts.write_csv(), file=regtest)


@pytest.mark.parametrize(
    ("df", "column", "regex"),
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
            id="realistic",
        ),
    ],
)
def test_extract_hour_and_replicate(df: pl.DataFrame, column: str, regex: Pattern, regtest) -> None:
    """Test extraction of hour and replicate from filename."""
    summary_counts = summary.extract_day_hour_and_replicate(df, column, regex)
    print(summary_counts.write_csv(), file=regtest)


@pytest.mark.parametrize(
    ("df", "groupby", "converted"),
    [
        pytest.param(
            pl.DataFrame(
                {
                    "transcript_id": ["a", "a", "b", "b", "c", "d", "d", "e"],
                    "chr": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr4", "chr4", "chr5"],
                    "conversion": [True, False, True, False, False, False, True, True],
                    "n": [1, 2, 3, 4, 5, 6, 7, 8],
                }
            ),
            ["transcript_id", "chr"],
            "conversion",
            id="missing conversion True for transcript_id on chr 3",
        ),
        pytest.param(
            "sample_data_summary_counts",
            [
                "Transcript_id",
                "Strand",
                "Start",
                "End",
                "Assignment",
                "day",
                "hour",
                "replicate",
            ],
            "one_or_more_conversion",
            id="real data",
        ),
    ],
)
def test_aggregate_conversions(
    df: pl.DataFrame | str,
    groupby: list[str],
    converted: str,
    request: pytest.FixtureRequest,
    regtest,
) -> None:
    """Test derivation on non-captured dataset."""
    df = request.getfixturevalue(df) if isinstance(df, str) else df
    aggregated_conversions = summary._aggregate_conversions(df, groupby, converted)
    print(aggregated_conversions.write_csv(), file=regtest)


@pytest.mark.parametrize(
    ("df", "groupby", "converted"),
    [
        pytest.param(
            pl.DataFrame(
                {
                    "transcript_id": ["a", "b", "c", "d", "e"],
                    "chr": ["chr1", "chr2", "chr3", "chr4", "chr5"],
                    "len": [2, 2, 1, 2, 1],
                    "conversion": [True, True, False, True, True],
                },
                schema={
                    "transcript_id": pl.datatypes.String,
                    "chr": pl.datatypes.String,
                    "len": pl.datatypes.UInt32,
                    "conversion": pl.datatypes.Boolean,
                },
            ),
            ["transcript_id", "chr"],
            "conversion",
            id="missing conversion True for transcript_id on chr 3",
        ),
        pytest.param(
            "sample_data_summary_counts",
            None,
            "one_or_more_conversion",
            id="real data",
        ),
    ],
)
def test_filter_no_conversions(
    df: pl.DataFrame | str,
    groupby: list[str],
    converted: str,
    request: pytest.FixtureRequest,
    regtest,
) -> None:
    """Test filtering of non conversions."""
    df = request.getfixturevalue(df) if isinstance(df, str) else df
    filtered_no_conversions = summary._filter_no_conversions(df, groupby, converted)
    print(filtered_no_conversions.write_csv(), file=regtest)


@pytest.mark.parametrize(
    ("df", "groupby", "converted"),
    [
        pytest.param(
            "sample_data_summary_counts",
            None,
            "one_or_more_conversion",
            id="real data groupby none",
        ),
    ],
)
def test_get_one_or_more_conversion(
    df: pl.DataFrame, groupby: list[str], converted: str, request: pytest.FixtureRequest, regtest
) -> None:
    """Test that inner_join_no_conversions returns the correct subset."""
    df = request.getfixturevalue(df) if isinstance(df, str) else df
    one_or_more_conversions = summary._get_one_or_more_conversion(df, groupby, converted)
    print(one_or_more_conversions.write_csv(), file=regtest)


@pytest.mark.parametrize(
    ("df", "groupby", "average"),
    [
        pytest.param(
            pl.DataFrame(
                {
                    "day": [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
                    "hour": [0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 16, 16, 16, 16],
                    "replicate": [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4],
                    "conversion_total": [8, 8, 8, 8, 1, 2, 3, 4, 10, 20, 30, 40, 0, 0, 0, 1],
                    "conversion_percent": [
                        10.0,
                        10.0,
                        10.0,
                        10.0,
                        5.0,
                        6.0,
                        7.0,
                        8.0,
                        10.0,
                        20.0,
                        30.0,
                        40.0,
                        0.5,
                        0.5,
                        0.8,
                        0.9,
                    ],
                }
            ),
            ["day", "hour"],
            "mean",
            id="simple mean",
        ),
        pytest.param(
            pl.DataFrame(
                {
                    "day": [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
                    "hour": [0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 16, 16, 16, 16],
                    "replicate": [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4],
                    "conversion_total": [8, 8, 8, 8, 1, 2, 3, 4, 10, 20, 30, 40, 0, 0, 0, 1],
                    "conversion_percent": [
                        10.0,
                        10.0,
                        10.0,
                        10.0,
                        5.0,
                        6.0,
                        7.0,
                        8.0,
                        10.0,
                        20.0,
                        30.0,
                        40.0,
                        0.5,
                        0.5,
                        0.8,
                        0.9,
                    ],
                }
            ),
            ["day", "hour"],
            "median",
            id="simple median",
        ),
        pytest.param(
            "test_average_replicates",
            None,
            "mean",
            id="real data mean",
        ),
        pytest.param(
            "test_average_replicates",
            None,
            "median",
            id="real data median",
        ),
    ],
)
def test_average_replicates(
    df: pl.DataFrame | str, groupby: list[str], average: str, request: pytest.FixtureRequest, regtest
) -> None:
    """Test the _average_replicate() function."""
    df = request.getfixturevalue(df) if isinstance(df, str) else df
    df_average = summary._average_replicates(df, groupby, average)
    print(df_average.write_csv(), file=regtest)


def test_average_replicates_valueerror() -> None:
    """Test raising of ValueError if invalid ''average'' parameter is passed."""
    with pytest.raises(ValueError, match="Invalid value for average"):
        summary._average_replicates(pl.DataFrame({"a": [1, 2]}), groupby=["a"], average="mode")


@pytest.mark.parametrize(
    ("df", "base_day", "base_hour"),
    [
        pytest.param(
            pl.DataFrame(
                {
                    "day": [0, 0, 1, 1],
                    "hour": [0, 8, 0, 16],
                    "conversion_total": [8.0, 2.5, 25.0, 0.25],
                    "conversion_percent": [10.0, 6.5, 25.0, 0.675],
                }
            ),
            0,
            0,
            id="simple",
        ),
        pytest.param("test_baseline_mean", 0, 0, id="real"),
    ],
)
def test_select_base_levels(
    df: pl.DataFrame | str,
    base_day: int,
    base_hour: int,
    request: pytest.FixtureRequest,
    regtest,
) -> None:
    """Test subsetting of data for baseline values."""
    df = request.getfixturevalue(df) if isinstance(df, str) else df
    baseline = summary._select_base_levels(df, base_day, base_hour)
    print(baseline.write_csv(), file=regtest)


@pytest.mark.parametrize(
    (
        "file_ext",
        "directory",
        "groupby",
        "conversions_var",
        "conversions_threshold",
        "test_file",
        "regex",
        "expected_numbers",
    ),
    [
        pytest.param(
            ".tsv",
            RESOURCES / "tsv",
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
            r"^d(\w+)_(\w+)hr(\w+)",
            {
                "shape": (6108, 14),
                "count_max": 5,
                "count_min": 1,
                "total_max": 6,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.16666666666666666,
                "unique_files": 32,
            },
            id="tsv",
        ),
        pytest.param(
            ".parquet",
            RESOURCES / "parquet",
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
            r"^d(\w+)_(\w+)hr(\w+)",
            {
                "shape": (6108, 14),
                "count_max": 5,
                "count_min": 1,
                "total_max": 6,
                "total_min": 1,
                "percent_max": 1.0,
                "percent_min": 0.16666666666666666,
                "unique_files": 32,
            },
            id="parquet",
        ),
    ],
)
def test_statistics_class(
    file_ext: str,
    directory: str,
    groupby: list[str],
    conversions_var: str,
    conversions_threshold: int,
    test_file: str,
    regex: Pattern,
    expected_numbers: dict[str, Any],
) -> None:
    """Test instantiation of statistics class."""
    statistics = summary.Statistics(
        file_ext=file_ext,
        directory=directory,
        groupby=groupby,
        conversions_var=conversions_var,
        conversions_threshold=conversions_threshold,
        test_file=test_file,
        regex=regex,
    )
    assert statistics.file_ext == file_ext
    assert statistics.directory == directory
    assert statistics.regex == regex
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
