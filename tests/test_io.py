"""Test the io.py module."""

import argparse
from collections.abc import Callable
from datetime import datetime
from io import TextIOWrapper
from pathlib import Path
from typing import Any, TextIO

import pandas as pd
import polars as pl
import pysam
import pytest

from isoslam import io

# pylint: disable=protected-access

BASE_DIR = Path.cwd()
RESOURCES = BASE_DIR / "tests" / "resources"


CONFIG = {
    "this": "is",
    "a": "test",
    "yaml": "file",
    "int": 123,
    "float": 3.1415,
    "logical": True,
    "nested": {"something": "else"},
    "a_list": [1, 2, 3],
}


def test_get_date_time() -> None:
    """Test the fetching of a formatted date and time string."""
    assert datetime.strptime(io._get_date_time(), "%Y-%m-%d %H:%M:%S")


def test_str_to_path(tmp_path: Path) -> None:
    """Test string objects are converted to Path objects."""
    test_dir = str(tmp_path)
    converted_path = io._str_to_path(test_dir)

    assert isinstance(converted_path, Path)
    assert tmp_path == converted_path


def test_path_to_str(tmp_path: Path) -> None:
    """Test that Path objects in dictionaries are converted to strings."""
    CONFIG_PATH = {
        "this": "is",
        "a": "test",
        "with": tmp_path,
        "and": {"nested": tmp_path / "nested"},
    }
    CONFIG_STR = io._path_to_str(CONFIG_PATH)

    assert isinstance(CONFIG_STR, dict)
    assert isinstance(CONFIG_STR["with"], str)
    assert CONFIG_STR["with"] == str(tmp_path)
    assert isinstance(CONFIG_STR["and"]["nested"], str)
    assert CONFIG_STR["and"]["nested"] == str(tmp_path / "nested")


def test_read_yaml() -> None:
    """Test reading of YAML file."""
    sample_config = io.read_yaml(RESOURCES / "test.yaml")

    assert sample_config == CONFIG


# def test_create_config(tmp_path: Path) -> None:
#     """Test creation of configuration file from default."""


@pytest.mark.parametrize(
    ("filename", "config", "expected_filename"),
    [
        ("test_config_with_comments.yaml", None, "test_config_with_comments.yaml"),
        ("test_config_with_comments", None, "test_config_with_comments.yaml"),
        (None, "default", "config.yaml"),
        (None, None, "config.yaml"),
    ],
)
def test_write_config_with_comments(tmp_path: Path, filename: str, config: str, expected_filename: str) -> None:
    """Test writing of config file with comments.

    If and when specific configurations for different sample types are introduced then the parametrisation can be
    extended to allow these adding their names under "config" and introducing specific parameters that may differe
    between the configuration files.
    """
    # Setup argparse.Namespace with the tests parameters
    args = argparse.Namespace()
    args.filename = filename
    args.output_dir = tmp_path
    args.config = config
    args.simple = False

    # Write default config with comments to file
    io.create_config(args)

    # Read the written config
    with Path.open(tmp_path / expected_filename, encoding="utf-8") as f:
        written_config = f.read()

    # Validate that the written config has comments in it
    assert "Config file generated" in written_config
    assert "For more information on configuration and how to use it" in written_config
    # Validate some of the parameters are present
    assert "log_level:" in written_config
    assert "output_dir: output" in written_config
    assert "bam_file: data/bam/.bam" in written_config
    assert "vcf_file: data/vcf/" in written_config


def test_write_yaml(tmp_path: Path) -> None:
    """Test reading of YAML file."""
    io.write_yaml(
        config=CONFIG,
        output_dir=tmp_path,
        config_file="test.yaml",
        header_message="This is a test YAML configuration file",
    )
    outfile = tmp_path / "test.yaml"
    assert outfile.is_file()


@pytest.mark.parametrize(
    ("file_path", "object_type", "description", "compression"),
    [
        pytest.param(
            RESOURCES / "bam" / "d0_0hr1_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            "BAM version 1 compressed sequence data",
            "BGZF",
            id="file 0hr1 as Path",
        ),
        pytest.param(
            "tests/resources/bam/d0_0hr1_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            "BAM version 1 compressed sequence data",
            "BGZF",
            id="file 0hr1 as str",
        ),
        pytest.param(
            RESOURCES / "bam" / "d0_no4sU_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            "BAM version 1 compressed sequence data",
            "BGZF",
            id="file no4sU as Path",
        ),
        pytest.param(
            "tests/resources/bam/d0_no4sU_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            "BAM version 1 compressed sequence data",
            "BGZF",
            id="file no4sU as str",
        ),
    ],
)
def test_load_bam(
    file_path: str | Path,
    object_type: pysam.libcalignmentfile.AlignmentFile,
    description: str,
    compression: str,
) -> None:
    """Test loading of bam file."""
    bam_file = io._load_bam(file_path)
    assert isinstance(bam_file, object_type)
    assert bam_file.description == description
    assert bam_file.compression == compression


@pytest.mark.parametrize(
    ("file_path", "object_type"),
    [
        pytest.param(RESOURCES / "bed" / "test_coding_introns.bed", TextIOWrapper, id="bed file as Path"),
        pytest.param("tests/resources/bed/test_coding_introns.bed", TextIOWrapper, id="bed file as str"),
    ],
)
def test_load_bed(file_path: str | Path, object_type: TextIO) -> None:
    """Test loading of bed file."""
    bed_file = io._load_bed(file_path)
    assert isinstance(bed_file, object_type)


@pytest.mark.parametrize(
    ("file_path", "object_type"),
    [
        pytest.param(
            RESOURCES / "gtf" / "test_wash1.gtf", pysam.libctabix.tabix_generic_iterator, id="gtf file as Path"
        ),
        pytest.param(
            "tests/resources/gtf/test_wash1.gtf", pysam.libctabix.tabix_generic_iterator, id="gtf file as str"
        ),
    ],
)
def test_load_gtf(file_path: str | Path, object_type: str) -> None:
    """Test loading of gtf file."""
    gtf_file = io._load_gtf(file_path)
    assert isinstance(gtf_file, object_type)


@pytest.mark.parametrize(
    ("file_path", "object_type", "compression", "is_remote"),
    [
        pytest.param(RESOURCES / "vcf" / "d0.vcf.gz", pysam.libcbcf.VariantFile, "BGZF", False, id="d0 as Path"),
        pytest.param("tests/resources/vcf/d0.vcf.gz", pysam.libcbcf.VariantFile, "BGZF", False, id="d0 as str"),
    ],
)
def test_load_vcf(
    file_path: str | Path, object_type: pysam.libcbcf.VariantFile, compression: str, is_remote: bool
) -> None:
    """Test loading of vcf file."""
    vcf_file = io._load_vcf(file_path)
    assert isinstance(vcf_file, object_type)
    assert vcf_file.compression == compression
    assert vcf_file.is_remote == is_remote


@pytest.mark.parametrize(
    ("file_ext", "function"),
    [
        pytest.param(".bam", io._load_bam, id="bam extension"),
        pytest.param(".bed", io._load_bed, id="bed extension"),
        pytest.param(".gtf", io._load_gtf, id="gtf extension"),
        pytest.param(".vcf", io._load_vcf, id="vcf extension"),
        pytest.param(".vcf.gz", io._load_vcf, id="vcf.gz extension"),
    ],
)
def test_get_loader(file_ext: str, function: Callable) -> None:
    """Test io._get_loader returns the correct function."""
    assert io._get_loader(file_ext) == function


@pytest.mark.parametrize(
    ("file_ext"),
    [
        pytest.param(".csv", id="csv extension"),
        pytest.param(".txt", id="txt extension"),
        pytest.param(".tar.gz", id="tar.gz extension"),
    ],
)
def test_get_loader_value_error(file_ext: str) -> None:
    """Test io._get_loader() raises ValueError with incorrect file extension."""
    with pytest.raises(ValueError):  # noqa: PT011
        io._get_loader(file_ext)


@pytest.mark.parametrize(
    ("file_path", "object_type"),
    [
        pytest.param(
            RESOURCES / "bam" / "d0_0hr1_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            id="bam file 0hr1 as Path",
        ),
        pytest.param(
            "tests/resources/bam/d0_0hr1_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            id="bam file 0hr1 as str",
        ),
        pytest.param(
            RESOURCES / "bam" / "d0_no4sU_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            id="bam file no4sU as Path",
        ),
        pytest.param(
            "tests/resources/bam/d0_no4sU_filtered_remapped_sorted.bam",
            pysam.libcalignmentfile.AlignmentFile,
            id="bam file no4sU as str",
        ),
        # pytest.param(RESOURCES / "bed" / "test_coding_introns.bed", id="bed file as Path"),
        # pytest.param("tests/resources/bed/test_coding_introns.bed", id="bed file as str"),
        pytest.param(
            RESOURCES / "gtf" / "test_wash1.gtf", pysam.libctabix.tabix_generic_iterator, id="gtf file as Path"
        ),
        pytest.param(
            "tests/resources/gtf/test_wash1.gtf", pysam.libctabix.tabix_generic_iterator, id="gtf file as str"
        ),
        pytest.param(RESOURCES / "vcf" / "d0.vcf.gz", pysam.libcbcf.VariantFile, id="d0 as Path"),
        pytest.param("tests/resources/vcf/d0.vcf.gz", pysam.libcbcf.VariantFile, id="d0 as str"),
    ],
)
def test_load_file(file_path: str | Path, object_type: Any) -> None:
    """Test loading of files."""
    file = io.load_file(file_path)
    assert isinstance(file, object_type)


@pytest.mark.parametrize(
    ("pattern", "expected"),
    [
        pytest.param(
            "tests/resources/results/*.tsv",
            {
                Path(RESOURCES / "results" / "d0_0hr1_filtered_remapped_sorted.tsv"),
                Path(RESOURCES / "results" / "d0_no4sU_filtered_remapped_sorted.tsv"),
            },
            id="tsv files",
        ),
        pytest.param(
            "tests/resources/csv/output/*.csv",
            {
                Path(RESOURCES / "csv" / "output" / "d0_0hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_0hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_0hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_0hr4.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_12hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_12hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_12hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_12hr4.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_3hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_3hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_3hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_3hr4.csv"),
                Path(RESOURCES / "csv" / "output" / "d0_no4sU.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_0hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_0hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_0hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_12hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_12hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_3hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_3hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_3hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d16_no4sU.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_0hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_0hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_0hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_0hr4.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_12hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_12hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_12hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_12hr4.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_3hr1.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_3hr2.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_3hr3.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_3hr4.csv"),
                Path(RESOURCES / "csv" / "output" / "d2_no4sU.csv"),
            },
            id="csv files",
        ),
        pytest.param(
            "tests/**/*.gtf",
            {
                Path(RESOURCES / "gtf", "test_wash1.gtf"),
            },
            id="gtf files",
        ),
    ],
)
def test_find_files(pattern: str, expected: list[Path]) -> None:
    """Test finding of files."""
    files_found = set(io._find_files(pattern))
    assert files_found == expected


@pytest.mark.parametrize(
    ("pattern", "directory", "columns", "expected", "expected_shape"),
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
            (200, 11),
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
            (200, 11),
            id="parquet",
        ),
        pytest.param(
            ".csv",
            RESOURCES / "csv",
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
            (200, 11),
            id="csv",
        ),
        pytest.param(
            ".parquet",
            RESOURCES / "parquet",
            [
                "Read_UID",
                "Transcript_id",
                "Chr",
                "Strand",
                "Assignment",
                "Conversions",
                "Convertible",
                "Coverage",
            ],
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
            (200, 9),
            id="parquet subset of columns",
        ),
        pytest.param(
            ".csv",
            RESOURCES / "csv",
            [
                "Read_UID",
                "Transcript_id",
                "Conversions",
                "Convertible",
                "Coverage",
            ],
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
            (200, 6),
            id="csv subset of columns",
        ),
    ],
)
def test_load_output_files(
    pattern: str, directory: str | Path, columns: list[str], expected: set, expected_shape: tuple[int, int]
) -> None:
    """Test loading of files."""
    files_found = io.load_output_files(pattern, directory, columns)
    assert set(files_found.keys()) == expected
    for _, dataframe in files_found.items():
        assert isinstance(dataframe, pl.DataFrame)
        assert dataframe.shape == expected_shape


@pytest.mark.parametrize(
    ("data", "outfile", "separator"),
    [
        pytest.param(pd.DataFrame({"test1": [0, 1], "test2": [2, 3]}), "results.parquet", ",", id="Pandas to parquet"),
        pytest.param(pd.DataFrame({"test1": [0, 1], "test2": [2, 3]}), "results.tsv", "\t", id="Pandas to tsv"),
        pytest.param(pd.DataFrame({"test1": [0, 1], "test2": [2, 3]}), "results.csv", ",", id="Pandas to csv"),
        pytest.param(pl.DataFrame({"test1": [0, 1], "test2": [2, 3]}), "results.parquet", ",", id="Polars to parquet"),
        pytest.param(pl.DataFrame({"test1": [0, 1], "test2": [2, 3]}), "results.tsv", "\t", id="Polars to tsv"),
        pytest.param(pl.DataFrame({"test1": [0, 1], "test2": [2, 3]}), "results.csv", ",", id="Polars to csv"),
    ],
)
def test_data_frame_to_file(data: pd.DataFrame | pl.DataFrame, outfile: str, separator: str, tmp_path: Path) -> None:
    """Test data_frame_to_file() for different formats."""
    io.data_frame_to_file(data, output_dir=tmp_path, outfile=outfile, sep=separator)
    outdir_file = tmp_path / outfile
    assert outdir_file.is_file()


@pytest.mark.parametrize(
    ("not_dataframe"),
    [
        pytest.param("string", id="string"),
        pytest.param(10, id="int"),
        pytest.param([1, 2, 3], id="list"),
        pytest.param((1, 2, 3), id="tuple"),
    ],
)
def test_data_frame_to_file_typeerror(not_dataframe: Any, tmp_path: Path) -> None:
    """Test TypeError is raised by data_frame_to_file() when data is not Polars or Pandas DataFrame."""
    with pytest.raises(TypeError):
        io.data_frame_to_file(data=not_dataframe, output_dir=tmp_path, outfile="test.csv", sep=",")


@pytest.mark.parametrize(
    ("assigned_conversions", "coverage_counts", "read_uid", "assignment", "delim", "expected"),
    [
        pytest.param(
            {("read1", (1, 2, 9, "-")), ("read2", (40, 60, 9, "+"))},
            {"converted_position": 4, "convertible": 8, "coverage": 16},
            10,
            "Ret",
            "\t",
            pd.DataFrame.from_dict(
                {
                    0: [10, 10],
                    1: ["read2", "read1"],
                    2: [40, 1],
                    3: [60, 2],
                    4: [9, 9],
                    5: ["+", "-"],
                    6: ["Ret", "Ret"],
                    7: [4, 4],
                    8: [8, 8],
                    9: [16, 16],
                }
            ),
            id="dummy constructs",
        ),
    ],
)
def test_write_assigned_conversions(  # pylint: disable=too-many-positional-arguments
    assigned_conversions: dict[str, set[Any]],
    coverage_counts: dict[str, int],
    read_uid: int | str,
    assignment: str,
    delim: str,
    expected: pd.DataFrame,
    tmp_path: Path,
) -> None:
    """Test files are written correctly."""
    with Path.open(tmp_path / "test_output.tsv", "w", encoding="utf-8") as outfile:
        io.write_assigned_conversions(
            assigned_conversions,
            coverage_counts,
            read_uid,
            assignment,
            outfile,
            delim,
        )

        output = tmp_path / "test_output.tsv"
        assert output.is_file
    result = pd.read_csv(output, sep=delim, header=None)
    # Ensure order is consistent with expected order
    result.sort_values(1, ascending=False, inplace=True, ignore_index=True)
    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    ("raw_dict", "expected"),
    [
        pytest.param({"string": "str"}, {"string": str}, id="string"),
        pytest.param({"int": "int"}, {"int": int}, id="int"),
        pytest.param({"float": "float"}, {"float": float}, id="float"),
        pytest.param({"dict": "dict"}, {"dict": dict}, id="dict"),
        pytest.param({"list": "list"}, {"list": list}, id="list"),
        pytest.param({"tuple": "tuple"}, {"tuple": tuple}, id="tuple"),
        pytest.param(
            {"string": "str", "int": "int", "float": "float", "dict": "dict", "list": "list", "tuple": "tuple"},
            {"string": str, "int": int, "float": float, "dict": dict, "list": list, "tuple": tuple},
            id="multiple",
        ),
    ],
)
def test_type_schema(raw_dict: dict[str, str], expected: dict[str, type]) -> None:
    """Test conversion of schema with strings to types."""
    converted = io._type_schema(raw_dict)
    assert converted == expected
