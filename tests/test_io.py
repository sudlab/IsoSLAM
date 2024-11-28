"""Test the io.py module."""

import argparse
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import Any

import pysam
import pytest

from isoslam import io

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
    """Test that Path objects are converted to strings."""
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


def test_create_config(tmp_path: Path) -> None:
    """Test creation of configuration file from default."""


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


@pytest.mark.skip
@pytest.mark.parametrize(
    ("file_path", "object_type"),
    [
        pytest.param(RESOURCES / "bed" / "test_coding_introns.bed", str, id="bed file as Path"),
        pytest.param("tests/resources/bed/test_coding_introns.bed", str, id="bed file as str"),
    ],
)
def test_load_bed(file_path: str | Path, object_type: str) -> None:
    """Test loading of bed file."""
    bed_file = io._load_bam()
    assert isinstance(bed_file, object_type)


@pytest.mark.skip
@pytest.mark.parametrize(
    ("file_path", "object_type"),
    [
        pytest.param(RESOURCES / "gtf" / "test_wash1.gtf", str, id="gtf file as Path"),
        pytest.param("tests/resources/gtf/test_wash1.gtf", str, id="gtf file as str"),
    ],
)
def test_load_gtf(file_path: str | Path, object_type: str) -> None:
    """Test loading of gtf file."""
    gtf_file = io._load_bed()
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
        # pytest.param(RESOURCES / "gtf" / "test_wash1.gtf", id="gtf file as Path"),
        # pytest.param("tests/resources/gtf/test_wash1.gtf", id="gtf file as str"),
        pytest.param(RESOURCES / "vcf" / "d0.vcf.gz", pysam.libcbcf.VariantFile, id="d0 as Path"),
        pytest.param("tests/resources/vcf/d0.vcf.gz", pysam.libcbcf.VariantFile, id="d0 as str"),
    ],
)
def test_load_file(file_path: str | Path, object_type: Any) -> None:
    """Test loading of files."""
    file = io.load_file(file_path)
    assert isinstance(file, object_type)
