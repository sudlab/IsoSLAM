"""Test the io.py module."""

import argparse
from datetime import datetime
from pathlib import Path

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


def test_load_bam() -> None:
    """Test loading of bam file."""
    # bam = io.load_bam()
    pass


def test_load_bed() -> None:
    """Test loading of bed file."""
    # bed = io.load_bam()
    pass


def test_load_gtf() -> None:
    """Test loading of gtf file."""
    # gtf = io.load_bed()
    pass


def test_load_vcf() -> None:
    """Test loading of vcf file."""
    # vcf = io.load_vcf()
    pass
