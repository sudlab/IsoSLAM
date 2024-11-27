"""Tests of the processing modules entry point, argument parsing and ability to correctly select programs."""

from collections.abc import Callable
from pathlib import Path

import pytest

from isoslam import io, processing


@pytest.mark.parametrize(
    ("option", "help_message"),
    [
        pytest.param(["-h"], "usage:", id="Single letter help"),
        pytest.param(["--help"], "program", id="Multiple letter help"),
        pytest.param(
            ["process"],
            "Process all files and run all summary plotting and statistics.",
            id="process module",
            marks=pytest.mark.xfail(reason="work in progress"),
        ),
    ],
)
def test_entry_point_help(option: str, help_message: str, capsys: pytest.CaptureFixture) -> None:
    """Test the entry_point()'s help argument."""
    try:
        processing.entry_point(manually_provided_args=option)
    except SystemExit:
        pass
    output = capsys.readouterr().out
    assert help_message in output


@pytest.mark.parametrize(
    ("options", "expected_function", "expected_args"),
    [
        pytest.param(
            ["--config", "dummy/config/dir/config.yaml", "process"],
            processing.process,
            {"config_file": Path("dummy/config/dir/config.yaml")},
            id="dummy config and process",
        ),
        pytest.param(
            ["--config", "dummy/config/dir/config.yaml", "process", "--bam-file", "data/bam/some.bam"],
            processing.process,
            {"config_file": Path("dummy/config/dir/config.yaml"), "bam_file": Path("data/bam/some.bam")},
            id="dummy config, bam file and process",
        ),
        pytest.param(
            [
                "--config",
                "dummy/config/dir/config.yaml",
                "--output-dir",
                "output",
                "process",
                "--bam-file",
                "data/bam/some.bam",
            ],
            processing.process,
            {
                "config_file": Path("dummy/config/dir/config.yaml"),
                "output_dir": Path("output"),
                "bam_file": Path("data/bam/some.bam"),
            },
            id="dummy config and output, process with bam file",
        ),
        pytest.param(
            [
                "--config",
                "dummy/config/dir/config.yaml",
                "--output-dir",
                "output",
                "process",
                "--bam-file",
                "data/bam/some.bam",
                "--gtf-file",
                "data/gtf/some.gtf",
                "--bed-file",
                "data/bed/some.bed",
                "--vcf-file",
                "data/vcf/some.vcf.gz",
            ],
            processing.process,
            {
                "config_file": Path("dummy/config/dir/config.yaml"),
                "output_dir": Path("output"),
                "bam_file": Path("data/bam/some.bam"),
                "gtf_file": Path("data/gtf/some.gtf"),
                "bed_file": Path("data/bed/some.bed"),
                "vcf_file": Path("data/vcf/some.vcf.gz"),
            },
            id="dummy config and output, process with all files",
        ),
        pytest.param(
            [
                "create-config",
                "--filename",
                "test_config.yaml",
                "--output-dir",
                "custom_configs",
            ],
            io.create_config,
            {
                "filename": Path("test_config.yaml"),
                "output_dir": Path("custom_configs"),
            },
            id="create config file",
        ),
    ],
)
def test_entry_point_sub_parsers(options: str, expected_function: Callable, expected_args: dict) -> None:
    """Test the correct function is returned by each subparser."""
    returned_args = processing.entry_point(manually_provided_args=options, testing=True)
    assert returned_args.func == expected_function
    returned_args_dict = vars(returned_args)
    for argument, value in expected_args.items():
        assert returned_args_dict[argument] == value
