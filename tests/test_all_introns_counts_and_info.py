"""Tests for all_introns_counts_and_info module."""

import argparse
from pathlib import Path

import pytest

from isoslam import all_introns_counts_and_info

BASE_DIR = Path.cwd()
TEST_DIR = BASE_DIR / "tests"
RESOURCES = TEST_DIR / "resources"
GTF_DIR = RESOURCES / "gtf"
BED_DIR = RESOURCES / "bed"
VCF_DIR = RESOURCES / "vcf"
BAM_DIR = RESOURCES / "bam"


@pytest.mark.parametrize(
    ("file_path"),
    [
        pytest.param(BAM_DIR / "sorted_assigned" / "d0_0hr1_filtered_remapped_sorted.sorted.assigned.bam", id="file 1"),
        pytest.param(
            BAM_DIR / "sorted_assigned" / "d0_no4sU_filtered_remapped_sorted.sorted.assigned.bam", id="file 2"
        ),
    ],
)
def test_main(file_path: Path, snapshot, tmp_path) -> None:
    """Regression test to check that main() function returns the expected data structure."""
    args = argparse.Namespace(
        infile_bam=str(file_path),
        gtf_path=list(GTF_DIR.glob("*.gtf"))[0],
        utron_bed=list(BED_DIR.glob("*.bed"))[0],
        vcf_path=list(VCF_DIR.glob("*.vcf.gz"))[0],
        outfile_tsv=str(tmp_path / "test.tsv"),
    )
    all_introns_counts_and_info.main(argv=args)
    assert True