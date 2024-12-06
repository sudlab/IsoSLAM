"""Fixtures for pytest."""

from pathlib import Path

import pytest
from pysam import AlignedSegment, AlignmentFile

from isoslam import io

BASE_DIR = Path.cwd()
TEST_DIR = BASE_DIR / "tests"
RESOURCES = TEST_DIR / "resources"
GTF_DIR = RESOURCES / "gtf"
BED_DIR = RESOURCES / "bed"
VCF_DIR = RESOURCES / "vcf"
BAM_DIR = RESOURCES / "bam"


# pylint: disable=redefined-outer-name


@pytest.fixture()
def bam_file() -> AlignmentFile:
    """Load a bam file."""
    return io.load_file(BAM_DIR / "d0_no4sU_filtered_remapped_sorted.bam")


@pytest.fixture()
def aligned_segment_28584(bam_file: AlignmentFile) -> AlignedSegment:
    """Extract an individual AlignedSegment from a ``.bam`` file."""
    return next(bam_file.fetch(contig="chr9", start=28592, end=28593))


@pytest.fixture()
def aligned_segment_17416(bam_file: AlignmentFile) -> AlignedSegment:
    """Extract an individual AlignedSegment from a ``.bam`` file."""
    # @ns-rse : I have no idea why the generator doesn't return AlignedSegments that match the start/end here
    return next(bam_file.fetch(contig="chr9", start=17804, end=18126))


@pytest.fixture()
def aligned_segment_18029(bam_file: AlignmentFile) -> AlignedSegment:
    """Extract an individual AlignedSegment from a ``.bam`` file."""
    # @ns-rse : I have no idea why the generator doesn't return AlignedSegments that match the start/end here
    return next(bam_file.fetch(contig="chr9", start=18156, end=24870))
