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
BAM_SORTED_ASSIGNED_DIR = BAM_DIR / "sorted_assigned"

# pylint: disable=redefined-outer-name


@pytest.fixture()
def bam_file1() -> AlignmentFile:
    """Load a bam file."""
    return io.load_file(BAM_SORTED_ASSIGNED_DIR / "d0_no4sU_filtered_remapped_sorted.sorted.assigned.bam")


@pytest.fixture()
def bam_file2() -> AlignmentFile:
    """Load a bam file."""
    return io.load_file(BAM_SORTED_ASSIGNED_DIR / "d0_0hr1_filtered_remapped_sorted.sorted.assigned.bam")


@pytest.fixture()
def bam_unaligned_file1() -> AlignmentFile:
    """Load an unsorted and unassigned ``.bam`` file (ignore the filename!)."""
    return io.load_file(BAM_DIR / "d0_no4sU_filtered_remapped_sorted.bam")


@pytest.fixture()
def aligned_segment_unassigned_28584(bam_unaligned_file1: AlignmentFile) -> AlignedSegment:
    """Extract an individual AlignedSegment from a ``.bam`` file."""
    # NB : .fetch() returns AlignedSegments that span the start and end region, not just those within
    return next(bam_unaligned_file1.fetch(contig="chr9", start=28592, end=28593))


@pytest.fixture()
def aligned_segment_unassigned_17416(bam_unaligned_file1: AlignmentFile) -> AlignedSegment:
    """Extract an individual AlignedSegment from a ``.bam`` file."""
    # NB : .fetch() returns AlignedSegments that span the start and end region, not just those within
    return next(bam_unaligned_file1.fetch(contig="chr9", start=17804, end=18126))


@pytest.fixture()
def aligned_segment_unassigned_18029(bam_unaligned_file1: AlignmentFile) -> AlignedSegment:
    """Extract an individual AlignedSegment from a ``.bam`` file."""
    # NB : .fetch() returns AlignedSegments that span the start and end region, not just those within
    return next(bam_unaligned_file1.fetch(contig="chr9", start=18156, end=24870))


@pytest.fixture()
def bam_sorted_assigned_file() -> list[AlignedSegment]:
    """Return a list of a ``AlignedSegment`` from a ``.bam`` file that has been sorted and assigned."""
    return list(io.load_file(BAM_SORTED_ASSIGNED_DIR / "d0_no4sU_filtered_remapped_sorted.sorted.assigned.bam"))


@pytest.fixture()
def aligned_segment_assigned_15967(bam_sorted_assigned_file: list[AlignedSegment]) -> AlignedSegment:
    """Return a single assigned AlignedSegment where start is 15967."""
    return [x for x in bam_sorted_assigned_file if x.reference_start == 15967][0]


@pytest.fixture()
def aligned_segment_assigned_14770(bam_sorted_assigned_file: list[AlignedSegment]) -> AlignedSegment:
    """Return a single assigned AlignedSegment where start is 14770."""
    return [x for x in bam_sorted_assigned_file if x.reference_start == 14770][0]


@pytest.fixture()
def aligned_segment_assigned_17814(bam_sorted_assigned_file: list[AlignedSegment]) -> AlignedSegment:
    """Return a single assigned AlignedSegment where start is 17814."""
    return [x for x in bam_sorted_assigned_file if x.reference_start == 17814][0]
