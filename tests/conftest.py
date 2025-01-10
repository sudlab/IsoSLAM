"""Fixtures for pytest."""

from pathlib import Path
from typing import Any

import pysam
import pytest
from pysam import AlignedSegment, AlignmentFile

from isoslam import io, isoslam

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
def gtf_file() -> pysam.tabix_generic_iterator:
    """Load a ``.gtf`` file."""
    return io.load_file(GTF_DIR / "test_wash1.gtf")


@pytest.fixture()
def bed_file() -> dict:
    """Load a ``.gtf`` file."""
    return io.load_file(BED_DIR / "test_coding_introns.bed")


@pytest.fixture()
def extract_strand_transcript() -> tuple[dict, dict]:
    """Tuple of dictionaries from ``.gtf`` file."""
    return isoslam.extract_strand_transcript(GTF_DIR / "test_wash1.gtf")


@pytest.fixture()
def extract_transcript() -> dict:
    """Extract dictionary of transcripts from ``.bed`` file."""
    return isoslam.extract_transcripts(BED_DIR / "test_coding_introns.bed")


@pytest.fixture()
def default_config() -> dict[str:Any]:
    """Load the default configuration."""
    return io.load_and_update_config(args=None)


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
def bam_sorted_assigned_file_no4sU() -> list[AlignedSegment]:
    """Return a list of a ``AlignedSegment`` from a ``.bam`` file that has been sorted and assigned."""
    return list(io.load_file(BAM_SORTED_ASSIGNED_DIR / "d0_no4sU_filtered_remapped_sorted.sorted.assigned.bam"))


@pytest.fixture()
def aligned_segment_assigned_15967(bam_sorted_assigned_file_no4sU: list[AlignedSegment]) -> AlignedSegment:
    """Return a single assigned AlignedSegment where start is 15967."""
    return [x for x in bam_sorted_assigned_file_no4sU if x.reference_start == 15967][0]


@pytest.fixture()
def aligned_segment_assigned_14770(bam_sorted_assigned_file_no4sU: list[AlignedSegment]) -> AlignedSegment:
    """Return a single assigned AlignedSegment where start is 14770."""
    return [x for x in bam_sorted_assigned_file_no4sU if x.reference_start == 14770][0]


@pytest.fixture()
def aligned_segment_assigned_17814(bam_sorted_assigned_file_no4sU: list[AlignedSegment]) -> AlignedSegment:
    """Return a single assigned AlignedSegment where start is 17814."""
    return [x for x in bam_sorted_assigned_file_no4sU if x.reference_start == 17814][0]


@pytest.fixture()
def bam_sorted_assigned_file_0hr1() -> list[AlignedSegment]:
    """Return a list of a ``AlignedSegment`` from a ``.bam`` file that has been sorted and assigned."""
    return list(io.load_file(BAM_SORTED_ASSIGNED_DIR / "d0_0hr1_filtered_remapped_sorted.sorted.assigned.bam"))


@pytest.fixture()
def aligned_segment_assigned_21051(bam_sorted_assigned_file_0hr1: list[AlignedSegment]) -> AlignedSegment:
    """
    Return a single assigned AlignedSegment where start is 21051.

    Along with the next fixture (20906) this was the 484th pair of reads that caused some early problems with regression
    tests which is why they have been pulled out and added to the parameterised tests.

    The reason these problems occurred is because initially 9 alignments are found at this location but 6 of these
    should be filtered out subsequently, a mistake in refactoring used the wrong length. But this will serve as a useful
    aligned segment to be testing because of this filtering.
    """
    return [x for x in bam_sorted_assigned_file_0hr1 if x.reference_start == 21051][0]


@pytest.fixture()
def aligned_segment_assigned_20906(bam_sorted_assigned_file_0hr1: list[AlignedSegment]) -> AlignedSegment:
    """
    Return a single assigned AlignedSegment where start is 20906.

    Along with the previous fixture (21051) this was the 484th pair of reads that cause some early problems with
    regression tests which is why they have been pulled out and added to the parameterised tests.

    The reason these problems occurred is because initially 9 alignments are found at this location but 6 of these
    should be filtered out subsequently, a mistake in refactoring used the wrong length. But this will serve as a useful
    aligned segment to be testing because of this filtering.
    """
    return [x for x in bam_sorted_assigned_file_0hr1 if x.reference_start == 20906][0]
