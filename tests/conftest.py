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


@pytest.fixture()
def feature_pair_14770_17814(
    bam_sorted_assigned_file_no4sU: list[AlignedSegment],
    extract_transcript: dict[str, list[int | str]],
    extract_strand_transcript: tuple[dict[str, list[Any]], dict[str, str]],
) -> tuple[dict["str", Any], list, list, list, list]:
    """Feature pair for 17814 and 14770 augment with utrons and return along with block start/ends."""
    aligned_segment_14770 = [x for x in bam_sorted_assigned_file_no4sU if x.reference_start == 14770][0]
    aligned_segment_17814 = [x for x in bam_sorted_assigned_file_no4sU if x.reference_start == 17814][0]
    pair_features = isoslam.extract_features_from_pair(
        [
            aligned_segment_14770,
            aligned_segment_17814,
        ]
    )
    # ...and get the utron for both reads and add them to the dictionary
    _, gene_transcript = extract_strand_transcript
    pair_features["read1"]["utron"] = isoslam.extract_utron(
        pair_features["read1"], gene_transcript, coordinates=extract_transcript
    )
    pair_features["read2"]["utron"] = isoslam.extract_utron(
        pair_features["read2"], gene_transcript, coordinates=extract_transcript
    )
    # ...then we can get the block_starts/ends
    block_starts1, block_ends1 = isoslam.zip_blocks(aligned_segment_14770)
    block_starts2, block_ends2 = isoslam.zip_blocks(aligned_segment_17814)
    blocks = {
        "read1": {"starts": block_starts1, "ends": block_ends1},
        "read2": {"starts": block_starts2, "ends": block_ends2},
    }
    # return (pair_features, block_starts1, block_ends1, block_starts2, block_ends2)
    return (pair_features, blocks)


@pytest.fixture()
def feature_pair_20906_21051(
    bam_sorted_assigned_file_0hr1: list[AlignedSegment],
    extract_transcript: dict[str, list[int | str]],
    extract_strand_transcript: tuple[dict[str, list[Any]], dict[str, str]],
) -> tuple[dict["str", Any], list, list, list, list]:
    """Feature pair for 21051 and 20906 augment with utrons and return along with block start/ends."""
    aligned_segment_20906 = [x for x in bam_sorted_assigned_file_0hr1 if x.reference_start == 20906][0]
    aligned_segment_21051 = [x for x in bam_sorted_assigned_file_0hr1 if x.reference_start == 21051][0]
    pair_features = isoslam.extract_features_from_pair(
        [
            aligned_segment_20906,
            aligned_segment_21051,
        ]
    )
    # ...and get the utron for both reads and add them to the dictionary
    _, gene_transcript = extract_strand_transcript
    pair_features["read1"]["utron"] = isoslam.extract_utron(
        pair_features["read1"], gene_transcript, coordinates=extract_transcript
    )
    pair_features["read2"]["utron"] = isoslam.extract_utron(
        pair_features["read2"], gene_transcript, coordinates=extract_transcript
    )
    # ...then we can get the block_starts/ends
    block_starts1, block_ends1 = isoslam.zip_blocks(aligned_segment_20906)
    block_starts2, block_ends2 = isoslam.zip_blocks(aligned_segment_21051)
    blocks = {
        "read1": {"starts": block_starts1, "ends": block_ends1},
        "read2": {"starts": block_starts2, "ends": block_ends2},
    }
    # return (pair_features, block_starts1, block_ends1, block_starts2, block_ends2)
    return (pair_features, blocks)
