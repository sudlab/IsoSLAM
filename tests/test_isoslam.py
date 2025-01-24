"""Tests for the isoslam module."""

from pathlib import Path
from types import GeneratorType
from typing import Any

import pytest  # type: ignore[import-not-found]

from isoslam import isoslam

BASE_DIR = Path.cwd()
RESOURCES = BASE_DIR / "tests" / "resources"
BAM_DIR = RESOURCES / "bam"
BED_DIR = RESOURCES / "bed"
GTF_DIR = RESOURCES / "gtf"
VCF_DIR = RESOURCES / "vcf"
BAM_SORTED_ASSIGNED_DIR = BAM_DIR / "sorted_assigned"

# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments


@pytest.mark.parametrize(
    ("bed_file", "expected_transcript"),
    [
        pytest.param(  # type: ignore[misc]
            BED_DIR / "test_coding_introns.bed",
            {
                "ENST00000442898": [
                    ("9", 14940, 15080, "ENST00000442898", "-"),
                    ("9", 15149, 15908, "ENST00000442898", "-"),
                    ("9", 16061, 16717, "ENST00000442898", "-"),
                    ("9", 16876, 16964, "ENST00000442898", "-"),
                    ("9", 17166, 17343, "ENST00000442898", "-"),
                    ("9", 17479, 17718, "ENST00000442898", "-"),
                    ("9", 17855, 18027, "ENST00000442898", "-"),
                    ("9", 18174, 18380, "ENST00000442898", "-"),
                    ("9", 18492, 24850, "ENST00000442898", "-"),
                    ("9", 25004, 29601, "ENST00000442898", "-"),
                ]
            },
            id="bed coding introons",
        ),
    ],
)
def test_isoslam_extract_transcripts(
    bed_file: str | Path,
    expected_transcript: dict[Any, list[tuple[Any, int, int, Any, Any]]],
) -> None:
    """Test extraction of tanscript data from bed file using extract_transcripts()."""
    assert isoslam.extract_transcripts(bed_file) == expected_transcript


@pytest.mark.parametrize(
    ("gtf_file", "expected_strand", "expected_transcript"),
    [
        pytest.param(  # type: ignore[misc]
            GTF_DIR / "test_wash1.gtf",
            {"MSTRG.63147": "-"},
            {"MSTRG.63147": ["ENST00000442898"]},
            id="gtf file as Path",
        ),
    ],
)
def test_extract_strand_transcript(
    gtf_file: str | Path, expected_strand: dict[Any, Any], expected_transcript: dict[Any, Any]
) -> None:
    """Test extraction of strand and transcript from gtf file using extract_strand_transcript()."""
    strand, transcript = isoslam.extract_strand_transcript(gtf_file)
    assert strand == expected_strand
    assert transcript == expected_transcript


@pytest.mark.parametrize(
    ("bam_file", "expected_length"),
    [
        pytest.param(  # type: ignore[misc]
            BAM_SORTED_ASSIGNED_DIR / "d0_no4sU_filtered_remapped_sorted.sorted.assigned.bam",
            710,
            id="file 1",
        ),
        pytest.param(
            BAM_SORTED_ASSIGNED_DIR / "d0_0hr1_filtered_remapped_sorted.sorted.assigned.bam",
            1845,
            id="file 2",
        ),
    ],
)
def test_extract_segment_pairs(bam_file: str | Path, expected_length: int) -> None:
    """Test extraction of paired segments from a sorted and assigned ``.bam`` file."""
    alignment_file = isoslam.extract_segment_pairs(bam_file)
    assert isinstance(alignment_file, GeneratorType)
    assert len(list(alignment_file)) == expected_length


@pytest.mark.parametrize(
    ("aligned_segment", "start", "end", "length", "status", "transcript", "block_start", "block_end", "reverse"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_unassigned_28584",
            28584,
            28733,
            150,
            None,
            None,
            (28584, 28704),
            (28704, 28733),
            True,
            id="28584 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_unassigned_17416",
            17416,
            17805,
            150,
            None,
            None,
            (17416, 17718),
            (17479, 17805),
            False,
            id="17416 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_unassigned_18029",
            18029,
            18385,
            150,
            None,
            None,
            (18029, 18380),
            (18174, 18385),
            True,
            id="18029 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_assigned_17814",
            17814,
            18136,
            150,
            "Assigned",
            "MSTRG.63147",
            (17814, 18027),
            (17855, 18136),
            True,
            id="17814 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_14770",
            14770,
            14876,
            150,
            "Assigned",
            "MSTRG.63147",
            (14770,),
            (14876,),
            False,
            id="14770 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_15967",
            15967,
            16117,
            150,
            "Assigned",
            "MSTRG.63147",
            (15967,),
            (16117,),
            False,
            id="15967 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_21051",
            21051,
            21201,
            150,
            "Unassigned_NoFeatures",
            None,
            (21051,),
            (21201,),
            True,
            id="21051 - Unassigned",
        ),
        pytest.param(
            "aligned_segment_assigned_20906",
            20906,
            21056,
            150,
            "Unassigned_NoFeatures",
            None,
            (20906,),
            (21056,),
            False,
            id="20906 - Unassigned",
        ),
    ],
)
def test_extract_features_from_read(
    aligned_segment: str,
    start: int,
    end: int,
    length: int,
    status: str,
    transcript: str,
    block_start: tuple[int, int],
    block_end: tuple[int, int],
    reverse: bool,
    request: pytest.FixtureRequest,
) -> None:
    """Test extract of features from an unassigned and assigned segment reads."""
    segment = isoslam.extract_features_from_read(request.getfixturevalue(aligned_segment))
    assert isinstance(segment, dict)
    assert segment["start"] == start
    assert segment["end"] == end
    assert segment["length"] == length
    assert segment["status"] == status
    assert segment["transcript"] == transcript
    assert segment["block_start"] == block_start
    assert segment["block_end"] == block_end
    assert segment["reverse"] == reverse


@pytest.mark.parametrize(
    ("aligned_segment1", "aligned_segment2", "expected"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_assigned_17814",
            "aligned_segment_assigned_14770",
            {
                "read1": {
                    "start": 17814,
                    "end": 18136,
                    "length": 150,
                    "status": "Assigned",
                    "transcript": "MSTRG.63147",
                    "block_start": (17814, 18027),
                    "block_end": (17855, 18136),
                    "reverse": True,
                },
                "read2": {
                    "start": 14770,
                    "end": 14876,
                    "length": 150,
                    "status": "Assigned",
                    "transcript": "MSTRG.63147",
                    "block_start": (14770,),
                    "block_end": (14876,),
                    "reverse": False,
                },
            },
            id="28584 and 17416 - Assigned and transcript MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_21051",
            "aligned_segment_assigned_20906",
            {
                "read1": {
                    "start": 21051,
                    "end": 21201,
                    "length": 150,
                    "status": "Unassigned_NoFeatures",
                    "transcript": None,
                    "block_start": (21051,),
                    "block_end": (21201,),
                    "reverse": True,
                },
                "read2": {
                    "start": 20906,
                    "end": 21056,
                    "length": 150,
                    "status": "Unassigned_NoFeatures",
                    "transcript": None,
                    "block_start": (20906,),
                    "block_end": (21056,),
                    "reverse": False,
                },
            },
            id="21051 and 21056 - Assignment and Transcript are None",
        ),
    ],
)
def test_extract_features_from_pair(
    aligned_segment1: str,
    aligned_segment2: str,
    expected: dict[str, dict[str, int | None | str | tuple[int, int]]],
    request: pytest.FixtureRequest,
) -> None:
    """Test extract of features from a list of pairs."""
    read_pair = isoslam.extract_features_from_pair(
        [
            request.getfixturevalue(aligned_segment1),
            request.getfixturevalue(aligned_segment2),
        ]
    )
    assert isinstance(read_pair, dict)
    assert read_pair == expected


@pytest.mark.parametrize(
    ("aligned_segment", "transcript_id", "length"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_unassigned_28584",
            "",
            0,
            id="28584 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_assigned_17814",
            "ENST00000442898",
            10,
            id="17814 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_14770",
            "ENST00000442898",
            10,
            id="14770 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_15967",
            "ENST00000442898",
            10,
            id="15967 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_21051",
            None,
            0,
            id="21051 - Unassigned",
        ),
        pytest.param(
            "aligned_segment_assigned_20906",
            None,
            0,
            id="20906 - Unassigned",
        ),
    ],
)
def test_extract_utron(
    aligned_segment: str,
    transcript_id: str,
    length: int,
    extract_transcript: dict[str, list[int | str]],
    extract_strand_transcript: tuple[dict[str, list[Any]], dict[str, str]],
    request: pytest.FixtureRequest,
) -> None:
    """Test extraction of the untranslated regions using extract_utron()."""
    segment = isoslam.extract_features_from_read(request.getfixturevalue(aligned_segment))
    _, gene_transcript = extract_strand_transcript
    untranslated_region = isoslam.extract_utron(segment, gene_transcript, coordinates=extract_transcript)
    assert isinstance(untranslated_region, list)
    assert len(untranslated_region) == length
    if len(untranslated_region):
        assert untranslated_region[0][3] == transcript_id  # type: ignore[misc]


@pytest.mark.parametrize(
    ("aligned_segment", "expected_start", "expected_end"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_unassigned_28584",
            (28584, 28704),
            (28704, 28733),
            id="28584 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_assigned_17814",
            (17814, 18027),
            (17855, 18136),
            id="17814 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_14770",
            (14770,),
            (14876,),
            id="14770 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_15967",
            (15967,),
            (16117,),
            id="15967 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_21051",
            (21051,),
            (21201,),
            id="21051 - Unassigned",
        ),
        pytest.param(
            "aligned_segment_assigned_20906",
            (20906,),
            (21056,),
            id="20906 - Unassigned",
        ),
    ],
)
def test_zip_blocks(
    aligned_segment: str,
    expected_start: tuple[int],
    expected_end: tuple[int],
    request: pytest.FixtureRequest,
) -> None:
    """Test that block starts and ends are zipped correctly."""
    start_blocks, end_blocks = isoslam.zip_blocks(request.getfixturevalue(aligned_segment))
    assert start_blocks == expected_start
    assert end_blocks == expected_end


@pytest.mark.parametrize(
    (
        "feature_pair",
        "read",
        "expected",
    ),
    [
        pytest.param(  # type: ignore[misc]
            "feature_pair_within_15967_24715",
            "read1",
            {"ENST00000442898": [(16061, 16717, "9", "-")]},
            id="15967 and 24715 feature pair from no4sU",
        ),
        pytest.param(
            "feature_pair_within_17790_18093",
            "read1",
            {"ENST00000442898": [(17855, 18027, "9", "-")]},
            id="17790 and 18093 feature pair from 0hr1",
        ),
        pytest.param(
            "feature_pair_within_17709_18290",
            "read1",
            {"ENST00000442898": [(17479, 17718, "9", "-"), (17855, 18027, "9", "-")]},
            id="17709 and 18290 feature pair from 0hr1",
        ),
        pytest.param(
            "feature_pair_within_14770_17814",
            "read1",
            {},
            id="14770 and 17814 feature pair from no4sU",
        ),
    ],
)
def test_filter_within_introns(
    feature_pair: tuple[dict["str", Any], list[int], list[int], list[int], list[int]],
    read: str,
    expected: dict[str, tuple[Any]],
    request: pytest.FixtureRequest,
) -> None:
    """
    Test filtering within introns.

    Tests are very basic and not all conditional statements are covered, in particular instances where transcripts span
    the whole region are not covered. This may be challenging as the length is _always_ 150.
    """
    # We first extract the required pair_features and block start/ends from the fixture
    pair_features, blocks = request.getfixturevalue(feature_pair)
    within_introns = isoslam.filter_within_introns(
        pair_features,
        blocks,
        read,
    )
    assert within_introns == expected


@pytest.mark.parametrize(
    (
        "feature_pair",
        "read",
        "expected",
    ),
    [
        pytest.param(  # type: ignore[misc]
            "feature_pair_spliced_17739_17814",
            "read1",
            {"ENST00000442898": [(17855, 18027, "9", "-")]},
            id="14770 and 17814 for read1 from no4sU",
        ),
        pytest.param(
            "feature_pair_spliced_17430_18155",
            "read2",
            {"ENST00000442898": [(18174, 18380, "9", "-"), (18492, 24850, "9", "-")]},
            id="17430 and 18155 for from 0hr1",
        ),
        pytest.param(
            "feature_pair_within_14770_17814",
            "read2",
            {"ENST00000442898": [(17855, 18027, "9", "-")]},
            id="17430 and 18155 for from no4sU",
        ),
    ],
)
def test_filter_spliced_utrons(
    feature_pair: tuple[dict["str", Any], list[int], list[int], list[int], list[int]],
    read: str,
    expected: dict[str, list[Any]],
    request: pytest.FixtureRequest,
) -> None:
    """Test filtering spliced utrons.

    Tests are very basic and not all conditional statements are covered.
    """
    # We first extract the required pair_features and block start/ends from the fixture
    pair_features, blocks = request.getfixturevalue(feature_pair)
    within_introns = isoslam.filter_spliced_utrons(
        pair_features,
        blocks,
        read,
    )
    assert within_introns == expected


@pytest.mark.parametrize(
    ("reads1", "reads2", "expected"),
    [
        pytest.param(  # type: ignore[misc]
            {
                "transcript1": [
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                ]
            },
            {
                "transcript1": [
                    (1, 2, "9", "Ret"),
                    (5, 6, "11", "Ret"),
                ]
            },
            {
                ("transcript1", (1, 2, "9", "Ret")),
                ("transcript1", (3, 4, "10", "Ret")),
                ("transcript1", (5, 6, "11", "Ret")),
            },
            id="one transcripts across two reads with one duplicate",
        ),
        pytest.param(
            {
                "transcript1": [
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                ],
                "transcript2": [
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                ],
            },
            {
                "transcript1": [
                    (1, 2, "9", "Ret"),
                    (5, 6, "11", "Ret"),
                ],
                "transcript2": [
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                ],
            },
            {
                ("transcript1", (1, 2, "9", "Ret")),
                ("transcript1", (3, 4, "10", "Ret")),
                ("transcript1", (5, 6, "11", "Ret")),
                ("transcript2", (5, 6, "20", "Ret")),
                ("transcript2", (5, 6, "20", "Spl")),
                ("transcript2", (7, 8, "20", "Ret")),
                ("transcript2", (9, 10, "21", "Ret")),
            },
            id="two transcripts across two reads with one duplicate one very similar",
        ),
        pytest.param(
            {
                "transcript1": [
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                    (1, 2, "9", "Ret"),
                    (3, 4, "10", "Ret"),
                ],
                "transcript2": [
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                ],
            },
            {
                "transcript1": [
                    (1, 2, "9", "Ret"),
                    (5, 6, "11", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                    (5, 6, "20", "Ret"),
                    (7, 8, "20", "Ret"),
                ],
                "transcript2": [
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                    (5, 6, "20", "Spl"),
                    (9, 10, "21", "Ret"),
                ],
            },
            {
                ("transcript1", (1, 2, "9", "Ret")),
                ("transcript1", (3, 4, "10", "Ret")),
                ("transcript1", (5, 6, "11", "Ret")),
                ("transcript1", (5, 6, "20", "Ret")),
                ("transcript1", (7, 8, "20", "Ret")),
                ("transcript2", (5, 6, "20", "Ret")),
                ("transcript2", (5, 6, "20", "Spl")),
                ("transcript2", (7, 8, "20", "Ret")),
                ("transcript2", (9, 10, "21", "Ret")),
            },
            id="BIG dictionaries",
        ),
    ],
)
def test_unique_conversions(
    reads1: dict[str, list[tuple[Any]]], reads2: dict[str, list[tuple[Any]]], expected: set[list[Any]]
) -> None:
    """Test the combining of two read sets into one unique set."""
    unique_conversions = isoslam.unique_conversions(reads1, reads2)
    assert unique_conversions == expected


@pytest.mark.parametrize(
    ("set1", "set2", "expected_set1", "expected_set2"),
    [
        pytest.param(  # type: ignore[misc]
            {"a", "b", "c"},
            {"d", "e", "f"},
            {"a", "b", "c"},
            {"d", "e", "f"},
            id="No overlap",
        ),
        pytest.param({"a", "b", "c"}, {"b", "d", "e", "f"}, {"a", "c"}, {"d", "e", "f"}, id="b common"),
        pytest.param({"a", "b", "c"}, {"b", "c", "d", "e", "f"}, {"a"}, {"d", "e", "f"}, id="b and ccommon"),
        pytest.param({"a", "b", "c"}, {"a", "b", "c"}, set(), set(), id="identical"),
    ],
)
def test_remove_common_reads(set1: set[str], set2: set[str], expected_set1: set[str], expected_set2: set[str]) -> None:
    """Test common elements are removed from two sets."""
    set1, set2 = isoslam.remove_common_reads(set1, set2)  # type: ignore[arg-type, assignment]
    assert set1 == expected_set1
    assert set2 == expected_set2


@pytest.mark.parametrize(
    (
        "read",
        "vcf_file_fixture",
        "conversion_from",
        "conversion_to",
        "len_convertible",
        "len_converted_position",
        "len_coverage",
    ),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_assigned_14770", "vcf_file", "T", "C", 23, 0, 106, id="no4sU 14770 T>C"
        ),
        pytest.param("aligned_segment_assigned_14770", "vcf_file", "A", "G", 14, 0, 106, id="no4sU 14770 A>G"),
        pytest.param("aligned_segment_assigned_21051", "vcf_file", "T", "C", 31, 0, 150, id="0hr1 21051 T>C"),
        pytest.param("aligned_segment_assigned_21051", "vcf_file", "A", "G", 23, 1, 150, id="0hr1 21051 A>G"),
        pytest.param("aligned_segment_assigned_20906", "vcf_file", "T", "C", 17, 0, 150, id="0hr1 20906 T>C"),
        pytest.param("aligned_segment_assigned_20906", "vcf_file", "A", "G", 50, 1, 150, id="0hr1 20906 A>G"),
    ],
)
def test_conversions_per_read(
    read: str,
    vcf_file_fixture: str,
    conversion_from: str,
    conversion_to: str,
    len_convertible: int,
    len_converted_position: int,
    len_coverage: int,
    request: pytest.FixtureRequest,
) -> None:
    """Test the building of sets of conversion for a given read."""
    convertible, converted_position, coverage = isoslam.conversions_per_read(
        request.getfixturevalue(read),
        conversion_from,
        conversion_to,
        convertible=set(),
        converted_position=set(),
        coverage=set(),
        vcf_file=request.getfixturevalue(vcf_file_fixture),
    )
    assert isinstance(convertible, set)
    assert isinstance(converted_position, set)
    assert isinstance(coverage, set)
    assert len(convertible) == len_convertible
    assert len(converted_position) == len_converted_position
    assert len(coverage) == len_coverage


@pytest.mark.parametrize(
    ("forward_read", "reverse_read", "forward_conversion", "reverse_conversion", "vcf_file_fixture", "expected"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_assigned_14770",
            "aligned_segment_assigned_17814",
            {"from": "A", "to": "G"},
            {"from": "T", "to": "C"},
            "vcf_file",
            {"convertible": 46, "converted_position": 0, "coverage": 256},
            id="no4sU 14770/17814",
        ),
        pytest.param(
            "aligned_segment_assigned_14770",
            "aligned_segment_assigned_14770",
            {"from": "A", "to": "G"},
            {"from": "T", "to": "C"},
            "vcf_file",
            {"convertible": 37, "converted_position": 0, "coverage": 106},
            id="no4sU 14770/14770 (checks sets are unique)",
        ),
        pytest.param(
            "aligned_segment_assigned_21051",
            "aligned_segment_assigned_20906",
            {"from": "A", "to": "G"},
            {"from": "T", "to": "C"},
            "vcf_file",
            {"convertible": 40, "converted_position": 1, "coverage": 295},
            id="0hr1 21051/20906",
        ),
    ],
)
def test_count_conversions_across_pairs(
    forward_read: str,
    reverse_read: str,
    forward_conversion: dict[str, str],
    reverse_conversion: dict[str, str],
    vcf_file_fixture: str,
    expected: dict[str, int],
    request: pytest.FixtureRequest,
) -> None:
    """Test that pairs of conversions are counted correctly."""
    counts = isoslam.count_conversions_across_pairs(
        forward_read=request.getfixturevalue(forward_read),
        reverse_read=request.getfixturevalue(reverse_read),
        vcf_file=request.getfixturevalue(vcf_file_fixture),
        forward_conversion=forward_conversion,
        reverse_conversion=reverse_conversion,
    )
    for count in ["convertible", "converted_position", "coverage"]:
        assert counts[count] == expected[count]
