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


<<<<<<< HEAD
<<<<<<< HEAD
@pytest.mark.parametrize(
    ("aligned_segment1", "aligned_segment2", "expected"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_28584",
            "aligned_segment_17416",
            {
                "read1": {
                    "start": 28584,
                    "end": 28733,
                    "length": 149,
                    "status": None,
                    "transcript": None,
                    "block_start": (28584, 28704),
                    "block_end": (28704, 28733),
                },
                "read2": {
                    "start": 17416,
                    "end": 17805,
                    "length": 389,
                    "status": None,
                    "transcript": None,
                    "block_start": (17416, 17718),
                    "block_end": (17479, 17805),
                },
            },
            id="28584 and 17416 - Assignment and Transcript are None",
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
||||||| parent of 3878e7a (fixup! feature(isoslam): Adds functions to extract features from reads)
# @pytest.mark.parametrize(
#     ("aligned_segment", "start", "end", "length", "status", "transcript", "block_start", "block_end"),
#     [
#         pytest.param(  # type: ignore[misc]
#             "aligned_segment_28584",
#             28584,
#             28733,
#             149,
#             None,
#             None,
#             (28584, 28704),
#             (28704, 28733),
#             id="28584 - Assignment and Transcript are None",
#         ),
#         pytest.param(
#             "aligned_segment_17416",
#             17416,
#             17805,
#             389,
#             "Assigned",
#             "Something",
#             (17416, 17718),
#             (17479, 17805),
#             id="17416 - Assignment and Transcript are None",
#             marks=pytest.mark.xfail(reasion="Not getting Assigned or Transcript tags for some reason?"),
#         ),
#         pytest.param(
#             "aligned_segment_18029",
#             18029,
#             18385,
#             356,
#             "Assigned",
#             "something",
#             (18029, 18380),
#             (18174, 18385),
#             id="17416 - Assignment and Transcript are None",
#             marks=pytest.mark.xfail(reasion="Not getting Assigned tag for some reason?"),
#         ),
#     ],
# )
# def test_extract_features_from_read_param(
#     aligned_segment: str,
#     start: int,
#     end: int,
#     length: int,
#     status: str,
#     transcript: str,
#     block_start: tuple[int, int],
#     block_end: tuple[int, int],
#     request: pytest.FixtureRequest,
# ) -> None:
#     """Test extract of features from an aligned segment read."""
#     segment = isoslam.extract_features_from_read(request.getfixturevalue(aligned_segment))
#     assert isinstance(segment, dict)
#     assert segment["start"] == start
#     assert segment["end"] == end
#     assert segment["length"] == length
#     assert segment["status"] == status
#     assert segment["transcript"] == transcript
#     assert segment["block_start"] == block_start
#     assert segment["block_end"] == block_end


# @pytest.mark.parametrize(
#     ("aligned_segment1", "aligned_segment2", "expected"),
#     [
#         pytest.param(  # type: ignore[misc]
#             "aligned_segment_28584",
#             "aligned_segment_17416",
#             {
#                 "read1": {
#                     "start": 28584,
#                     "end": 28733,
#                     "length": 149,
#                     "status": None,
#                     "transcript": None,
#                     "block_start": (28584, 28704),
#                     "block_end": (28704, 28733),
#                 },
#                 "read2": {
#                     "start": 17416,
#                     "end": 17805,
#                     "length": 389,
#                     "status": None,
#                     "transcript": None,
#                     jls_extract_var: (17416, 17718),
#                     "block_end": (17479, 17805),
#                 },
#             },
#             id="28584 and 17416 - Assignment and Transcript are None",
#         ),
#     ],
# )
# def test_extract_features_from_pair(
#     aligned_segment1: str,
#     aligned_segment2: str,
#     expected: dict[str, dict[str, int | None | str | tuple[int, int]]],
#     request: pytest.FixtureRequest,
# ) -> None:
#     """Test extract of features from a list of pairs."""
#     read_pair = isoslam.extract_features_from_pair(
#         [
#             request.getfixturevalue(aligned_segment1),
#             request.getfixturevalue(aligned_segment2),
#         ]
#     )
#     assert isinstance(read_pair, dict)
#     assert read_pair == expected
=======
@pytest.mark.parametrize(
    ("aligned_segment", "start", "end", "length", "status", "transcript", "block_start", "block_end"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_unassigned_28584",
            28584,
            28733,
            149,
            None,
            None,
            (28584, 28704),
            (28704, 28733),
            id="28584 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_unassigned_17416",
            17416,
            17805,
            389,
            None,
            None,
            (17416, 17718),
            (17479, 17805),
            id="17416 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_unassigned_18029",
            18029,
            18385,
            356,
            None,
            None,
            (18029, 18380),
            (18174, 18385),
            id="18029 - Assignment and Transcript are None",
        ),
        pytest.param(
            "aligned_segment_assigned_17814",
            17814,
            18136,
            322,
            "Assigned",
            "MSTRG.63147",
            (17814, 18027),
            (17855, 18136),
            id="17814 - Assigned to MSTRG.63147",
        ),
        pytest.param(
            "aligned_segment_assigned_14770",
            14770,
            14876,
            106,
            "Assigned",
            "MSTRG.63147",
            (14770,),
            (14876,),
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
            id="15967 - Assigned to MSTRG.63147",
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
                    "length": 322,
                    "status": "Assigned",
                    "transcript": "MSTRG.63147",
                    "block_start": (17814, 18027),
                    "block_end": (17855, 18136),
                },
                "read2": {
                    "start": 14770,
                    "end": 14876,
                    "length": 106,
                    "status": "Assigned",
                    "transcript": "MSTRG.63147",
                    "block_start": (14770,),
                    "block_end": (14876,),
                },
            },
            id="28584 and 17416 - Assignment and Transcript are None",
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
>>>>>>> 3878e7a (fixup! feature(isoslam): Adds functions to extract features from reads)
||||||| parent of 3ee6aa7 (feature(isoslam): Adds extract_segment_pair() function)
@pytest.mark.parametrize(
    ("aligned_segment1", "aligned_segment2", "expected"),
    [
        pytest.param(  # type: ignore[misc]
            "aligned_segment_28584",
            "aligned_segment_17416",
            {
                "read1": {
                    "start": 28584,
                    "end": 28733,
                    "length": 149,
                    "status": None,
                    "transcript": None,
                    "block_start": (28584, 28704),
                    "block_end": (28704, 28733),
                },
                "read2": {
                    "start": 17416,
                    "end": 17805,
                    "length": 389,
                    "status": None,
                    "transcript": None,
                    "block_start": (17416, 17718),
                    "block_end": (17479, 17805),
                },
            },
            id="28584 and 17416 - Assignment and Transcript are None",
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
=======
# @pytest.mark.parametrize(
#     ("aligned_segment", "start", "end", "length", "status", "transcript", "block_start", "block_end"),
#     [
#         pytest.param(  # type: ignore[misc]
#             "aligned_segment_28584",
#             28584,
#             28733,
#             149,
#             None,
#             None,
#             (28584, 28704),
#             (28704, 28733),
#             id="28584 - Assignment and Transcript are None",
#         ),
#         pytest.param(
#             "aligned_segment_17416",
#             17416,
#             17805,
#             389,
#             "Assigned",
#             "Something",
#             (17416, 17718),
#             (17479, 17805),
#             id="17416 - Assignment and Transcript are None",
#             marks=pytest.mark.xfail(reasion="Not getting Assigned or Transcript tags for some reason?"),
#         ),
#         pytest.param(
#             "aligned_segment_18029",
#             18029,
#             18385,
#             356,
#             "Assigned",
#             "something",
#             (18029, 18380),
#             (18174, 18385),
#             id="17416 - Assignment and Transcript are None",
#             marks=pytest.mark.xfail(reasion="Not getting Assigned tag for some reason?"),
#         ),
#     ],
# )
# def test_extract_features_from_read_param(
#     aligned_segment: str,
#     start: int,
#     end: int,
#     length: int,
#     status: str,
#     transcript: str,
#     block_start: tuple[int, int],
#     block_end: tuple[int, int],
#     request: pytest.FixtureRequest,
# ) -> None:
#     """Test extract of features from an aligned segment read."""
#     segment = isoslam.extract_features_from_read(request.getfixturevalue(aligned_segment))
#     assert isinstance(segment, dict)
#     assert segment["start"] == start
#     assert segment["end"] == end
#     assert segment["length"] == length
#     assert segment["status"] == status
#     assert segment["transcript"] == transcript
#     assert segment["block_start"] == block_start
#     assert segment["block_end"] == block_end


# @pytest.mark.parametrize(
#     ("aligned_segment1", "aligned_segment2", "expected"),
#     [
#         pytest.param(  # type: ignore[misc]
#             "aligned_segment_28584",
#             "aligned_segment_17416",
#             {
#                 "read1": {
#                     "start": 28584,
#                     "end": 28733,
#                     "length": 149,
#                     "status": None,
#                     "transcript": None,
#                     "block_start": (28584, 28704),
#                     "block_end": (28704, 28733),
#                 },
#                 "read2": {
#                     "start": 17416,
#                     "end": 17805,
#                     "length": 389,
#                     "status": None,
#                     "transcript": None,
#                     jls_extract_var: (17416, 17718),
#                     "block_end": (17479, 17805),
#                 },
#             },
#             id="28584 and 17416 - Assignment and Transcript are None",
#         ),
#     ],
# )
# def test_extract_features_from_pair(
#     aligned_segment1: str,
#     aligned_segment2: str,
#     expected: dict[str, dict[str, int | None | str | tuple[int, int]]],
#     request: pytest.FixtureRequest,
# ) -> None:
#     """Test extract of features from a list of pairs."""
#     read_pair = isoslam.extract_features_from_pair(
#         [
#             request.getfixturevalue(aligned_segment1),
#             request.getfixturevalue(aligned_segment2),
#         ]
#     )
#     assert isinstance(read_pair, dict)
#     assert read_pair == expected
>>>>>>> 3ee6aa7 (feature(isoslam): Adds extract_segment_pair() function)
