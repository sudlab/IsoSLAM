"""IsoSLAM module."""

from collections import defaultdict
from collections.abc import Generator
from pathlib import Path
from typing import Any

from loguru import logger
from pysam import AlignedSegment

from isoslam import io


def extract_transcripts(bed_file: str | Path) -> dict[Any, list[tuple[Any, int, int, Any, Any]]]:
    """
    Extract features from ``.bed`` file and return as a dictionary indexed by ``transcript_id``.

    Parameters
    ----------
    bed_file : str | Path
        Path, as string or pathlib Path, to a ``.bed`` file.

    Returns
    -------
    dict[Any, list[tuple[Any, int, int, Any, Any]]]
        Dictionary of ``chromosome``, ``start``, ``end``, ``transcript_id`` and ``bedstrand`` indexed by
        ``transcript_id``.
    """
    coordinates = defaultdict(list)
    for line in io.load_file(bed_file):
        contents = line.strip().split("\t")
        transcript_id = contents[3].replace("_intron", "")
        coordinates[transcript_id].append(
            (
                contents[0],
                int(contents[1]),
                int(contents[2]),
                transcript_id,
                contents[5],
            )
        )
    logger.info(f"Extracted features from : {bed_file}")
    return coordinates


def extract_strand_transcript(gtf_file: str | Path) -> tuple[defaultdict[Any, Any], defaultdict[Any, list[Any]]]:
    """
    Extract strand and transcript ID data from ``.gtf`` file.

    Parameters
    ----------
    gtf_file : Path | str
        Path to a 'gtf' file.

    Returns
    -------
    tuple[dict[str, tuple[str]], dict[str, tuple[str]]]
        Two dictionaries are returned, one of the ``strand`` the other of the ``transcript_id`` both using the
        ``gene_id`` as the key.
    """
    strand = defaultdict(str)
    transcript = defaultdict(list)
    for entry in io.load_file(gtf_file):
        if not entry.feature == "transcript":
            continue
        strand[entry.gene_id] = entry.strand
        transcript[entry.gene_id].append(entry.transcript_id)
    logger.info(f"Extracted features from : {gtf_file}")
    return (strand, transcript)


def extract_segment_pairs(bam_file: str | Path) -> Generator[AlignedSegment]:
    """
    Extract pairs of AlignedSegments from a ``.bam`` file.

    When there are two adjacent ``AlignedSegments`` with the same ``query_name`` only the first is paired, subsequent
    segments are dropped.

    Parameters
    ----------
    bam_file : str | Path
        Path to a ``.bam`` file.

    Yields
    ------
    Generator
        Itterable of paired segments.
    """
    previous_read: str | None = None
    pair: list[AlignedSegment] = []
    for read in io.load_file(bam_file):
        # Return pairs of reads, i.e. not on first pass, nor if query_name matches the previous read
        if previous_read is not None and previous_read != read.query_name:
            yield pair
            pair = []
            previous_read = read.query_name
        previous_read = read.query_name
        pair.append(read)
    # Don't forget to return the last pair!
    yield pair


def extract_features_from_read(read: AlignedSegment) -> dict[str, int | str | None | tuple[int, int]]:
    """
    Extract start, end and length from an aligned segment read.

    Parameters
    ----------
    read : AlignedSegment
        An aligned segment read from ``pysam``.

    Returns
    -------
    dict[str, Any]
        Dictionary of ``start``, ``end`` and ``length`` of the segment.
    """
    block_start, block_end = zip(*read.get_blocks())
    try:
        status = read.get_tag("XS")
    except KeyError:
        status = None
    try:
        transcript = read.get_tag("XT")
    except KeyError:
        transcript = None
    return {
        "start": read.reference_start,
        "end": read.reference_end,
        "length": read.query_length,
        "status": status,
        "transcript": transcript,
        "block_start": block_start,
        "block_end": block_end,
    }


def extract_features_from_pair(pair: list[AlignedSegment]) -> dict[str, dict[str, Any]]:
    """
    Extract features from a pair of reads.

    Parameters
    ----------
    pair : list[AlignedSegment]
        A list of two aligned segments from ``pysam``.

    Returns
    -------
    dic[str, dict[str, Any]]
        Returns a nested dictionaries of the ``start``, ``end`` and ``length`` of each read.
    """
    return {
        "read1": extract_features_from_read(pair[0]),
        "read2": extract_features_from_read(pair[1]),
    }


def extract_utron(features: dict[str, Any], gene_transcript: Any, coordinates: Any) -> list[tuple[int | str]] | None:
    """
    Extract and sum the utrons based on tag.

    ACTION : This function needs better documentation, my guess is that its extracting the transcripts to genes and
    then getting some related information (what I'm not sure) from the .bed file and adding these up.

    Parameters
    ----------
    features : str
        A tag from an assigned read.
    gene_transcript : TextIO
        Transcript to gene from a ``.gtf`` file.
    coordinates : Any
        Untranslated region coordinates from a ``.bed`` file.

    Returns
    -------
    list | None
        List of the length of assigned regions.
    """
    if features["status"] == "Assigned":
        untranslated_regions = [coordinates[transcript] for transcript in gene_transcript[features["transcript"]]]
        return sum(untranslated_regions, [])
    return []
