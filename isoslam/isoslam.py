"""IsoSLAM module."""

from collections import defaultdict
from pathlib import Path
from typing import Any

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
        Two dictionaries are returned, one of the ``strand`` the other of the ``transcript_id`` both using the ``gene_id`` as
        the key.
    """
    strand = defaultdict(str)
    transcript = defaultdict(list)
    for entry in io.load_file(gtf_file):
        if not entry.feature == "transcript":
            continue
        strand[entry.gene_id] = entry.strand
        transcript[entry.gene_id].append(entry.transcript_id)
    return (strand, transcript)
