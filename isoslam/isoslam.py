"""IsoSLAM module."""

from collections import defaultdict
from pathlib import Path
from typing import Any

from isoslam import io


def extract_transcripts(bed_file: str | Path) -> dict[Any, list[tuple[Any, int, int, Any, Any]]]:
    """
    Extract features from `.bed` file and return as a dictionary indexed by transcript_id.

    Parameters
    ----------
    bed_file : str | Path
        Path, as string or pathlib Path, to a `.bed` file.

    Returns
    -------
    dict[Any, list[tuple[Any, int, int, Any, Any]]]
        Nested dictionary of chromosome, start, end and bedstrand indexed by transcript_id.
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
