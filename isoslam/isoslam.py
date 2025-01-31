"""IsoSLAM module."""

from collections import defaultdict
from collections.abc import Generator, Iterator
from pathlib import Path
from typing import Any

import polars as pl
from loguru import logger
from pysam import AlignedSegment, VariantFile

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
    try:
        reverse = read.is_reverse
    except KeyError:
        reverse = None
    return {
        "start": read.reference_start,
        "end": read.reference_end,
        "length": read.query_length,
        "status": status,
        "transcript": transcript,
        "block_start": block_start,
        "block_end": block_end,
        "reverse": reverse,
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


def zip_blocks(read: AlignedSegment) -> Iterator[tuple[Any, ...]]:
    """
    Zip the block starts and ends into two lists.

    Parameters
    ----------
    read : AlignedSegment
        An individual aligned segment read from a ''.bam'' file.

    Returns
    -------
    tuple[list[int], list[int]]
        Tuple of two lists of integers the first is start location, the second is the end location.
    """
    return zip(*read.get_blocks())


def filter_within_introns(
    pair_features: dict[str, dict[str, Any]],
    blocks: dict[str, dict[str, set[int]]],
    read: str = "read1",
) -> dict[str, tuple[Any]]:
    """
    Filter utrons that are within introns.

    Parameters
    ----------
    pair_features : dict[str, dict]
        Dictionary of extracted features and utron in both read directions.
    blocks : dic[str: dict[str, set]]
        Nested dictionary of start and ends for each read. Top level is read, with a dictionary of start and end.
    read : str
        Direction of read to filter on, default is ''read1'' but can also use ''read2''.

    Returns
    -------
    dict[str, tuple(Any)]
        Dictionary of the chromosome, start, end and strand of transcripts that are within introns.
    """
    within_intron: dict[str, Any] = {}
    for chromosome, start, end, transcript_id, strand in pair_features[read]["utron"]:
        start_end_within_intron = (
            start <= pair_features[read]["start"] <= end or start <= pair_features[read]["end"] <= end
        )
        spans_intron = (
            pair_features[read]["start"] < start
            and pair_features[read]["end"] > end
            and (end - start) < pair_features[read]["length"]
        )
        if (  # pylint: disable=too-many-boolean-expressions
            (start_end_within_intron or spans_intron)
            # Start should not be in ends and ends should not be in start, can we combine the start and end block
            # sets I wonder?
            and start not in blocks["read1"]["ends"]
            and end not in blocks["read1"]["starts"]
            and start not in blocks["read2"]["ends"]
            and end not in blocks["read2"]["starts"]
        ):
            # Why add an empty list and append a tuple?
            if transcript_id not in within_intron:
                within_intron[transcript_id] = []
            within_intron[transcript_id].append((start, end, chromosome, strand))
    return within_intron


def filter_spliced_utrons(
    pair_features: dict[str, dict[str, Any]],
    blocks: dict[str, dict[str, set[int]]],
    read: str = "read1",
) -> dict[str, list[Any]]:
    """
    Filter utrons where start is in the block ends or end is in the block start.

    Parameters
    ----------
    pair_features : dict[str, dict]
        Dictionary of extracted features and utron in both read directions.
    blocks : dic[str: dict[str, set]]
        Nested dictionary of start and ends for each read. Top level is read, with a dictionary of start and end.
    read : str
        Direction of read to filter on, default is ''read1'' but can also use ''read2''.

    Returns
    -------
    dict[str, tuple(Any)]
        Dictionary of the chromosome, start, end and strand of transcripts that are within introns.
    """
    spliced_3ui: dict[str, list[Any]] = {}
    for chromosome, start, end, transcript_id, strand in pair_features[read]["utron"]:
        if start in blocks[read]["ends"] and end in blocks[read]["starts"]:
            # Why add an empty list and append a tuple?
            if transcript_id not in spliced_3ui:
                spliced_3ui[transcript_id] = []
            spliced_3ui[transcript_id].append((start, end, chromosome, strand))
    return spliced_3ui


def unique_conversions(
    reads1: dict[str, list[Any]],
    reads2: dict[str, list[Any]],
) -> frozenset[list[Any]]:
    """
    Create a unique set of conversions that are to be retained.

    Parameters
    ----------
    reads1 : dict[str, list[tuple[Any]]]
        A dictionary of reads mapped to transcripts (key) which overlap introns. Each read has the ''start'',  ''end'',
       ''chromsome'' and ''strand'' recorded.
    reads2 : dict[str, list[tuple[Any]]]
        A dictionary of reads mapped to transcripts (key) which overlap introns. Each read has the ''start'',  ''end'',
       ''chromsome'' and ''strand'' recorded.

    Returns
    -------
    set[list[Any]]
        Combines the two sets of observations and de-duplicates them, returning only the unique assigned conversions.
    """
    flat1 = [(key, nested_list) for key, values in reads1.items() for nested_list in values]
    flat2 = [(key, nested_list) for key, values in reads2.items() for nested_list in values]
    # logger.debug("Extracted unique conversions in both reads, combining to unique set.")
    return frozenset(flat1 + flat2)  # type: ignore[arg-type]


def remove_common_reads(retained: set[list[Any]], spliced: set[list[Any]]) -> tuple[set[list[Any]], set[list[Any]]]:
    """
    Remove reads that are common to both retained and spliced sets.

    Parameters
    ----------
    retained : set[list[Any]]
        Set of retained reads. Each item is a tuple with ''transcript_id'' and a list of ''start'', ''end'',
        ''chromosome'' and ''strand''.
    spliced : set[list[Any]]
        Set of retained reads. Each item is a tuple with ''transcript_id'' and a list of ''start'', ''end'',
        ''chromosome'' and ''strand''.

    Returns
    -------
    tuple[set[list[Any]], set[list[Any]]]
        A tuple of the ''retained'' (first) and ''spliced'' reads with common items removed.
    """
    common = retained & spliced
    retained -= common
    spliced -= common
    # logger.debug("Removed common elements from retained and spliced.")
    return (retained, spliced)


def conversions_per_read(  # pylint: disable=too-many-positional-arguments
    read: AlignedSegment,
    conversion_from: str,
    conversion_to: str,
    convertible: set[str],
    converted_position: set[str],
    coverage: set[str],
    vcf_file: VariantFile,
) -> tuple[set[str], set[str], set[str]]:
    """
    Build sets of genome position for conversions, converted positions and coverage for a given read.

    Parameters
    ----------
    read : dict[str, dict[str, Any]]
        Aligned read.
    conversion_from : str
        The base pair the conversion is from, typically either ''T'' or ''C''.
    conversion_to : str
        The base pair the conversion is to, typically the opposite pairing of ''from'', i.e. ''A'' or ''C''
        respectively.
    convertible : set
        Set, possibly empty, to which the genome position is added if the sequence at a given location matches
        ''conversion_from''.
    converted_position : set
        Set, possibly, empty, to which the genome position is added if a conversion has occurred.
    coverage : set
        Set, possibly empty, to which the genome position is added for all aligned pairs of a read.
    vcf_file : VariantFile
        VCF file.

    Returns
    -------
    tuple[set[str], set[str], set[str]]
        Three sets of the ''convertible'', ''converted_position'' and ''coverage''.
    """
    # Ensure we have upper case conversions to compare
    conversion_from = conversion_from.upper()
    conversion_to = conversion_to.upper()
    for read_position, genome_position, genome_sequence in read.get_aligned_pairs(with_seq=True):
        if None in (read_position, genome_position, genome_sequence):
            continue
        coverage.add(genome_position)
        if genome_sequence.upper() == conversion_from:
            convertible.add(genome_position)

        # If the sequence at this position has been converted compared to the genome sequence...
        if read.query_sequence[read_position].upper() == conversion_to and genome_sequence.upper() == conversion_from:
            # ...check that this is a new variant at this position? Question : Is this the correctinterpretation?
            variants_at_position = list(vcf_file.fetch(read.reference_name, genome_position, genome_position + 1))
            if variants_at_position:
                if any(variant.alts[0].upper() == conversion_to.upper() for variant in variants_at_position):
                    pass
                else:
                    converted_position.add(genome_position)
            else:
                converted_position.add(genome_position)
    # logger.debug(f"convertible : {convertible}\nconverted_position : {converted_position}\ncoverage : {coverage}")
    return (convertible, converted_position, coverage)


def count_conversions_across_pairs(
    forward_read: dict[str, dict[str, Any]],
    reverse_read: dict[str, dict[str, Any]],
    vcf_file: VariantFile,
    forward_conversion: dict[str, str] | None = None,
    reverse_conversion: dict[str, str] | None = None,
) -> dict[str, int]:
    """
    Count conversions across paired reads.

    Parameters
    ----------
    forward_read : dict[str, dict[str, Any]]
        Aligned segment for forward read.
    reverse_read : dict[str, dict[str, Any]]
        Aligned segment for reversed read.
    vcf_file : VariantFile
        Variant File.
    forward_conversion : dict, optional
        Forward conversion dictionary typically ''{"from": "A", "to": "G"}''.
    reverse_conversion : dict, optional
        Reverse conversion, typically ''{"from": "T", "to": "C"}''.

    Returns
    -------
    tuple[int, int, int]
        Tuple of the number of convertible base pairs, the number of conversions and the coverage of the paired
      alignments.

    Raises
    ------
    ValueError
        ValueError is raised if either ''forward_conversion'' or ''reverse_conversion'' is ''None''.
    """
    if forward_conversion is None:
        raise ValueError("forward_conversion can not be empty.")
    if reverse_conversion is None:
        raise ValueError("reverse_conversion can not be empty.")

    # Count conversions on the forward read
    convertible, converted_position, coverage = conversions_per_read(
        forward_read,
        forward_conversion["from"],
        forward_conversion["to"],
        convertible=set(),
        converted_position=set(),
        coverage=set(),
        vcf_file=vcf_file,
    )
    # Count conversions on the reverse read
    convertible, converted_position, coverage = conversions_per_read(
        reverse_read,
        reverse_conversion["from"],
        reverse_conversion["to"],
        convertible,
        converted_position,
        coverage,
        vcf_file,
    )
    # logger.debug("Counted conversions paired reads")
    return {"convertible": len(convertible), "converted_position": len(converted_position), "coverage": len(coverage)}


def append_data(  # pylint: disable=too-many-positional-arguments
    assigned_conversions: set[list[Any]],
    coverage_counts: dict[str, int],
    read_uid: int,
    assignment: str,
    results: pl.DataFrame,
    schema: dict[str, type],
) -> pl.DataFrame:
    """
    Create a Polars dataframe combining the ''assigned_conversions'' and ''coverage_counts''.

    Adds ''assignment'' to the resulting dataframe.

    Parameters
    ----------
    assigned_conversions : set[list[Any]]
        A set of assigned conversions. Each element of the set is a list of key features (CHECK WHAT THESE ARE).
    coverage_counts : dict[str, int] dest_dir: str | Path
        A dictionary of coverage counts indexed by CHECK.
    read_uid : int
        Integer representing the unique read ID.
    assignment : str
        Type of assignment, either ''Rep'' or ''Spl'' (for Splice).
    results : pl.DataFrame
        Polars DataFrame to append data to. This will initially be empty but the schema matches the variables that are
        added.
    schema : dict[str, type]
        Schema dictionary for data frame.

    Returns
    -------
    pl.DataFrame
        Returns a Polars DataFrame of the data structure.
    """
    if results is None:
        results = pl.DataFrame(schema=schema)
    for transcript_id, position in assigned_conversions:
        start, end, chromosome, strand = position
        row = pl.DataFrame(
            data={
                "read_uid": read_uid,
                "transcript_id": transcript_id,
                "start": start,
                "end": end,
                "chr": chromosome,
                "strand": strand,
                "assignment": assignment,
                "conversions": coverage_counts["converted_position"],
                "convertible": coverage_counts["convertible"],
                "coverage": coverage_counts["coverage"],
            },
            schema=schema,
        )
        results = pl.concat([results, row])
    return results.sort(by=["read_uid", "transcript_id", "chr", "start", "end"])
