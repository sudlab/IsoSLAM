"""Entry point, sub-parsers and arguments and processing functions."""

import argparse as arg
import sys
from pathlib import Path
from typing import Any

import polars as pl
from loguru import logger

from isoslam import __version__, io, isoslam, logging, summary, validation

# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-statements


def create_parser() -> arg.ArgumentParser:
    """
    Create a parser for reading options.

    Parser is created with multiple sub-parsers for eading options to run ``isoslam``.

    Returns
    -------
    arg.ArgumentParser
        Argument parser.
    """
    parser = arg.ArgumentParser(
        description="Run various programs related to IsoSLAM. Add the name of the program you wish to run."
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Installed version of IsoSlam : {__version__}",
        help="Report the installed version of IsoSLAM.",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        dest="config_file",
        type=Path,
        required=False,
        help="Path to a YAML configuration file.",
    )
    parser.add_argument(
        "-b",
        "--base-dir",
        dest="base_dir",
        type=Path,
        required=False,
        help="Base directory to run isoslam on.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=Path,
        required=False,
        help="Output directory to write results to.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        dest="log_level",
        type=str,
        required=False,
        help="Logging level to use, default is 'info' for verbose output use 'debug'.",
    )

    subparsers = parser.add_subparsers(title="program", description="Available programs listed below:", dest="program")

    # Create sub-parsers for different stages

    # Process processes all files
    process_parser = subparsers.add_parser(
        "process",
        description="Process all files and run all summary plotting and statistics.",
        help="Process all files and run all summary plotting and statistics.",
    )
    process_parser.add_argument(
        "-b",
        "--bam-file",
        dest="bam_file",
        type=Path,
        required=False,
        help="Path to '.bam' file that has undergone read assignment with 'featureCount'.",
    )
    process_parser.add_argument(
        "-g", "--gtf-file", dest="gtf_file", type=Path, required=False, help="Path to '.gtf' transcript assembly file."
    )
    process_parser.add_argument(
        "-d",
        "--bed-file",
        dest="bed_file",
        type=Path,
        required=False,
        help="Path to '.bed' utron file. Must be bed6 format.",
    )
    process_parser.add_argument(
        "-v", "--vcf-file", dest="vcf_file", type=Path, required=False, help="Path to '.vcf.gz' file."
    )
    process_parser.add_argument(
        "-u",
        "--upper-pairs-limit",
        dest="upper_pairs_limit",
        type=int,
        required=False,
        help="Upper limit of pairs to be processed.",
    )
    process_parser.add_argument(
        "-f",
        "--first-matched-limit",
        dest="first_matched_limit",
        type=int,
        required=False,
        help="Limit of matches.",
    )
    process_parser.add_argument(
        "--delim",
        dest="delim",
        type=str,
        required=False,
        help="Delimiter to use in output.",
    )
    process_parser.add_argument(
        "--output-file",
        dest="output_file",
        type=str,
        required=False,
        help="File to write results to.",
    )
    process_parser.set_defaults(func=process)

    # Create configuration sub-parser
    create_config_parser = subparsers.add_parser(
        "create-config",
        description="Create a configuration file using the defaults.",
        help="Create a configuration file using the defaults.",
    )
    create_config_parser.add_argument(
        "-f",
        "--filename",
        dest="filename",
        type=Path,
        required=False,
        default="config.yaml",
        help="Name of YAML file to save configuration to (default 'config.yaml').",
    )
    create_config_parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=Path,
        required=False,
        default="./",
        help="Path to where the YAML file should be saved (default './' the current directory).",
    )
    create_config_parser.set_defaults(func=io.create_config)

    # Summarise counts sub-parser
    summary_counts_parser = subparsers.add_parser(
        "summary-counts",
        description="Summarise the counts.",
        help="Summarise the counts.",
    )
    summary_counts_parser.add_argument(
        "--file-ext",
        dest="file_ext",
        type=str,
        required=False,
        default=".tsv",
        help="File extension of summarized files to process.",
    )
    summary_counts_parser.add_argument(
        "--directory",
        dest="directory",
        type=Path,
        required=False,
        help="Directory to search for input files.",
    )
    summary_counts_parser.add_argument(
        "--conversions-var",
        dest="conversions_var",
        type=str,
        required=False,
        help="Name of column that holds details of conversions.",
    )
    summary_counts_parser.add_argument(
        "--conversions-threshold",
        dest="conversions_threshold",
        type=int,
        required=False,
        help="Minimum number of conversions.",
    )
    summary_counts_parser.add_argument(
        "--test-file",
        dest="test_file",
        type=str,
        required=False,
        help="Pattern used in test file names.",
    )
    summary_counts_parser.add_argument(
        "--filename-var",
        dest="filename_var",
        type=str,
        required=False,
        help="Name of column that holds file names.",
    )
    summary_counts_parser.add_argument(
        "--regex",
        dest="regex",
        type=str,
        required=False,
        help="Regular expression for extracting day/hour/replication from filenames.",
    )
    summary_counts_parser.add_argument(
        "--outfile",
        dest="outfile",
        type=Path,
        required=False,
        default="summary_counts.tsv",
        help="Output filename to save results to, will be nested under 'output_dir'.",
    )
    summary_counts_parser.add_argument(
        "--separator",
        dest="sep",
        type=str,
        required=False,
        default="\t",
        help="Field separator to use in output file, default is '\t' but other values (e.g. ',' are allowed).",
    )
    summary_counts_parser.set_defaults(func=summarise_counts)

    # Additional parsers for future functionality
    # summarize_counts_parser = subparsers.add_parser(
    #     "summarize",
    #     description="Summarize counts.",
    #     help="Summarize counts.",
    # )
    # summarize_counts_parser.set_defaults(func=summarize_counts)
    # plot_conversions_parser = subparsers.add_parser(
    #     "plot_conversions",
    #     description="Plot conversions.",
    #     help="Plot conversions.",
    # )
    # plot_conversions_parser.set_defaults(func=plot_conversions)
    return parser


def process(
    args: arg.Namespace | None,
) -> pl.DataFrame:
    """
    Process a set of files.

    Parameters
    ----------
    args : arg.Namespace | None
        Arguments function was invoked with.

    Returns
    -------
    pl.DataFrame
        Polars Dataframe of results.
    """
    config = io.load_and_update_config(args)
    logger.remove()
    if "log_level" in vars(args) and vars(args)["log_level"] is not None:
        logging.setup(level=vars(args)["log_level"])
    else:
        logging.setup(level=config["log_level"])
    validation.validate_config(config=config, schema=validation.DEFAULT_CONFIG_SCHEMA, config_type="configuration")

    # Load files...
    vcffile = io.load_file(config["vcf_file"])
    utron_coords = isoslam.extract_transcripts(config["bed_file"])
    strand_dict, tx2gene = isoslam.extract_strand_transcript(config["gtf_file"])
    # Setup Polars dataframe
    results = pl.DataFrame(schema=config["schema"])
    pairs_processed = 0
    first_matched = 0
    read_uid = 1
    # Process the BAM file
    for pair in isoslam.extract_segment_pairs(config["bam_file"]):
        # If there are an excessive number of pairs and/or matches per file we break out
        if pairs_processed >= config["upper_pairs_limit"]:
            break
        if first_matched > config["first_matched_limit"]:
            break
        # Perform a number of checks as to whether we should proceed with processing
        # 1. If we don't have a pair skip
        #    @ns-rse : I tried unpacking and found there were indeed instances where the "pairs" could got upto 8 in
        #    length?
        if len(pair) != 2:
            continue
        # 2. If either reads of the pair are unmapped skip.
        read1, read2 = pair
        if read1.is_unmapped or read2.is_unmapped:
            continue
        # Use the intersection of sets, this skips if either read1 or read2 aren't assigned, - or +
        pair_features = isoslam.extract_features_from_pair(pair)
        if not {"Assigned", "+", "-"} & {pair_features["read1"]["status"], pair_features["read2"]["status"]}:
            continue
        # Check that only Assigned, + and - status are included
        status_list = [pair_features["read1"]["status"], pair_features["read2"]["status"]]
        if all(status in ["Assigned", "+", "-"] for status in status_list):
            # If strands are equal assign one to strand
            if strand_dict[pair_features["read1"]["transcript"]] == strand_dict[pair_features["read2"]["transcript"]]:
                strand = strand_dict[pair_features["read1"]["transcript"]]
            # ...if not we skip this pair
            else:
                continue
        # If pairs are not both Assigned/+/- we set strand to the transcript from read1, otherwise its read2
        else:
            strand = (
                strand_dict[pair_features["read1"]["transcript"]]
                if pair_features["read1"]["status"] == "Assigned"
                else strand_dict[pair_features["read2"]["transcript"]]
            )
        # Set forward and reverse reads, if they don't match then we skip
        if read1.is_reverse and not read2.is_reverse:
            reverse_read = read1
            forward_read = read2
        elif read2.is_reverse and not read1.is_reverse:
            reverse_read = read2
            forward_read = read1
        else:
            # Not proper pair
            continue

        # We _haven't_ skipped any pairs so increment the pair counter
        pairs_processed += 1
        # Processing now begins...
        # Extract utron for the gene
        pair_features["read1"]["utron"] = isoslam.extract_utron(
            features=pair_features["read1"], gene_transcript=tx2gene, coordinates=utron_coords
        )
        pair_features["read2"]["utron"] = isoslam.extract_utron(
            features=pair_features["read2"], gene_transcript=tx2gene, coordinates=utron_coords
        )
        # Get blocks
        block_starts1, block_ends1 = isoslam.zip_blocks(read1)
        block_starts2, block_ends2 = isoslam.zip_blocks(read2)
        blocks = {
            "read1": {"starts": block_starts1, "ends": block_ends1},
            "read2": {"starts": block_starts2, "ends": block_ends2},
        }
        # Retain within introns
        read1_within_intron = isoslam.filter_within_introns(pair_features, blocks, read="read1")
        read2_within_intron = isoslam.filter_within_introns(pair_features, blocks, read="read2")
        # Retain spliced
        read1_spliced_3UI = isoslam.filter_spliced_utrons(pair_features, blocks, read="read1")
        read2_spliced_3UI = isoslam.filter_spliced_utrons(pair_features, blocks, read="read2")
        # List of all dictionaries
        retained_reads = [read1_within_intron, read2_within_intron, read1_spliced_3UI, read2_spliced_3UI]
        # Check that there are some retained regions (all dictionaries would be empty if there are none retained and
        # not {} evaluates to True, hence wrapping in all())
        if all(not contents for contents in retained_reads):
            continue
        # We have got a match so increment counts
        first_matched += 1

        # Unique conversions within introns to be retained
        assign_conversions_to_retained = isoslam.unique_conversions(read1_within_intron, read2_within_intron)
        # Unique conversions within 3UI to be retained
        assign_conversions_to_spliced = isoslam.unique_conversions(read1_spliced_3UI, read2_spliced_3UI)
        ## If there are any events in both we want to remove them - this should be rare
        assign_conversions_to_retained, assign_conversions_to_spliced = isoslam.remove_common_reads(
            assign_conversions_to_retained, assign_conversions_to_spliced
        )
        # If we are mapped to a +ve stranded transcript, then count T>C in the forward read and A>G in the reverse read.
        if strand == "+":
            coverage_counts = isoslam.count_conversions_across_pairs(
                forward_read=forward_read,
                reverse_read=reverse_read,
                vcf_file=vcffile,
                forward_conversion=config["forward_reads"],
                reverse_conversion=config["reverse_reads"],
            )
        elif strand == "-":
            # If we are mapped to a -ve stranded transcript, count T>C in the reverse read and A>G in the forward read.
            # NB - can either flip the reads or the conversion that are passed in, here we flip the read
            coverage_counts = isoslam.count_conversions_across_pairs(
                forward_read=reverse_read,
                reverse_read=forward_read,
                vcf_file=vcffile,
                forward_conversion=config["forward_reads"],
                reverse_conversion=config["reverse_reads"],
            )
        else:
            # should not be possible - but just in case
            pass

        results = isoslam.append_data(
            assigned_conversions=assign_conversions_to_retained,
            coverage_counts=coverage_counts,  # pylint: disable=possibly-used-before-assignment
            read_uid=read_uid,
            assignment="Ret",
            results=results,
            schema=config["schema"],
        )
        results = isoslam.append_data(
            assigned_conversions=assign_conversions_to_spliced,
            coverage_counts=coverage_counts,  # pylint: disable=possibly-used-before-assignment
            read_uid=read_uid,
            assignment="Spl",
            results=results,
            schema=config["schema"],
        )
        read_uid += 1

    results = results.sort(by=["read_uid", "transcript_id", "chr", "start", "end"])
    io.data_frame_to_file(data=results, output_dir=config["output_dir"], outfile=config["output_file"])

    return results


def summarise_counts(args: arg.Namespace | None) -> None:
    """
    Take a set of output files and summarise the number of conversions.

    Counts are made within file, chromosome, transcript, start, end, assignment and whether there is one or more
    conversion observed.

    Parameters
    ----------
    args : arg.Namespace | None
        Arguments function was invoked with.

    Returns
    -------
    None
        Function does not return anything.
    """
    # Load the configuration file (default or user) and update with supplied flags
    config = io.load_and_update_config(args)
    logger.remove()
    if vars(args)["log_level"] is not None:
        logging.setup(level=vars(args)["log_level"])
    else:
        logging.setup(level=config["log_level"])
    summary_counts_config = config["summary_counts"]
    output_config = summary_counts_config.pop("output")
    output_config["output_dir"] = config["output_dir"]
    summary_counts = summary.summary_counts(**summary_counts_config)
    summary_counts = summary_counts.sort(by=["Chr", "Transcript_id", "Start"])
    io.data_frame_to_file(summary_counts, **output_config)
    logger.info(f"Summary counts file written to : {output_config['output_dir']}/{output_config['outfile']}")


def entry_point(manually_provided_args: list[Any] | None = None, testing: bool = False) -> None | arg.Namespace:
    """
    Entry point for all IsoSLAM programs.

    Main entry point for running 'isoslam' which allows the different processing, plotting and testing modules to be
    run.

    Parameters
    ----------
    manually_provided_args : None
        Manually provided arguments.
    testing : bool
        Whether testing is being carried out.

    Returns
    -------
    None
        Function does not return anything.
    """
    parser = create_parser()
    args = parser.parse_args() if manually_provided_args is None else parser.parse_args(manually_provided_args)

    # If no module has been specified print help and exit
    if not args.program:
        parser.print_help()
        sys.exit()

    if testing:
        return args

    # Run the specified module(s)
    args.func(args)

    return None
