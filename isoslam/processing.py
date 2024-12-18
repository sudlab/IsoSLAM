"""Entry point, sub-parsers and arguments and processing functions."""

import argparse as arg
import sys
from pathlib import Path
from typing import Any

from loguru import logger

from isoslam import __version__, io, logging, summary


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
        "--file-pattern",
        dest="file_pattern",
        type=str,
        required=False,
        default="*_summarized.tsv",
        help="Regular expression for summarized files to process.",
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


def process(args: arg.Namespace | None) -> None:  # pylint: disable=unused-argument
    """
    Process a set of files.

    Parameters
    ----------
    args : arg.Namespace | None
        Arguments function was invoked with.

    Returns
    -------
    None
        Function does not return anything.
    """
    # config = io.read_yaml() if args.config is None else io.read_yaml(args.config) # type: ignore[call-arg,union-attr]
    return


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
    summary_counts.sort_values(by=["Chr", "Transcript_id", "Start"], inplace=True)
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
