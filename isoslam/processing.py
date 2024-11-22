"""Functions for Processing."""

import argparse as arg
import sys
from pathlib import Path
from typing import Any

# from isoslam import __version__


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
        "--version",
        # version=f"Installed version of IsoSlam : {__version__}",
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


def process(args: arg.Namespace | None) -> None:
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
    return


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
    print(f"{args.program=}")
    if not args.program:
        parser.print_help()
        sys.exit()

    if testing:
        return args

    # Run the specified module(s)
    args.func(args)

    return None
