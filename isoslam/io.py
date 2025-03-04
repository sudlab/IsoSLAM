"""Module for reading and writing files."""

import argparse
import gzip
import re
from collections.abc import Callable, Generator
from datetime import datetime
from importlib import resources
from io import TextIOWrapper
from pathlib import Path
from typing import Any, TextIO

import pandas as pd
import polars as pl
import pysam
from loguru import logger
from ruamel.yaml import YAML, YAMLError

from isoslam import utils

CONFIG_DOCUMENTATION_REFERENCE = """# For more information on configuration and how to use it see:
# https://sudlab.github.io/IsoSLAM\n"""


def _get_date_time() -> str:
    """
    Get a date and time for adding to generated files or logging.

    Returns
    -------
    str
        A string of the current date and time, formatted appropriately.
    """
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _str_to_path(path: str | Path) -> Path:
    """
    Ensure path is a Path object.

    Returns the current directory of passed './'.

    Parameters
    ----------
    path : str | Path
        Path to be converted.

    Returns
    -------
    Path
        Pathlib object of supplied path.
    """
    return Path().cwd() if path == "./" else Path(path).expanduser()


def _path_to_str(config: dict[str, Any]) -> dict[str, Any]:
    """
    Recursively traverse a dictionary and convert any Path() objects to strings for writing to YAML.

    Parameters
    ----------
    config : dict
        Dictionary to be converted.

    Returns
    -------
    Dict:
        The same dictionary with any Path() objects converted to string.
    """
    for key, value in config.items():
        if isinstance(value, dict):
            _path_to_str(value)
        elif isinstance(value, Path):
            config[key] = str(value)
    return config


def read_yaml(filename: str | Path | None = None) -> dict[str, Any] | None:
    """
    Read a YAML file.

    Parameters
    ----------
    filename : Union[str, Path]
        YAML file to read.

    Returns
    -------
    Dict
        Dictionary of the file.
    """
    if filename is None:
        filename = resources.files(__package__) / "default_config.yaml"  # type: ignore[assignment]
    with Path(filename).open(encoding="utf-8") as f:  # type: ignore[arg-type]
        try:
            yaml_file = YAML(typ="safe")
            return yaml_file.load(f)  # type: ignore[no-any-return]
        except YAMLError as exception:
            logger.error(exception)
            return {}


def write_yaml(
    config: dict,  # type: ignore[type-arg]
    output_dir: str | Path,
    config_file: str = "config.yaml",
    header_message: str | None = None,
) -> None:
    """
    Write a configuration (stored as a dictionary) to a YAML file.

    Parameters
    ----------
    config : dict
        Configuration dictionary.
    output_dir : Union[str, Path]
        Path to save the dictionary to as a YAML file (it will be called 'config.yaml').
    config_file : str
        Filename to write to.
    header_message : str
        String to write to the header message of the YAML file.
    """
    output_config = Path(output_dir) / config_file
    # Revert PosixPath items to string
    config = _path_to_str(config)

    if header_message:
        header = f"# {header_message} : {_get_date_time()}\n" + CONFIG_DOCUMENTATION_REFERENCE
    else:
        header = f"# Configuration from IsoSLAM run completed : {_get_date_time()}\n" + CONFIG_DOCUMENTATION_REFERENCE
    output_config.write_text(header, encoding="utf-8")

    yaml = YAML(typ="safe")
    with output_config.open("a", encoding="utf-8") as f:
        try:
            yaml.dump(config, f)
        except YAMLError as exception:
            logger.error(exception)


def load_and_update_config(args: argparse.Namespace | None) -> dict[str, Any]:
    """
    Load a configuration file to dictionary and update entries with user supplied arguments.

    If ''args'' does not contain any value for ''args.config_file'' the default configuration
    (''isoslam/default_config.yaml'') is loaded, otherwise the user specified configuration is loaded.

    Once the configuration is loaded any user specified options update the dictionary.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments supplied by user.

    Returns
    -------
    dict[str: Any]
        Dictionary of configuration optionsupdated with user specified options.
    """
    config = read_yaml() if vars(args)["config_file"] is None else read_yaml(vars(args)["config_file"])
    config["schema"] = _type_schema(config["schema"])  # type: ignore[index]
    return utils.update_config(config, vars(args))  # type: ignore[arg-type]


def _type_schema(schema: dict[str, str]) -> dict[str, type]:
    """
    Convert schemas to types.

    When read from a YAML file the schema's values are strings, they need converting to types to work with Polars.

    Parameters
    ----------
    schema : dict[str, str]
        Dictionary of schema types with values as strings.

    Returns
    -------
    dict[str, Type]
        Returns dictionary with schema types as such.
    """
    return {key: eval(value) for key, value in schema.items()}  # pylint: disable=eval-used  # noqa: S307


def create_config(args: argparse.Namespace | None = None) -> None:
    """
    Write the default configuration file to disk.

    Parameters
    ----------
    args : argparse.Namespace | None
        Optional arguments to parse.
    """
    filename = "config" if args.filename is None else args.filename  # type: ignore [union-attr]
    output_dir = Path("./") if args.output_dir is None else Path(args.output_dir)  # type: ignore [union-attr]
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = resources.files(__package__) / "default_config.yaml"
    config = config_path.read_text()

    if ".yaml" not in str(filename) and ".yml" not in str(filename):
        create_config_path = output_dir / f"{filename}.yaml"
    else:
        create_config_path = output_dir / filename

    with create_config_path.open("w", encoding="utf-8") as f:
        f.write(f"# Config file generated {_get_date_time()}\n")
        f.write(f"{CONFIG_DOCUMENTATION_REFERENCE}")
        f.write(config)
    logger.info(f"A sample configuration file has been written to : {str(create_config_path)}")
    logger.info(CONFIG_DOCUMENTATION_REFERENCE)


def load_file(file_path: str | Path) -> Any:
    """
    Load files of different types.

    Supports the following file types...

    * ``.bam`` - The sequence data that is to be analysed.
    * ``.bed`` - The locations of introns/splice junctions.
    * ``.gtf`` - Transcript structures from which the ``.bed`` file is derived.
    * ``.vcf`` - Locations of known sequences difference from the reference sequence.

    Parameters
    ----------
    file_path : str | Path
        Path to file to load.

    Returns
    -------
    Any
        Returns the loaded file as an object.
    """
    file_suffix = Path(file_path).suffix
    if file_suffix == ".gz":
        file_suffix = "".join(Path(file_path).suffixes)
    loader = _get_loader(file_suffix)
    return loader(file_path)


def _get_loader(file_ext: str = "bam") -> Callable:  # type: ignore[type-arg]
    """
    Creator component which determines which file loader to use.

    Parameters
    ----------
    file_ext : str
        File extension of file to be loaded.

    Returns
    -------
    function
        Returns the function appropriate for the required file type to be loaded.

    Raises
    ------
    ValueError
        Unsupported file extension results in ValueError.
    """
    if file_ext == ".bam":
        return _load_bam
    if file_ext in (".bed", ".bed.gz"):
        return _load_bed
    if file_ext == ".gtf":
        return _load_gtf
    if file_ext in (".vcf", ".vcf.gz"):
        return _load_vcf
    raise ValueError(file_ext)


def _load_bam(bam_file: str | Path) -> pysam.libcalignmentfile.AlignmentFile:
    """
    Load ``.bam`` file.

    ``.bam`` files are the sequence data that is to be analysed.

    Parameters
    ----------
    bam_file : str | Path
        Path, as string or pathlib Path, to a '.bam' file that is to be loaded.

    Returns
    -------
    pysam.libcalignmentfile.AlignmentFile
        Loads the specified alignment file.
    """
    try:
        return pysam.AlignmentFile(bam_file)
    except FileNotFoundError as e:
        raise e


def _load_bed(bed_file: str | Path) -> TextIO:
    """
    Open ``.bed`` file for reading, supports gzip compressed formats.

    ``.bed`` files contain the locations of introns/splice junctions.

    Parameters
    ----------
    bed_file : str | Path
        Path, as string or pathlib Path, to a '.bed' or '.bed.gz' file that is to be loaded.

    Returns
    -------
    TextIO
        Returns a connection to an open file object.
    """
    try:
        if Path(bed_file).suffix == ".gz":
            return gzip.open(bed_file, "rt", encoding="utf-8")
        return Path(bed_file).open(mode="r", encoding="utf-8")
    except OSError as e:
        raise e


def _load_gtf(gtf_file: str | Path) -> pysam.libctabix.tabix_generic_iterator:
    """
    Load ``.gtf`` file and return as an iterable.

    ``.gtf`` files contain the transcript structures from which the ``.bed`` file is derived.

    Parameters
    ----------
    gtf_file : str | Path
        Path, as string or pathlib Path, to a '.gtf' file that is to be loaded.

    Returns
    -------
    pysam.libctabix.tabix_generic_iterator
        Iterator of GTF file.
    """
    try:
        return pysam.tabix_iterator(Path(gtf_file).open(encoding="utf8"), parser=pysam.asGTF())
    except FileNotFoundError as e:
        raise e


def _load_vcf(vcf_file: str | Path) -> pysam.libcbcf.VariantFile:
    """
    Load ``.vcf`` file.

    ``.vcf`` files contain the locations of known sequences difference from the reference sequence. Any ``T > C``
    (e.g. SNPs) conversions that match to the location of known seuqence variations will be removed. These can be
    obtained from a reference collection of variation data (such as `dbSNP
    <https://www.ncbi.nlm.nih.gov/projects/SNP/get_html.cgi?whichHtml=overview>`_) or derived directly from the RNAseq
    reads.

    Parameters
    ----------
    vcf_file : str | Path
        Path, as string or pathlib Path, to a '.vcf' file that is to be loaded.

    Returns
    -------
    pysam.libcbcf.VariantFile
        Loads the specified VCF file.
    """
    try:
        return pysam.VariantFile(vcf_file)
    except FileNotFoundError as e:
        raise e


def _find_files(pattern: str = "**/*.tsv", directory: str | Path | None = None) -> Generator:  # type: ignore[type-arg]
    """
    Find files that match the given pattern.

    Parameters
    ----------
    pattern : str
        Pattern (regular expression) of files to search for.
    directory : str | Path | None
        Directory to search for files, if ''None'' the current working directory is used.

    Returns
    -------
    Generator[_P, None, None]
        A generator of files found that match the given pattern.
    """
    if directory is None:
        directory = Path.cwd()
    return Path(directory).rglob(pattern)


def load_output_files(
    file_ext: str = ".tsv", directory: str | Path | None = None, columns: list[str] | None = None
) -> dict[str, pl.DataFrame]:
    """
    Read a set of files into a list of Polars DataFrames.

    Supports reading ''.parquet'', ''.tsv'' and ``.csv``.

    Parameters
    ----------
    file_ext : str
        File name pattern to search for.
    directory : str | Path | None
        Directory to search for files.
    columns : list[str]
        List of column names to load, defaults will be set if ''None''.

    Returns
    -------
    list[pl.DataFrame]
        A list of Polars DataFrames of each file found.
    """
    if columns is None:
        columns = [
            "Read_UID",
            "Transcript_id",
            "Start",
            "End",
            "Chr",
            "Strand",
            "Assignment",
            "Conversions",
            "Convertible",
            "Coverage",
        ]
    pattern = f"*{file_ext}"
    if file_ext[file_ext.rfind(".") :] == ".parquet":
        results = {_file.stem: pl.read_parquet(_file, columns=columns) for _file in _find_files(pattern, directory)}
    else:
        if file_ext == ".tsv":
            separator = "\t"
        if file_ext == ".csv":
            separator = ","
        results = {
            _file.stem: pl.read_csv(_file, columns=columns, separator=separator)
            for _file in _find_files(pattern, directory)
        }
    return {key: df.with_columns(filename=pl.lit(key)) for key, df in results.items()}


def data_frame_to_file(
    data: pd.DataFrame | pl.DataFrame,
    output_dir: str | Path = "./output/",
    outfile: str = "summary_counts.tsv",
    sep: str = "\t",
    **kwargs: dict[Any, Any],
) -> None:
    """
    Write a Pandas DataFrame to disk.

    Parameters
    ----------
    data : pd.DataFrame | pl.DataFrame
        Pandas DataFrame to write to disk.
    output_dir : str | Path
        Location to write the output to, default is ''./output''.capitalize.
    outfile : str
        Filename to write data to.
    sep : str
        Separator to use in output file.
    **kwargs
        Dictionary of keyword arguments to pass to ''pandas.DataFrame.to_csv()''.
    """
    outdir_file = Path(output_dir) / f"{outfile}"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if isinstance(data, pl.DataFrame):
        try:
            if re.search(r"parquet$", str(outfile)):
                data.write_parquet(outdir_file, **kwargs)
            elif re.search(r"\..sv$", str(outfile)):
                data.write_csv(outdir_file, separator=sep, **kwargs)
            logger.debug(f"File written to : {outdir_file}")
        except Exception as e:
            raise e
    elif isinstance(data, pd.DataFrame):
        try:
            if re.search(r"parquet", str(outfile)):
                data.to_parquet(outdir_file, **kwargs)
            elif re.search(r"\..sv$", str(outfile)):
                data.to_csv(outdir_file, sep=sep, **kwargs)
            logger.debug(f"File written to : {outdir_file}")
        except Exception as e:
            raise e
    else:
        raise TypeError(f"Can not write output Pandas or Polar Dataframe object not supplied = {type(data)=}")


def write_assigned_conversions(  # pylint: disable=too-many-positional-arguments
    assigned_conversions: set[list[Any]],
    coverage_counts: dict[str, int],
    read_uid: int,
    assignment: str,
    outfile: TextIOWrapper,
    delim: str,
) -> None:
    r"""
    Write assigned conversions to files.

    Combines the ''coverage_counts'' with the ''assigned_conversions'' and outputs to disk at the specified location and
    filename with configurable delimiter.

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
    outfile : Any
        Open connection to write results to.
    delim : str
        Delimiter to be used between fields, typically '','' for ''.csv'' or ''\t'' for ''.tsv'' output.
    """
    for transcript_id, position in assigned_conversions:
        start, end, chromosome, strand = position
        outfile.write(
            f"{read_uid}{delim}{transcript_id}{delim}"
            f"{start}{delim}{end}{delim}{chromosome}{delim}"
            f"{strand}{delim}{assignment}{delim}{coverage_counts['converted_position']}{delim}"
            f"{coverage_counts['convertible']}{delim}{coverage_counts['coverage']}\n"
        )
