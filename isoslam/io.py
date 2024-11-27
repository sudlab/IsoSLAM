"""Module for reading files."""

import argparse
from datetime import datetime
from importlib import resources
from pathlib import Path

# from cgatcore import iotools
from loguru import logger

# import pysam
from ruamel.yaml import YAML, YAMLError

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


def _path_to_str(config: dict) -> dict:  # type: ignore[type-arg]
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


def read_yaml(filename: str | Path) -> dict | None:  # type: ignore[type-arg]
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
    with Path(filename).open(encoding="utf-8") as f:
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


def load_bam() -> None:
    """Load '.bam' file."""
    return


def load_bed() -> None:
    """Load '.bed' file."""
    return


def load_gtf() -> None:
    """Load '.bed' file."""
    return


def load_vcf() -> None:
    """Load '.vcf' file."""
    return
