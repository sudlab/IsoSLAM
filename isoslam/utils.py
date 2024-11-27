"""Utilities and helper functionsfor IsoSLAM."""

from argparse import Namespace

from loguru import logger

from isoslam import io


def update_config(config: dict, args: dict | Namespace) -> dict:  # type: ignore[type-arg]
    """
    Update the configuration with any arguments.

    Parameters
    ----------
    config : dict
        Dictionary of configuration (typically read from YAML file specified with '-c/--config <filename>').
    args : Namespace
        Command line arguments.

    Returns
    -------
    dict
        Dictionary updated with command arguments.
    """
    args = vars(args) if isinstance(args, Namespace) else args

    config_keys = config.keys()
    for arg_key, arg_value in args.items():
        if isinstance(arg_value, dict):
            update_config(config, arg_value)
        else:
            if arg_key in config_keys and arg_value is not None:
                original_value = config[arg_key]
                config[arg_key] = arg_value
                logger.debug(f"Updated config config[{arg_key}] : {original_value} > {arg_value} ")
    if "base_dir" in config.keys():
        config["base_dir"] = io._str_to_path(config["base_dir"])  # pylint: disable=protected-access
    if "output_dir" in config.keys():
        config["output_dir"] = io._str_to_path(config["output_dir"])  # pylint: disable=protected-access
    return config
