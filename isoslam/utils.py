"""Utilities and helper functionsfor IsoSLAM."""

from argparse import Namespace
from typing import Any

from loguru import logger

from isoslam import io


def update_config(config: dict[str, Any], args: dict[str, Any]) -> dict[str, Any]:
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
    args_keys = args.keys()
    for config_key, config_value in config.items():
        if isinstance(config_value, dict):
            update_config(config_value, args)
        else:
            if config_key in args_keys and args[config_key] is not None and config_value is not args[config_key]:
                original_value = config[config_key]
                config[config_key] = args[config_key]
                logger.info(f"Updated config config[{config_key}] : {original_value} > {args[config_key]} ")
    if "base_dir" in config.keys():
        config["base_dir"] = io._str_to_path(config["base_dir"])  # pylint: disable=protected-access
    if "output_dir" in config.keys():
        config["output_dir"] = io._str_to_path(config["output_dir"])  # pylint: disable=protected-access
    return config
