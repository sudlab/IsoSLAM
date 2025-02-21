"""Validation of configuration."""

from pathlib import Path
from typing import Any

from loguru import logger
from schema import Or, Schema, SchemaError

# pylint: disable=eval-used


def validate_config(config: dict[str, Any], schema: Schema, config_type: str) -> None:
    """
    Validate configuration.

    Parameters
    ----------
    config : dict
        Config dictionary imported by read_yaml() and parsed through clean_config().
    schema : Schema
        A schema against which the configuration is to be compared.
    config_type : str
        Description of configuration being validated.
    """
    try:
        schema.validate(config)
        logger.info(f"The {config_type} is valid.")
    except SchemaError as schema_error:
        raise SchemaError(
            f"There is an error in your {config_type} configuration. "
            "Please refer to the first error message above for details"
        ) from schema_error


DEFAULT_CONFIG_SCHEMA = Schema(
    {
        "log_level": Or(
            "debug",
            "info",
            "warning",
            "error",
            error="Invalid value in config for 'log_level', valid values are"
            "'info' (default), 'debug', 'error' or 'warning",
        ),
        "output_dir": Path,
        "output_file": Or(str, error="Invalid value in config for output_file, should be 'str'"),
        "bam_file": Path,
        "gtf_file": Path,
        "bed_file": Path,
        "vcf_file": Path,
        "upper_pairs_limit": Or(int, error="Invalid value in config for upper_pairs_limit should be 'int'"),
        "first_matched_limit": Or(int, error="Invalid value in config for first_matched_limit should be 'int'"),
        "forward_reads": {
            "from": Or(
                "A",
                "C",
                "T",
                "G",
                error="Invalid value in config for forward_reads.from, should be 'A', 'C', 'G' or 'T'",
            ),
            "to": Or(
                "A", "C", "T", "G", error="Invalid value in config for forward_reads.to, should be 'A', 'C', 'G' or 'T'"
            ),
        },
        "reverse_reads": {
            "from": Or(
                "A",
                "C",
                "T",
                "G",
                error="Invalid value in config for reverse_reads.from, should be 'A', 'C', 'G' or 'T'",
            ),
            "to": Or(
                "A", "C", "T", "G", error="Invalid value in config for reverse_reads.to, should be 'A', 'C', 'G' or 'T'"
            ),
        },
        "delim": Or("\t", ",", ";", error="Invalid value in config for delim, should be '\t', ',' or ';'"),
        "schema": {
            # Whilst strings in default_config.yaml they are updated to the appropriate type on loading
            "read_uid": Or(type(eval("int")), error="Invalid value in config for schema.read_uid, should be 'int'"),
            "transcript_id": Or(
                type(eval("str")), error="Invalid value in config for schema.transcript_id, should be 'str'"
            ),
            "start": Or(type(eval("int")), error="Invalid value in config for schema.read_uid, should be 'int'"),
            "end": Or(type(eval("int")), error="Invalid value in config for schema.read_uid, should be 'int'"),
            "chr": Or(type(eval("str")), error="Invalid value in config for schema.transcript_id, should be 'str'"),
            "strand": Or(type(eval("str")), error="Invalid value in config for schema.strand, should be 'str'"),
            "assignment": Or(type(eval("str")), error="Invalid value in config for schema.assignment, should be 'str'"),
            "conversions": Or(
                type(eval("int")), error="Invalid value in config for schema.conversions, should be 'int'"
            ),
            "convertible": Or(
                type(eval("int")), error="Invalid value in config for schema.convertible, should be 'int'"
            ),
            "coverage": Or(type(eval("int")), error="Invalid value in config for schema.coverage, should be 'int'"),
        },
        "summary_counts": {
            "file_pattern": Or(
                str, error="Invalid value in config for summary_counts.file_pattern, should be str/glob"
            ),
            "separator": Or("\t", ",", ";", error="Invalid value in config for delim, should be '\t', ',' or ';'"),
            "groupby": Or([str], error="Invalid value in config for summary_counts.groupby, should be list of str"),
            "output": {
                "outfile": Or(str, error="Invalid value in config for summary_counts.output.outfile, should be str"),
                "sep": Or("\t", ",", ";", error="Invalid value in config for delim, should be '\t', ',' or ';'"),
                "index": bool,
            },
        },
        "plot": bool,
    }
)
