# Extending IsoSLAM

<!-- markdownlint-disable MD024 -->

The modular nature of IsoSLAM and its use of `ruffus` and `cgat` mean that it is relatively straight-forward to add
additional steps or processing.

## Overview

1. Add a new module to `isoslam/<module_name>.py` and write functions, include [numpydoc
   strings](https://numpydoc.readthedocs.io/en/latest/format.html) so the functions/classes are documented.
2. Add all options to `isoslam/default_config.yaml`.
3. Add a `sub-parser` to `isoslam/processing.py` with command line options for all arguments to your function.
4. Add a `process_<module_name>` function to `isoslam/processing.py`.
5. Add an entry to the documentation to build the API documentation automatically.

By way of example the implementation of the `isoslam summary` sub-command is explained.

## Adding a module

This is probably the most flexible part, you can add the module as you see fit. You can use Object Orientated approach
and write a class or classes or functional programming and a series.

However you will need a single function that takes the input and various options.

### Example

The `summary` module appends multiple files produced from running IsoSLAM on a series of inputs and appends the
data. These are then summarised by a set of variables to give the number of counts.

#### `isoslam/summary.py`

```python
"""Functions for summarising output."""

import pandas as pd

from isoslam import io


def append_files(pattern: str = "**/*.tsv", separator: str = "\t") -> pd.DataFrame:
    """
    Append a set of files into a Pandas DataFrames.

    Parameters
    ----------
    pattern : str
        File name pattern to search for.
    separator : str
        Separator/delimiter used in files.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrames of each file found.
    """
    _data = io.load_files(pattern, separator)
    all_data = [data.assign(filename=key) for key, data in _data.items()]
    return pd.concat(all_data)


def summary_counts(
    file_pattern: str = "**/*.tsv",
    separator: str = "\t",
    groupby: list[str] | None = None,
    dropna: bool = True,
) -> pd.DataFrame:
    """
    Count the number of assigned read pairs.

    Groups the data by

    Parameters
    ----------
    file_pattern : str
        File name pattern to search for.
    separator : str
        Separator/delimiter used in files.
    groupby : list[str]
        List of variables to group the counts by.
    dropna : book
        Whether to drop rows with ``NA`` values.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrames of each file found.
    """
    if groupby is None:
        groupby = [
            "Transcript_id",
            "Chr",
            "Strand",
            "Start",
            "End",
            "Assignment",
            "Conversions",
            "filename",
        ]
    _data = append_files(file_pattern, separator)
    _data["one_or_more_conversion"] = _data["Conversions"] >= 1
    groupby.append("one_or_more_conversion")
    return _data.value_counts(subset=groupby, dropna=dropna).reset_index()
```

This included writing a function to search for files with a given `pattern` and load them using the specified
`separator`. As this is an Input/Output operation the functions were added to the `isoslam/io.py` module.

### `io.py`

```python
def _find_files(pattern: str = "**/*.tsv") -> Generator:  # type: ignore[type-arg]
    """
    Find files that match the given pattern.

    Parameters
    ----------
    pattern : str
        Pattern (regular expression) of files to search for.

    Returns
    -------
    Generator[_P, None, None]
        A generator of files found that match the given pattern.
    """
    pwd = Path.cwd()
    return pwd.rglob(pattern)


def load_files(pattern: str = "**/*.tsv", sep: str = "\t") -> dict[str, pd.DataFrame]:
    """
    Read a set of files into a list of Pandas DataFrames.

    Parameters
    ----------
    pattern : str
        File name pattern to search for.
    sep : str
        Separator/delimiter used in files.

    Returns
    -------
    list[pd.DataFrame]
        A list of Pandas DataFrames of each file found.
    """
    return {x.stem: pd.read_csv(x, sep=sep) for x in _find_files(pattern)}
```

## Add options to `isoslam/default_config.yaml`

We want to be consistent across the configuration file, which resides in `isoslam/default_config.yaml` and is used when
generating configurations using `isoslam create-config`. To do so the function parameters, in this example
`summary_counts()`, should be used as entries in the `isoslam/default_config.yaml`.

### Example

The top level of a modules configuration should match the module name, here `summary_counts`. Each entry is a key/value
pair that corresponds to the arguments of the function, and so we have `file_pattern`, `separator`, `groupby` and
`output` with their various options.

```yaml
summary_counts:
  file_pattern: "**/*.tsv"
  separator: "\t"
  groupby:
    - Transcript_id
    - Chr
    - Strand
    - Start
    - End
    - Assignment
    - Conversions
    - filename
  output:
    outfile: summary_counts.tsv
    sep: "\t"
    index: false
```

## Add a sub-parser to `isoslam/processing.py`

The function `create_parser()` is responsible for creating the `isoslam` arguments and sub-parsers and their associated
arguments.

Define a sub-parser and `add_argument()` for each option that is available. Keep names consistent with the arguments of
the main function you have written above and in turn the configuration values in `isoslam/default_config.yaml`.

### Example

```python
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
```

The last line here `.set_defaults(func=summarise_counts)` is the function that will be called when running the
subcommand and corresponds to the function

## Documentation

To have the documentation automatically built from the docstrings you have written for your functions you need to add a
`docs/api/<module_name>.md` that corresponds to each of the modules you have introduced and add a title, short description
and `::: isoslam.<module_name>`. If you have introduced more than one module then you will have to add a corresponding
file for each module/sub-module you have introduced. If these are nested please mirror the nesting structure in the
documentation.

### Example

```markdown
# Summary

::: isoslam.summary
handler: python
options:
docstring_style:
numpy
rendering:
show_signature_annotations: true
```
