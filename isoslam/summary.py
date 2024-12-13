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
        groupby = ["Transcript_id", "Chr", "Strand", "Start", "End", "Assignment", "Conversions", "filename"]
    _data = append_files(file_pattern, separator)
    _data["one_or_more_conversion"] = _data["Conversions"] >= 1
    groupby.append("one_or_more_conversion")
    return _data.value_counts(subset=groupby, dropna=dropna).reset_index()
