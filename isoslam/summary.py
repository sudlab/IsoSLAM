"""Functions for summarising output."""

from pathlib import Path

import polars as pl

from isoslam import io


def append_files(
    file_ext: str = ".tsv", directory: str | Path | None = None, columns: list[str] | None = None
) -> pl.DataFrame:
    """
    Append a set of files into a Polars DataFrame.

    Parameters
    ----------
    file_ext : str
        File extension to search for results to summarise.
    directory : str | Path | None
        Path on which to search for files with ''file_ext'', if ''None'' then current working directory is used.
    columns : list[str]
        Columns to load from data files.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrames of each file found.
    """
    _data = io.load_output_files(file_ext, directory, columns)
    all_data = [data.with_columns(filename=pl.lit(key)) for key, data in _data.items()]
    return pl.concat(all_data)


def summary_counts(
    file_ext: str = ".tsv",
    directory: str | Path | None = None,
    columns: list[str] | None = None,
    groupby: list[str] | None = None,
) -> pl.DataFrame:
    """
    Count the number of conversions across multiple files.

    Parameters
    ----------
    file_ext : str
        File extension to search for results to summarise.
    directory : str | Path | None
        Path on which to search for files with ''file_ext'', if ''None'' then current working directory is used.
    columns : list[str]
        Columns to load from data files.
    groupby : list[str]
        List of variables to group the counts by, if ''None'' then groups the data by ''Transcript_id'', ''Chr'',
        ''Strand'', ''Start'', ''End'', ''Assignment'', ''Conversions'', and   ''filename''.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame of the number of reads with one or more conversion across multiple files.
    """
    if groupby is None:
        groupby = ["Transcript_id", "Chr", "Strand", "Start", "End", "Assignment", "Conversions", "filename"]
    df = append_files(file_ext, directory, columns)
    # df["one_or_more_conversion"] = df["Conversions"] >= 1
    df = df.with_columns([(pl.col("Conversions") >= 1).alias("one_or_more_conversion")])
    groupby.append("one_or_more_conversion")
    return df.group_by(groupby).len(name="count")
