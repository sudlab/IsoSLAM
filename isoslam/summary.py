"""Functions for summarising output."""

from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from isoslam import io

DEFAULT_GROUPBY = ["Transcript_id", "Strand", "Start", "End", "Assignment", "filename"]


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


def summary_counts(  # pylint: disable=too-many-positional-arguments
    file_ext: str = ".tsv",
    directory: str | Path | None = None,
    columns: list[str] | None = None,
    groupby: list[str] | None = None,
    conversions_var: str = "Conversions",
    conversions_threshold: int = 1,
    test_file: str | None = "no4sU",
) -> pl.DataFrame:
    """
    Group the data and count by various factors.

    Typically though we want to know whether conversions have happened or not and this is based on the ''Conversions  >=
    1'', but this is configurable via the ''conversions_var'' and ''conversions_threshold'' parameters.

    Parameters
    ----------
    file_ext : str
        File extension to search for results to summarise.
    directory : str | Path | None
        Path on which to search for files with ''file_ext'', if ''None'' then current working directory is used.
    columns : list[str]
        Columns to load from data files.
    groupby : list[str]
        List of variables to group the counts by, if ''None'' then groups the data by ''Transcript_id'',
        ''Strand'', ''Start'', ''End'', ''Assignment'', and   ''filename''.
    conversions_var : str
        The column name that holds conversions, default ''Conversions''.
    conversions_threshold : int
        Threshold for counting conversions, default ''1''.
    test_file : str | None
        Unique identifier for test file, files with this string in their names are removed.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame counting the total conversions, number by whether conversions happened and the percentage.
    """
    if groupby is None:
        groupby = DEFAULT_GROUPBY
    df = append_files(file_ext, directory, columns)
    if test_file is not None:
        df = df.filter(pl.col("filename") != test_file)
    df = df.with_columns([(pl.col(conversions_var) >= conversions_threshold).alias("one_or_more_conversion")])
    # Get counts by variables, including one_or_more_conversion
    groupby.append("one_or_more_conversion")
    df_count_conversions = df.group_by(groupby).len(name="conversion_count")
    # Aggregate again ignoring one_or_more_conversion to give total counts at site
    groupby.remove("one_or_more_conversion")
    df_count_total = df.group_by(groupby).len(name="conversion_total")
    # Combine counts and totals and calculate percent
    df_count_conversions = df_count_conversions.join(df_count_total, on=groupby)
    df_count_conversions = df_count_conversions.with_columns(
        (pl.col("conversion_count") / pl.col("conversion_total")).alias("conversion_percent")
    )
    return extract_day_hour_and_replicate(df_count_conversions)


def extract_day_hour_and_replicate(
    df: pl.DataFrame, column: str = "filename", regex: str = r"^d(\w+)_(\w+)hr(\w+)_"
) -> pl.DataFrame:
    r"""
    Extract the hour and replicate from the filename stored in a dataframes column.

    Parameters
    ----------
    df : pl.DataFrame
        Polars DataFrame.
    column : str
        The name of the column that holds the filename, default ''filename''.
    regex : Pattern
        Regular expression pattern to extract the hour and replicate from, default ''r"^d(\w+)_(\w+)hr(\w+)_"''.

    Returns
    -------
    pl.DataFrame
        Polars DataFrame augmented with the hour and replicate extracted from the filename.
    """
    return df.with_columns(
        (pl.col(column).str.extract(regex, group_index=1).alias("day")),
        (pl.col(column).str.extract(regex, group_index=2).alias("hour")),
        (pl.col(column).str.extract(regex, group_index=3).alias("replicate")),
    )


# mypy: disable-error-code="no-redef"


@dataclass()
class Statistics:  # pylint: disable=too-many-instance-attributes
    """Staistical summary of results."""

    # Initialised attributes
    file_ext: str
    directory: str | Path
    columns: list[str] | None
    groupby: list[str] | None
    conversions_var: str | None
    conversions_threshold: int
    test_file: str | None

    # Generated atrtibute
    data: pl.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        """After initialisation the files are loaded and prepared for analysis."""
        self.data = summary_counts(
            file_ext=self._file_ext,
            directory=self._directory,
            columns=self._columns,
            groupby=self._groupby,
            conversions_var=self._conversions_var,
            conversions_threshold=self._conversions_threshold,
            test_file=self._test_file,
        )

    @property
    def file_ext(self) -> str:
        """
        Getter method for ''file_ext''.

        Returns
        -------
        str
            File extension that is loaded.
        """
        return self._file_ext

    @file_ext.setter
    def file_ext(self, value: str) -> None:
        """
        Setter for the file extension.

        Parameters
        ----------
        value : str
            File extension to load data.
        """
        self._file_ext = value

    @property
    def directory(self) -> str:
        """
        Getter method for ''directory''.

        Returns
        -------
        str
            Directory from which output files are loaded.
        """
        return self._directory

    @directory.setter
    def directory(self, value: str) -> None:
        """
        Setter for the file extension.

        Parameters
        ----------
        value : str
            Directory from which files are loaded.
        """
        self._directory = value

    @property
    def columns(self) -> list[str]:
        """
        Getter method for ''columns''.

        Returns
        -------
        list[str]
            List of columns that are loaded from output.
        """
        return self._columns

    @columns.setter
    def columns(self, value: list[str]) -> None:
        """
        Setter for the file extension.

        Parameters
        ----------
        value : str
            Columns that to load from the data.
        """
        self._columns = value

    @property
    def groupby(self) -> list[str]:
        """
        Getter method for ''groupby''.

        Returns
        -------
        list[str]
            List of variables to groupby.
        """
        return self._groupby

    @groupby.setter
    def groupby(self, value: list[str]) -> None:
        """
        Setter for the file extension.

        Parameters
        ----------
        value : list[str]
            Variables to group data by.
        """
        self._groupby = value

    @property
    def conversions_var(self) -> str:
        """
        Getter method for ''conversions_var''.

        Returns
        -------
        str
            The conversions variable.
        """
        return self._conversions_var

    @conversions_var.setter
    def conversions_var(self, value: str) -> None:
        """
        Setter for the file extension.

        Parameters
        ----------
        value : list[str]
            Variables to group data by.
        """
        self._conversions_var = value

    @property
    def conversions_threshold(self) -> int:
        """
        Getter method for ''conversions_threshold''.

        Returns
        -------
        int
            The conversion threshold for counting.
        """
        return self._conversions_threshold

    @conversions_threshold.setter
    def conversions_threshold(self, value: int) -> None:
        """
        Setter for the ''conversions_threshold''.

        Parameters
        ----------
        value : int
            Threshold value for counting conversions.
        """
        self._conversions_threshold = value

    @property
    def test_file(self) -> str:
        """
        Getter method for ''test_file''.

        Returns
        -------
        str
            String pattern of test filename for excluding test file data.
        """
        return self._test_file

    @test_file.setter
    def test_file(self, value: str) -> None:
        """
        Setter for the ''test_file'' value.

        Parameters
        ----------
        value : str
            Value of ''test_file''.
        """
        self._test_file = value

    @property
    def shape(self) -> tuple[int, int]:
        """
        Getter for the shape of the dataframe.

        Returns
        -------
        tuple[int, int]
            Shape of the Polars dataframe.
        """
        return self.data.shape  # type: ignore[no-any-return]

    @property
    def unique(self) -> int:
        """
        Getter for the number of unique files loaded.

        Returns
        -------
        int
            Number of unique rows.
        """
        return self.unique_rows()

    def unique_rows(self, columns: list[str] | None = None) -> int:
        """
        Identify unique rows in the data for a given set of columns.

        Parameters
        ----------
        columns : list[str]
            Columns to use for identifying unique observations. If ''None'' defaults to ''filename'' which returns the
            number of unique files loaded from the ''directory'' with ''file_ext''.

        Returns
        -------
        int
            Number of unique rows for the given set of variables.
        """
        columns = ["filename"] if columns is None else columns
        return len(self.data.unique(subset=columns))
