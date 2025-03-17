"""Functions for summarising output."""

from dataclasses import dataclass, field
from pathlib import Path
from re import Pattern

import polars as pl

from isoslam import io

GROUPBY_FILENAME = ["Transcript_id", "Strand", "Start", "End", "Assignment", "filename"]
GROUPBY_DAY_HR_REP = ["Transcript_id", "Strand", "Start", "End", "Assignment", "day", "hour", "replicate"]


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
    filename_col: str | None = None,
    regex: Pattern | None = None,
) -> pl.DataFrame:
    r"""
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
    filename_col : str | NOne
        Column that holds filename.
    regex : Pattern
        Regular expression pattern to extract the hour and replicate from, default ''r"^d(\w+)_(\w+)hr(\w+)_"''.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame counting the total conversions, number by whether conversions happened and the percentage.
    """
    if groupby is None:
        groupby = GROUPBY_FILENAME
    if filename_col is None:
        filename_col = "filename"
    if regex is None:
        regex = r"^d(\w+)_(\w+)hr(\w+)_"
    df = append_files(file_ext, directory, columns)
    if test_file is not None:
        df = df.filter(pl.col(filename_col) != test_file)
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
    df_count_conversions = extract_day_hour_and_replicate(df_count_conversions, filename_col, regex)
    # Sort the data and remove tests (where day is null)
    sort = groupby + ["day", "hour", "replicate", "one_or_more_conversion"]
    df_count_conversions = df_count_conversions.sort(sort, maintain_order=True)
    return df_count_conversions.filter(~pl.col("day").is_null())


def extract_day_hour_and_replicate(
    df: pl.DataFrame, column: str = "filename", regex: Pattern = r"^d(\w+)_(\w+)hr(\w+)_"
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
        (pl.col(column).str.extract(regex, group_index=1).str.to_integer().alias("day")),
        (pl.col(column).str.extract(regex, group_index=2).str.to_integer().alias("hour")),
        (pl.col(column).str.extract(regex, group_index=3).str.to_integer().alias("replicate")),
    )


def _aggregate_conversions(
    df: pl.DataFrame, groupby: list[str] | None = None, converted: str = "Converted"
) -> pl.DataFrame:
    """
    Subset data where there have not been one or more conversions.

    NB : This needs a better description, I've failed to capture the essence of what is being done here.

    Parameters
    ----------
    df : pl.DataFrame
        Summary dataframe aggregated to give counts of one or more conversion.
    groupby : list[str], optional
        Variables to group the data by.
    converted : str
        Variable that contains whether conversions have been observed or not.

    Returns
    -------
    pl.DataFrame
        Aggregated dataframe.
    """
    if groupby is None:
        groupby = GROUPBY_DAY_HR_REP
    # Its important to ensure that the data is not just groupby but that within that it is then sorted by the converted
    # variable. This _should_ be the case if being passed data from summary_count() but to make sure we explicitly sort
    # the data so that pl.first(converted) will _always_ get 'False' first if pl.len() == 2
    # Making sure this was correct cause @ns-rse quite a few headaches as initially it appeared that the sorting was
    # retained from earlier steps but that True < False!
    sortby = groupby.copy()
    sortby.append(converted)
    df = df.sort(sortby)
    q = df.lazy().group_by(groupby, maintain_order=True).agg(pl.len(), pl.first(converted))
    non_captured = q.collect()
    return non_captured.sort(groupby)


def _filter_no_conversions(
    df: pl.DataFrame,
    groupby: list[str] | None = None,
    converted: str | None = "one_or_more_conversion",
    test: bool = False,
) -> pl.DataFrame:
    """
    Filter dataframe for instances where only no conversions have been observed.

    NB : This needs a better description, I've failed to capture the essence of what is being done here.

    Parameters
    ----------
    df : pl.DataFrame
        Summary dataframe aggregated to give counts of one or more conversion.
    groupby : list[str], optional
        Variables to group the data by.
    converted : str
        Variable that contains whether conversions have been observed or not.
    test : bool
        Whether the function is being tested or not. This will prevent a call to ''_aggregate_conversions()'' to
        aggregate the input and simply filter the data.

    Returns
    -------
    pl.DataFrame
        Aggregated dataframe.
    """
    if not test:
        df = _aggregate_conversions(df, groupby, converted)
    return df.filter((pl.col("len") == 1) & (pl.col(converted) == False)).drop(  # noqa: E712 # pylint: disable=singleton-comparison
        "len"
    )


def _inner_join_no_conversions(
    df: pl.DataFrame, groupby: list[str] | None = None, converted: str = "Converted", test: bool = False
) -> pl.DataFrame:
    """
    Make a dummy set of data where no conversions are observed setting the count and percent to zero.

    This function takes as input the results of ''summary_count()'' it will not work with intermediate files.

    NB : This needs a better description or renaming of function, I've failed to capture the essence of what is being
    done here.

    Parameters
    ----------
    df : pl.DataFrame
        Summary dataframe aggregated to give counts of one or more conversion.
    groupby : list[str], optional
        Variables to group the data by.
    converted : str
        Variable that contains whether conversions have been observed or not.

    Returns
    -------
    pl.DataFrame
        Aggregated dataframe.
    """
    if groupby is None:
        groupby = GROUPBY_DAY_HR_REP
    no_conversions = _filter_no_conversions(df, groupby, converted)
    groupby.append(converted)
    print(f"\n{groupby=}\n")
    no_conversions = df.join(no_conversions, on=groupby, how="inner", maintain_order="left")
    no_conversions = no_conversions.with_columns(
        conversion_count=0,
        # This is currently  hard coded we need to make it flexible to use the value of converted,
        # pl.cols(converted). Neither of the following work...
        #    pl.col(converted) = True,
        #    pl.lit(True).alias(converted),
        one_or_more_conversion=True,
        conversion_percent=0.0,
    )
    no_conversions = no_conversions.with_columns(pl.col("conversion_count").cast(pl.UInt32))
    try:
        # Not sure we should drop conversion_total, but if retained we would need a way of getting this value on a
        # group_by basis and retaining it
        df = df.drop(["filename", "conversion_total"])
    except:
        pass
    # Order no_conversions
    no_conversions = no_conversions.select(df.columns)
    print(f"\n{df.shape=}\n")
    print(f"\n{df.schema=}\n")
    print(f"\n{no_conversions.shape=}\n")
    print(f"\n{no_conversions.schema=}\n")
    df = pl.concat([df, no_conversions.select(df.columns)])
    return df.filter(pl.col(converted) == True).sort(groupby)  # noqa: E712 # pylint: disable=singleton-comparison


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
    regex: Pattern | None

    # Generated atrtibute
    data: pl.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        """After initialisation the files are loaded and prepared for analysis."""
        self.data = summary_counts(
            file_ext=self._file_ext,
            directory=self._directory,
            columns=self._columns,
            regex=self._regex,
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
    def regex(self) -> Pattern:
        """
        Getter method for ''regex''.

        Returns
        -------
        Pattern
            regex for extracting day/hour/replication from filename.
        """
        return self._regex

    @regex.setter
    def regex(self, value: Pattern) -> None:
        """
        Setter for regex used to extract day/hour/replication from filename..

        Parameters
        ----------
        Pattern
            Regex to use for extracting day/hour/replication from filename.
        """
        self._regex = value

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

    # def aggregate_data() -> None:
    #     pass
