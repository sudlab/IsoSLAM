"""Functions for summarising output."""

from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from isoslam import io

GROUPBY_FILENAME = ["Transcript_id", "Strand", "Start", "End", "Assignment", "filename"]
GROUPBY_DAY_HR_REP = ["Transcript_id", "Strand", "Start", "End", "Assignment", "day", "hour", "replicate"]

# pylint: disable=too-many-arguments
# pylint: disable=too-many-positional-arguments


def append_files(file_ext: str = ".tsv", directory: str | Path | None = None) -> pl.DataFrame:
    """
    Append a set of files into a Polars DataFrame.

    Parameters
    ----------
    file_ext : str
        File extension to search for results to summarise.
    directory : str | Path | None
        Path on which to search for files with ``file_ext``, if ``None`` then current working directory is used.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrames of each file found.
    """
    _data = io.load_output_files(file_ext, directory)
    all_data = [data.with_columns(filename=pl.lit(key)) for key, data in _data.items()]
    return pl.concat(all_data)


def summary_counts(
    file_ext: str = ".tsv",
    directory: str | Path | None = None,
    groupby: list[str] | None = None,
    conversions_var: str = "Conversions",
    conversions_threshold: int = 1,
    test_file: str | None = "no4sU",
    filename_var: str | None = None,
    regex: str | None = None,
) -> pl.DataFrame:
    r"""
    Group the data and count by various factors.

    Typically though we want to know whether conversions have happened or not and this is based on the ``Conversions  >=
    1``, but this is configurable via the ``conversions_var`` and ``conversions_threshold`` parameters.

    Parameters
    ----------
    file_ext : str
        File extension to search for results to summarise.
    directory : str | Path | None
        Path on which to search for files with ``file_ext``, if ``None`` then current working directory is used.
    groupby : list[str]
        List of variables to group the counts by, if ``None`` then groups the data by ``Transcript_id``,
        ``Strand``, ``Start``, ``End``, ``Assignment``, and   ``filename``.
    conversions_var : str
        The column name that holds conversions, default ``Conversions``.
    conversions_threshold : int
        Threshold for counting conversions, default ``1``.
    test_file : str | None
        Unique identifier for test file, files with this string in their names are removed.
    filename_var : str | NOne
        Column that holds filename.
    regex : str
        Regular expression pattern to extract the hour and replicate from, default ``r"^d(\w+)_(\w+)hr(\w+)_"``.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame counting the total conversions, number by whether conversions happened and the percentage.
    """
    if groupby is None:
        groupby = GROUPBY_FILENAME
    if filename_var is None:
        filename_var = "filename"
    if regex is None:
        regex = r"^d(\w+)_(\w+)hr(\w+)_"
    df = append_files(file_ext, directory)
    if test_file is not None:
        df = df.filter(pl.col(filename_var) != test_file)
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
    df_count_conversions = extract_day_hour_and_replicate(df_count_conversions, filename_var, regex)
    # Sort the data and remove tests (where day is null)
    sort = groupby + ["day", "hour", "replicate", "one_or_more_conversion"]
    df_count_conversions = df_count_conversions.sort(sort, maintain_order=True)
    return df_count_conversions.filter(~pl.col("day").is_null())


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
        The name of the column that holds the filename, default ``filename``.
    regex : str
        Regular expression pattern to extract the hour and replicate from, default ``r"^d(\w+)_(\w+)hr(\w+)_"``.

    Returns
    -------
    pl.DataFrame
        Polars DataFrame augmented with the hour and replicate extracted from the filename.
    """
    return df.with_columns(
        (pl.col(column).str.extract(regex, group_index=1).str.to_integer(strict=False).alias("day")),
        (pl.col(column).str.extract(regex, group_index=2).str.to_integer(strict=False).alias("hour")),
        (pl.col(column).str.extract(regex, group_index=3).str.to_integer(strict=False).alias("replicate")),
    )


def _aggregate_conversions(
    df: pl.DataFrame, groupby: list[str] | None = None, converted: str | None = "one_or_more_conversion"
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
    # Making sure this was correct caused @ns-rse quite a few headaches as initially it appeared that the sorting was
    # retained from earlier steps but that True < False!
    sortby = groupby.copy()
    sortby.append(converted)  # type: ignore[arg-type]
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
        Whether the function is being tested or not. This will prevent a call to ``_aggregate_conversions()`` to
        aggregate the input and simply filter the data.

    Returns
    -------
    pl.DataFrame
        Aggregated dataframe.
    """
    if not test:
        df = _aggregate_conversions(df, groupby, converted)
    # pylint: disable=singleton-comparison
    return df.filter((pl.col("len") == 1) & (pl.col(converted) == False)).drop("len")  # noqa: E712


def _get_one_or_more_conversion(
    df: pl.DataFrame, groupby: list[str] | None = None, converted: str = "one_or_more_conversion"
) -> pl.DataFrame:
    """
    Extract instances where one or more conversion has occurred.

    There are some cases where this isn't the case and for a given subset the ``converted`` variable, which indicates if
    one or more conversion has occurred will only be ``False`` For such instances dummy entries are created based on the
    ``groupby`` variable and appended to the subset of instances where this one or more conversions have been observed.

    This function takes as input the results of ``summary_count()`` it will not work with intermediate files.

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
    no_conversions = df.join(no_conversions, on=groupby, how="inner", maintain_order="left")
    no_conversions = no_conversions.with_columns(
        conversion_count=0,
        one_or_more_conversion=True,
        conversion_percent=0.0,
    )
    no_conversions = no_conversions.with_columns(pl.col("conversion_count").cast(pl.UInt32))
    df = pl.concat([df, no_conversions.select(df.columns)])
    keep = groupby + ["conversion_count", "conversion_total", "conversion_percent"]
    # pylint: disable=singleton-comparison
    return df.filter(pl.col(converted) == True).select(keep).sort(groupby)  # noqa: E712


def _percent_conversions_across_replicates(
    df: pl.DataFrame, groupby: list[str] | None, count: str = "conversion_count", total: str = "conversion_total"
) -> pl.DataFrame:
    """
    Percentage of conversions across replicates for each time point.

    The raw counts and total conversions for each replicate are available. These are summed and the percentage of
    conversions across replicates calculated. This is mathematically the same as taking the weighted mean of the
    percentage of conversions within each replicate.

    Parameters
    ----------
    df : pl.DataFrame
        Polars Dataframe of conversions.
    groupby : list[str], optional
        Variables to ``group_by`` the data, default is ``transcript_id, start, end, assignment, day, hour``.
    count : str
        Variable/column name holding the counts, default is ``conversion_count``.
    total : str
        Variable/column name holding the total number of conversions, default is ``conversion_total``.

    Returns
    -------
    pl.DataFrame
        Weighted mean of the percentage of conversions (weighted by total conversions) across replicates for the given
        transcript/assignment/strand/day/hour (as specified by ``groupby``).
    """
    if groupby is None:
        groupby = ["Transcript_id", "Strand", "Start", "End", "Assignment", "day", "hour"]
    _keep = groupby + [count, total]
    return (
        df.select(_keep)
        .group_by(groupby, maintain_order=True)
        .agg([pl.col(count).sum(), pl.col(total).sum()])
        .with_columns(((pl.col(count) / pl.col(total)) * 100).alias("conversion_percent"))
    )


def _select_base_levels(df: pl.DataFrame, base_day: int = 0, base_hour: int = 0) -> pl.DataFrame:
    """
    Select the base level reference across all data.

    This allows selecting the base level of totals and percents which are used for normalising values. Will drop the
    column ``replicate`` from the data frame.

    Parameters
    ----------
    df : pl.DataFrame
        Polars Dataframe of conversions.
    base_day : int
        Day to be used for reference, default is ``0`` and is unlikely to need changing.
    base_hour : int
        Hour to be used for reference, default is ``0`` and is unlikely to need changing.

    Returns
    -------
    pl.DataFrame
        Subset of data with values at baseline (default ``day == 0 & hour == 0``).
    """
    return (
        df.select(pl.all().name.map(lambda col_name: col_name.replace("conversion", "baseline")))
        .filter((pl.col("day") == base_day) & (pl.col("hour") == base_hour))
        .drop(["day", "hour"])
    )


def _merge_average_with_baseline(
    df_average: pl.DataFrame,
    df_baseline: pl.DataFrame,
    join_on: list[str] | None = None,
    remove_zero_baseline: bool = True,
) -> pl.DataFrame:
    """
    Merge a data frame with the baseline measurements.

    Typically for this workflow this involves merging the average data frame (across replicates at each of the
    transcripts/start/end/strand/assignments) with the average at the baseline to allow normalising the data.

    Parameters
    ----------
    df_average : pl.DataFrame
        Polars Dataframe of averaged data.
    df_baseline : pl.DataFrame
        Polars Dataframe of averaged baseline data.
    join_on : list[str] | None
        Variables to join the data frames on, if ``None`` (default) it is set to ``Transcript_id, Start, End,
        Assignment, Strand``.
    remove_zero_baseline : bool
        Remove instances where the baseline percentage conversion is zero.

    Returns
    -------
    pl.DataFrame
        Averaged and baseline data frame merged on ``join_on``.
    """
    if join_on is None:
        join_on = ["Transcript_id", "Start", "End", "Assignment", "Strand"]
    if remove_zero_baseline:
        df_baseline = df_baseline.filter(pl.col("baseline_percent") != 0.0)
    return df_average.join(df_baseline, on=join_on)


def _derive_weight_within_isoform(
    df: pl.DataFrame,
    groupby: list[str] | None,
    total: str = "conversion_total",
) -> pl.DataFrame:
    """
    Calculate weighting used for normalised percentages within each isoform across all time points.

    Where the number of total reads (across replications) is higher then we are more confident in the percentage of
    conversions observed and so we weight the percentages at each time point by the proportion of total counts which
    were calculated previously when deriving the percentage of conversions across replicates (with the
    ``_percent_conversions_across_replicates()`` function).

    Parameters
    ----------
    df : pl.DataFrame
        Dataframe for which weights are to be derived.
    groupby : list[str]
        Grouping for summation of total counts, defaults to ``["Transcript_id", "Strand", "Start", "End",
        "Assignment"]``.
    total : str
        Variable that nolds the total number of conversions (across all replicates), default is ``conversion_total`` and
        shouldn't need changing.

    Returns
    -------
    pl.DataFrame
        DataFrame with two new columns, the sum of total conversions across replicates and time points
        (``conversion_total_all_time_points``) and the weight of conversions at each time point (``conversion_weight``).
    """
    groupby = ["Transcript_id", "Strand", "Start", "End", "Assignment"] if groupby is None else groupby
    counts_across_isoform = df.group_by(groupby).agg([pl.col(total).sum().alias("conversion_total_all_time_points")])
    df = df.join(counts_across_isoform, on=groupby, how="inner")
    return df.with_columns((pl.col(total) / pl.col("conversion_total_all_time_points")).alias("conversion_weight"))


def _normalise(
    df: pl.DataFrame,
    to_normalise: str = "conversion_percent",
    baseline: str = "baseline_percent",
    normalised: str = "normalised_percent",
) -> pl.DataFrame:
    """
    Normalise variables based on the baseline measurement.

    Assumes that you have merged the averaged dataset with the averaged baseline variables so that the parameter of
    interest as its related baseline measurement paired with it. Values are normalised by dividing by the baseline value
    such that baseline will always start at ``1`` and subsequent values (time-points) are relative to this and show
    increases or decreases. Typically these will be relative changes in the (averaged) percentage of conversions.

    Parameters
    ----------
    df : pl.DataFrame
        Dataframe from ``_merge_average_with_baseline``.
    to_normalise : str
        Variable to be normalised, default is ``conversion_percent``.
    baseline : str
        Variable to use for normalising, default is ``baseline_percent``.
    normalised : str
        Variable name for normalised value, default is ``normalised_percent``.

    Returns
    -------
    pl.DataFrame
        Polars dataframe with normalised values.
    """
    return df.with_columns([(pl.col(to_normalise) / pl.col(baseline)).alias(normalised)])


def _find_read_pairs(
    df: pl.DataFrame, index_columns: set[str] | None = None, assignment: str | None = "Assignment"
) -> pl.DataFrame:
    """
    Find instances where there are conversions for both ``Return`` and ``Splice`` assignments.

    Parameters
    ----------
    df : pl.DataFrame
        Polars DataFrame.
    index_columns : set
        List of index columns to select from the dataframe. Should include the unique identifiers, typically
        (``Transcript_id``, ``Strand``, ``Start`` and ``End`` which are the defaults) but does not need to include the
        ''assignment'' column.
    assignment : str
        Column the defines assignment of events to ``Ret`` (``Return``) or ``Spl`` (``Splice``).

    Returns
    -------
    pl.DataFrame
        Polars DataFrame of the ``index_columns`` where both a ``Ret`` and ``Spl`` event have been observed.
    """
    if assignment is None:
        assignment = "Assignment"
    if index_columns is None:
        index_columns = {"Transcript_id", "Strand", "Start", "End"}
    index_columns.add(assignment)
    df_return = df.select(index_columns).filter(pl.col(assignment) == "Ret")
    df_splice = df.select(index_columns).filter(pl.col(assignment) == "Spl")
    index_columns.remove(assignment)
    # We use sorted(index_columns, reverse=True) so that the order is consistent for testing, the reverse option roughly
    # gets things close to the expected order of columns used in the data.
    return (
        df_return.join(df_splice, on=index_columns, how="inner")
        .select(sorted(index_columns, reverse=True))
        .unique()
        .sort(by=sorted(index_columns, reverse=True))
    )


# mypy: disable-error-code="no-redef"


@dataclass()
class Statistics:  # pylint: disable=too-many-instance-attributes
    """Staistical summary of results."""

    # Initialised attributes
    file_ext: str
    directory: str | Path
    groupby: list[str] | None
    conversions_var: str | None
    conversions_threshold: int
    test_file: str | None
    regex: str | None

    # Generated atrtibute
    data: pl.DataFrame = field(init=False)
    averages: pl.DataFrame = field(init=False)
    baseline: pl.DataFrame = field(init=False)
    normliased: pl.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        """After initialisation the files are loaded and prepared for analysis."""
        self.data = summary_counts(
            file_ext=self._file_ext,
            directory=self._directory,
            regex=self._regex,
            groupby=self._groupby,
            conversions_var=self._conversions_var,
            conversions_threshold=self._conversions_threshold,
            test_file=self._test_file,
        )
        _df = _aggregate_conversions(self.data, self.groupby, self._conversions_var)
        _df = _filter_no_conversions(_df, self.groupby, self._conversions_var, test=False)
        _df = _get_one_or_more_conversion(_df, self.groupby, self._conversions_var)
        self.averages = _percent_conversions_across_replicates(_df, self.groupby)
        self.baseline = _select_base_levels(self.averages)
        self.normalised = _merge_average_with_baseline(self.averages, self.baseline)
        # Normalise mean conversion percent change by baseline
        self.normalised = _normalise(
            self.normalised, to_normalise="conversion_percent", baseline="baseline_percent", normalised="normalised"
        )
        # Derive weights within transcript/isoform based on total counts
        self.normalised = _derive_weight_within_isoform(self.normalised, groupby=None, total="conversion_total")

    @property
    def file_ext(self) -> str:
        """
        Getter method for ``file_ext``.

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
        Getter method for ``directory``.

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
    def regex(self) -> str:
        """
        Getter method for ``regex``.

        Returns
        -------
        str
            Regex for extracting day/hour/replication from filename.
        """
        return self._regex

    @regex.setter
    def regex(self, value: str) -> None:
        """
        Setter for regex used to extract day/hour/replication from filename..

        Parameters
        ----------
        value : str
            Regex to use for extracting day/hour/replication from filename.
        """
        self._regex = value

    @property
    def groupby(self) -> list[str]:
        """
        Getter method for ``groupby``.

        Returns
        -------
        list[str]
            List of variables to groupby.
        """
        return self._groupby

    @groupby.setter
    def groupby(self, value: list[str]) -> None:
        """
        Setter for the ``groupby`` property.

        Parameters
        ----------
        value : list[str]
            Variables to group data by.
        """
        self._groupby = value

    @property
    def conversions_var(self) -> str:
        """
        Getter method for ``conversions_var``.

        Returns
        -------
        str
            The conversions variable.
        """
        return self._conversions_var

    @conversions_var.setter
    def conversions_var(self, value: str) -> None:
        """
        Setter for the ``conversions_var`` property.

        Parameters
        ----------
        value : list[str]
            Variables to group data by.
        """
        self._conversions_var = value

    @property
    def conversions_threshold(self) -> int:
        """
        Getter method for ``conversions_threshold``.

        Returns
        -------
        int
            The conversion threshold for counting.
        """
        return self._conversions_threshold

    @conversions_threshold.setter
    def conversions_threshold(self, value: int) -> None:
        """
        Setter for the ``conversions_threshold``.

        Parameters
        ----------
        value : int
            Threshold value for counting conversions.
        """
        self._conversions_threshold = value

    @property
    def test_file(self) -> str:
        """
        Getter method for ``test_file``.

        Returns
        -------
        str
            String pattern of test filename for excluding test file data.
        """
        return self._test_file

    @test_file.setter
    def test_file(self, value: str) -> None:
        """
        Setter for the ``test_file`` value.

        Parameters
        ----------
        value : str
            Value of ``test_file``.
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
            Columns to use for identifying unique observations. If ``None`` defaults to ``filename`` which returns the
            number of unique files loaded from the ``directory`` with ``file_ext``.

        Returns
        -------
        int
            Number of unique rows for the given set of variables.
        """
        columns = ["filename"] if columns is None else columns
        return len(self.data.unique(subset=columns))

    # def aggregate_data() -> None:
    #     pass
