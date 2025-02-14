# Adding modules

This document describes how to add additional post-processing modules to summarise, model or plot the results. As the
processing pipeline uses Python these steps are undertaken using [Pandas][pandas], [Polars][polars],
[plotnine][plotnine] and [statsmodels][statsmodels] but other frameworks can be used to extend functionality such as
[Scikit-learn][sklearn] (or any other Python module!).

You should set yourself up with a development environment as described in the [contributing](index.md) section so
that linting and pre-commit hooks will run locally, shortening the development feedback loop.

## Organising Modules

Typically the functionality that is likely to be added will be some aspect that summarises or plots the post-processing
results. To which end there is the [`isoslam.summary`][summary] and [`isoslam.plotting`][plotting] modules to which
functionality should be added. There is no harm in adding a new module, particularly if there are a large number of
functions that are required but you may wish to consider adding the entry point function to one of these modules.

### Parameters

The parameters should be clearly defined with [typehints][typing] as the [pre-commit][pre-commit] hooks will fail (if
not locally then in Continuous Integration which blocks merging). The nomenclature used for parameters should be the
basis of configuration options used (see next section). This makes it possible to leverage [`**kwargs`][kwargs] to pass
options loaded from the configuration dictionary and updated from the command line, through to the functions.

## Configuration

Configuration options should be added to `isoslam/default_config.yaml`. A section should be defined for the module you
are adding, in this worked example we are creating the `plot_conversions` and so we would add a section that
corresponds to the arguments required for plotting. The function name is the top-level and options for this module are
nested within.

```yaml
plot_conversions:
  group_by: "read"
  theme: "classic"
```

### Validation

There is a validation module in place [`isoslam.validation`][validation] which checks that the parameters in the
`default_config.yaml`, a user supplied configuration or command line options are of the expected type. You need to add
the options you have added to `isoslam/default_config.yaml` to the `DEFAULT_CONFIG_SCHEMA` that is defined in the
[`isoslam.validation`][validation] module. The examples there should be informative for writing/adding new dictionary
entries. The keys are the fields expected in the configuration, the values are the expected types or the `schema.Or()`
function which states the type(s)/values that are permitted and lists an `error="<error message>"` that is displayed if
the condition is not met. For the above additional configuration you would add the following to the
`DEFAULT_CONFIG_SCHEMA`, nesting the options as reflected in the configuration structure.

```python
DEFAULT_CONFIG_SCHEMA = Schema(
    {
        "plot_conversions": {
            "group_by": Or(
                "read",
                "pair",
                error="Invalid value in config for plot.conversions.group_by, valid values are 'read' or 'pair'",
            ),
            "theme": Or(
                "classic",
                "bw",
                error="Invalid value in config for plot.theme, valid values are 'classic' or 'bw'",
            ),
        }
    }
)
```

## Entry Points

To make the module available at the command line, and in turn possible to integrate into the CGAT pipeline with you need
what is known in Python packaging as an [entry point][setuptools_entrypoint]. This is a method of providing a simple
command line interface to access your program and sub-modules so that they do not need prefixing with `python -m`. The
module where this is setup is [`isoslam.processing`][processing] where you will see there is a `create_config()`
function which creates an argument parser along with sub-parsers. The new function you are adding will be added as a
sub-parser, in the example below we add a sub-parser for plotting the number of conversions per read.

There should be one argument for every configuration option defined in `default_config.yaml`, which in turn mirrors the
options used in the functions that you call, and each `dest` should match these names. The updating of the configuration
based on command line options is contingent on these aligning.

```python
def create_parser() -> arg.ArgumentParser:
    """
    Create a parser for reading options.

    Parser is created with multiple sub-parsers for eading options to run ``isoslam``.

    Returns
    -------
    arg.ArgumentParser
        Argument parser.
    """
    ...

    # Plot conversions per read
    plot_conversions = subparsers.add_parser(
        "plot-conversions-per-read",
        description="Plot the conversions per read.",
        help="Plot the conversions per read.",
    )
    plot_conversions_parser.add_argument(
        "--file-pattern",
        dest="file_pattern",
        type=str,
        required=False,
        default="*_summarized.tsv",
        help="Regular expression for summarized files to process.",
    )
    plot_conversions_parser.add_argument(
        "--outfile",
        dest="outfile",
        type=Path,
        required=False,
        default="conversions.png",
        help="Output filename to save results to, will be nested under 'output_dir'.",
    )
    plot_conversions_parser.add_argument(
        "--separator",
        dest="sep",
        type=str,
        required=False,
        default="\t",
        help="Field separator to use in output file, default is '\t' but other values (e.g. ',' are allowed).",
    )
    plot_conversions.set_defaults(func=plot_conversions)
```

This sets up the subparser `plot_conversions_parser` which has three optional arguments to specify the
`file_pattern` to be searched for and subsequently loaded, the `outfile` name which will be nested under the
`output_dir` (which is an argument to the `isoslam` entry point) and the `separator` that is used in the files. Finally
a default function is set, in this case `plot_conversions`.

The `plot_conversions()` function, which we will define within the processing module does the work of calling
the module you have added. The only argument it needs is `args` which will be the `arguments.Namespace` that is created
by the argument parser. These are used to update the default options which read from the `isoslam/default_config.yaml`
with values the user enters which is the validated with a call to `validation.validate_config()`.

The task of plotting conversions requires that we first load a series of files, combine them and summarise
them which are the first set of steps taken, including some subsetting of configuration options. This data is then
summarized and plotted using the functions defined and imported from the `isoslam.summary` module and the
`isoslam.plotting` modules.

```python
from isoslam import plotting as plot
from isoslam import summary, validation


def plot_conversions(args: arg.Namespace | None) -> None:
    """
    Take a set of output files and summarise the number of conversions.

    Counts are made within file, chromosome, transcript, start, end, assignment and whether there is one or more
    conversion observed.

    Parameters
    ----------
    args : arg.Namespace | None
        Arguments function was invoked with.

    Returns
    -------
    None
        Function does not return anything.
    """
    config = io.load_and_update_config(args)
    logger.remove()
    if vars(args)["log_level"] is not None:
        logging.setup(level=vars(args)["log_level"])
    else:
        logging.setup(level=config["log_level"])
    # Validate the configuration
    validation.validate_config(
        config=config,
        schema=validation.DEFAULT_CONFIG_SCHEMA,
        config_type="configuration",
    )
    # Load and summarise the data
    plot_conversions_config = config["plot_conversions"]
    output_config = summary_counts_config.pop("output")
    output_config["output_dir"] = config["output_dir"]
    plot.conversions(**plot_conversions_config)
    logger.info(
        f"Conversions per read plotted to : {output_config['output_dir]}/{output_config['outfile]}"
    )
```

[kwargs]: https://realpython.com/python-kwargs-and-args/
[pandas]: https://pandas.pydata.org/
[plotnine]: https://plotnine.org/
[polars]: https://pola.rs/
[pre-commit]: https://pre-commit.com/
[sklearn]: https://scikit-learn.org/stable/
[setuptools_entrypoint]: https://setuptools.pypa.io/en/latest/userguide/entry_point.html
[statsmodels]: https://www.statsmodels.org/stable/index.html
[typing]: https://docs.python.org/3/library/typing.html
