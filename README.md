# IsoSLAM

<div align="center">

<!-- [![PyPI version](https://badge.fury.io/py/isoslam.svg)](https://badge.fury.io/py/isoslam) -->
<!-- ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/isoslam) -->

[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Code style: flake8](https://img.shields.io/badge/code%20style-flake8-456789.svg)](https://github.com/psf/flake8)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

<!-- [![codecov](https://codecov.io/gh/sudlab/IsoSLAM/branch/dev/graph/badge.svg)]
(https://codecov.io/gh/sudlab/IsoSLAM) -->
<!-- [![pre-commit.ci -->
<!-- status](https://results.pre-commit.ci/badge/github/sudlab/IsoSLAM/main.svg)]
(https://results.pre-commit.ci/latest/github/sudlab/IsoSLAM/main) -->

</div>

<!-- <div align="center"> -->

<!-- [![Downloads](https://static.pepy.tech/badge/isoslam)](https://pepy.tech/project/isoslam) -->
<!-- [![Downloads](https://static.pepy.tech/badge/isoslam/month)](https://pepy.tech/project/isoslam) -->
<!-- [![Downloads](https://static.pepy.tech/badge/isoslam/week)](https://pepy.tech/project/isoslam) -->

<!-- </div> -->

IsoSLAM is a Python package for processing and working with [SLAMSeq][slamseq] data of RNA expression.

## Installation

IsoSLAM is available on the [Python Package Index (PyPI)][pypi]. It is recommended to use a Python Virtual Environment
and install within that.

```bash
pip install isoslam
```

You can also install `IsoSLAM` directly from this GitHub repository using `pip`...

```bash
pip install git+https://github.com/sudlab/IsoSLAM.git
```

Alternatively you can clone the repository and install it from there, optionally using the `-e` flag to make the code
editable should you wish to work on developing the code base.

```bash
cd ~/path/to/clone/to
git clone git@github.com:sudlab/IsoSLAM.git
cd IsoSLAM
pip install -e .
```

If you wish to include optional dependencies you may do so.

```bash
pip install -e .[dev,docs,tests]
```

### Conda/BioConda

If you use [Conda][conda] or [BioConda][bioconda] to manage your virtual environments and wish to document the
installation of IsoSLAM in a YAML file you can add it as a `pip` dependency as shown in the sample `isoslam.yaml` file
below which includes [cgatcore][cgat] and [ruffus][ruffus] as dependencies which will be installed from one of the
listed Conda channels.

```yaml
name: isoslam

channels:
  - conda-forge
  - bioconda
  - default

dependencies:
  - cgatcore
  - ruffus
  - pip
  - pip:
      - isoslam
```

You can then create an environment using the following.

```bash
conda env create --name isoslam --file isoslam.yaml
```

## Usage

On installation the `isoslam` entry point will be added to the `$PATH` of your Virtual Environment. You can then invoke
it with the `--help` flag to see the available options.

```bash
❱ isoslam --help
usage: isoslam [-h] [-v] [-c CONFIG_FILE] [-b BASE_DIR] [-o OUTPUT_DIR] [-l LOG_LEVEL] {process,create-config,summary-counts} ...

Run various programs related to IsoSLAM. Add the name of the program you wish to run.

options:
  -h, --help            show this help message and exit
  -v, --version         Report the installed version of IsoSLAM.
  -c, --config-file CONFIG_FILE
                        Path to a YAML configuration file.
  -b, --base-dir BASE_DIR
                        Base directory to run isoslam on.
  -o, --output-dir OUTPUT_DIR
                        Output directory to write results to.
  -l, --log-level LOG_LEVEL
                        Logging level to use, default is 'info' for verbose output use 'debug'.

program:
  Available programs listed below:

  {process,create-config,summary-counts}
    process             Process all files and run all summary plotting and statistics.
    create-config       Create a configuration file using the defaults.
    summary-counts      Summarise the counts.
```

Each sub-command includes help too, for example...

```bash
❱ isoslam process --help
usage: isoslam process [-h] [-b BAM_FILE] [-g GTF_FILE] [-d BED_FILE] [-v VCF_FILE] [-u UPPER_PAIRS_LIMIT] [-f FIRST_MATCHED_LIMIT] [--delim DELIM]
                       [--output-file OUTPUT_FILE]

Process all files and run all summary plotting and statistics.

options:
  -h, --help            show this help message and exit
  -b, --bam-file BAM_FILE
                        Path to '.bam' file that has undergone read assignment with 'featureCount'.
  -g, --gtf-file GTF_FILE
                        Path to '.gtf' transcript assembly file.
  -d, --bed-file BED_FILE
                        Path to '.bed' utron file. Must be bed6 format.
  -v, --vcf-file VCF_FILE
                        Path to '.vcf.gz' file.
  -u, --upper-pairs-limit UPPER_PAIRS_LIMIT
                        Upper limit of pairs to be processed.
  -f, --first-matched-limit FIRST_MATCHED_LIMIT
                        Limit of matches.
  --delim DELIM         Delimiter to use in output.
  --output-file OUTPUT_FILE
                        File to write results to.
```

### Including in Pipelines

Calls to this entry point can easily be incorporated into CGAT/`ruffus` or other pipelines.

## How it Works

## Contributing

Contributions are welcome. We have a [Code of Conduct](CODE_OF_CONDUCT.md) that we ask you respect.

If you have bugs or feature requests please [create an issue][isoslam_issue], there are templates for reporting
[bugs][isoslam_bug] and making [feature requests][isoslam_feature].

If you wish to contribute fixes or features yourself you can find detailed information in the
[contributing](docs/src/CONTRIBUTING.md) document (which is rendered on the [website][website]).

## Licence

**This software is licensed as specified by the [MIT License](LICENSE).**

## Citation

Please use the [Citation File Format](https://citation-file-format.github.io/) which is available in this repository. A
BibTex or APA formatted citation can be easily accessed from the "_Cite this repository_" link on the right hand side.

[bioconda]: https://bioconda.github.io/
[cgat]: https://cgat-developers.github.io/cgat-core/
[conda]: https://docs.conda.io/en/latest/
[isoslam_issue]: https://github.com/sudlab/IsoSLAM/issues/new/choose
[isoslam_bug]: https://github.com/sudlab/IsoSLAM/issues/new?assignees=&labels=bug&projects=&template=bug_report.yaml&title=%5BBug%5D%3A+
[isoslam_feature]: https://github.com/sudlab/IsoSLAM/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.yaml&title=%5Bfeature%5D+%3A+
[pypi]: https://pypi.org/
[ruffus]: http://www.ruffus.org.uk/
[slamseq]: https://www.lexogen.com/slamseq-metabolic-rna-labeling/
[website]: https://sudlab.github.io/IsoSLAM/
