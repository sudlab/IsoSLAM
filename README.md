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

For now you can install `IsoSLAM` directly from this GitHub repository using `pip`

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

## Usage

## How it Works

## Contributing

Contributions are welcome. We have a [Code of Conduct](CODE_OF_CONDUCT.md) that we ask you respect.

If you have bugs or feature requests please [create an issue][isoslam_issue], there are templates for reporting
[bugs][isoslam_bug] and making [feature requests][isoslam_feature].

If you wish to contribute fixes or features yourself you can find detailed information in the
[contributing](docs/src/CONTRIBUTING.md) document (which is rendered on the [website][contributing]).

## Licence

**This software is licensed as specified by the [MIT License](LICENSE).**

## Citation

Please use the [Citation File Format](https://citation-file-format.github.io/) which is available in this repository. A
BibTex or APA formatted citation can be easily accessed from the "_Cite this repository_" link on the right hand side.

[contributing]: https://sudlab.github.io/IsoSLAM/main/contributing.html
[isoslam_issue]: https://github.com/sudlab/IsoSLAM/issues/new/choose
[isoslam_bug]: https://github.com/sudlab/IsoSLAM/issues/new?assignees=&labels=bug&projects=&template=bug_report.yaml&title=%5BBug%5D%3A+
[isoslam_feature]: https://github.com/sudlab/IsoSLAM/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.yaml&title=%5Bfeature%5D+%3A+
[slamseq]: https://www.lexogen.com/slamseq-metabolic-rna-labeling/
