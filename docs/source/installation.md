# Installation

Ideally you should install IsoSLAM under a Python Virtual Environment. Details of how to work with and use these is
beyond the scope of this documentation but some advice can be found in the
[contributing](contributing#virtual-environments) section.

## GitHub

There are two methods of installing IsoSLAM from its [GitHub repository][isoslam].

### Cloning

You can clone the repository and install from the clone.

```bash
git clone git@github.com:sudlab/IsoSLAM.git
cd IsoSLAM
pip install -e .
```

By using the `-e` (editable) flag it means you can switch branches.

### `pip` from GitHub

The package installer for Python [pip][pip] can be used to install packages directly from their version control
homepage.

```bash
pip install git+ssh://git@github.com/IsoSLAM
```

If you want to install a specific branch you can

## PyPI

We intend to publish IsoSLAM to the [Python Package Index (PyPI)][pypi]. When available you will be able to install with

```bash
pip install IsoSLAM
```

**NB** IsoSLAM is **NOT** currently available on PyPI.

[isoslam]: https://github.com/sudlab/IsoSLAM
[pip]: https://pip.pypa.io/en/stable/
[pypi]: https://pypi.org/
