# Installation

Ideally you should install IsoSLAM under a Python Virtual Environment. Details of how to work with and use these is
beyond the scope of this documentation but some advice can be found in the
[contributing](contributing/index.md#virtual-environments) section.

## IsoSLAM

### PyPI

Whilst not complete (what software is?) IsoSLAM is available to install from the [Python Package Index (PyPI)][pypi].

```bash
pip install isoslam
```

#### Conda/BioConda

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

### GitHub

There are two methods of installing IsoSLAM from its [GitHub repository][isoslam]. Which you should use depends on
whether you wish to hack on the code/contribute to the development.

#### `pip` from GitHub

If you only wish to try out the latest features and have no intention of modifying any of the code then you can use the
package installer for Python [pip][pip] to install packages directly from their version control homepage.

```bash
pip install git+ssh://git@github.com/IsoSLAM
```

If you want to install a specific branch or commit you can do so.

```bash
pip install git+ssh://git@github.com/IsoSLAM@<branch-name>
pip install git+ssh://git@github.com/IsoSLAM@<commit-hash>
```

#### Cloning

If you wish to hack at the code and possibly contribute to development then you should clone the repository and install
from there in "editable" mode along with the `dev` / `docs` / `tests` optional dependencies.

```bash
git clone git@github.com:sudlab/IsoSLAM.git
cd IsoSLAM
pip install -e .[dev,docs,tests]
```

By using the `-e` (editable) flag it means you can switch branches or make changes to the code and they will be
reflected in the code base when you run commands.

**NB** - If you are not a member of the [sudlab organisation][sudlab] on GitHub you should first [fork][gh_fork] the
repository and make Pull Requests from your fork to the original repository.

### Bioconda

**NB** IsoSLAM is **NOT** currently available on [Bioconda][bioconda] but it is
[planned](https://github.com/sudlab/IsoSLAM/issues/41).

## Dependencies

There are a number of external dependencies required for running IsoSLAM as part of a [ruffus][ruffus] pipeline, which
is typically essential to prepare the files for processing with IsoSlam.

- `samtools` / `bcftools` are both required and can be downloaded from [htslib][htslib]
- [VarScan][varscan] ([documentation][varscan_docs])
- [subread][subread] ([documentation][subread_docs])

### Indirect Dependencies

The pipeline for running the various steps in processing data rely on the [cgat][cgat] tools. These have two external
dependencies themselves and if you wish to use such a pipeline will have to install these.

- [bedtools][bedtools]
- [wigToBigWig][wigtobigwig] - the [minimal tools][wigtobigwig-min] may suffice.

### GNU/Linux

#### Arch Linux

If you use [Arch Linux][arch] the packages are available in the [Arch Linux User Repository (AUR)][aur]

```bash
mkdir ~/aur && cd ~/aur
git clone https://aur.archlinux.org/htslib.git
git clone https://aur.archlinux.org/samtools.git
git clone https://aur.archlinux.org/bcftools.git
git clone https://aur.archlinux.org/subread.git
git clone https://aur.archlinux.org/bedtools.git
cd htslib
makepkg -sri
cd ../bcftools
makepkg -sri
cd ../samtools
makepkg -sri
cd ../subread
makepkg -sri
cd ../bedtools
makepkg -sri
```

[varscan][varscan] is written in Java, you need to download the [latest
release](https://github.com/dkoboldt/varscan/releases)

You can make a wrapper to run this, assuming you have saved the file to `~/.local/jar/VarScan.v2.4.6.jar` (adjust for
the version you have downloaded), you can create the following short script and make it executable, placing it in your
`$PATH` (the example below uses `~/.local/bin/`)

```bash
#!/bin/bash

java -jar ~/.local/jar/VarScan.v2.4.6.jar
```

Make the file executable and you can then run `varscan`

```bash
chmod 755 ~/.local/bin/varscan
varscan --version
VarScan v2.4.6

***NON-COMMERCIAL VERSION***

USAGE: java -jar VarScan.jar [COMMAND] [OPTIONS]

COMMANDS:
 pileup2snp  Identify SNPs from a pileup file
 pileup2indel  Identify indels a pileup file
 pileup2cns  Call consensus and variants from a pileup file
 mpileup2snp  Identify SNPs from an mpileup file
 mpileup2indel  Identify indels an mpileup file
 mpileup2cns  Call consensus and variants from an mpileup file

 somatic   Call germline/somatic variants from tumor-normal pileups
 mpileup2somatic  Call germline/somatic variants in multi-tumor-normal mpileup (beta feature in v2.4.5)
 copynumber  Determine relative tumor copy number from tumor-normal pileups
 readcounts  Obtain read counts for a list of variants from a pileup file

 filter   Filter SNPs by coverage, frequency, p-value, etc.
 somaticFilter  Filter somatic variants for clusters/indels
 fpfilter  Apply the false-positive filter

 processSomatic  Isolate Germline/LOH/Somatic calls from output
 copyCaller  GC-adjust and process copy number changes from VarScan copynumber output
 compare   Compare two lists of positions/variants
 limit   Restrict pileup/snps/indels to ROI positions
```

#### Gentoo

`samtools` and `bcftools` are available in Portage, to install.

```bash
emerge --sync && emerge -av samtools bcftools bedtools ucsc-genome-browser
```

#### Ubuntu/Debian

All three packages are available for Debian based repositories.

```bash
sudo apt-get update
sudo apt-get install samtools bcftools bedtools

```

#### UCSC Genome Browser

This needs installing from source across all distributions. The [minimal][wigtobigwig-min] should suffice.

```bash
git clone git@github.com:ucscGenomeBrowser/kent-core.git
cd kent-core
sudo make
```

#### Source Install

The releases pages includes instructions on how to build the package from source, but note that you will then have to
manually update the packages when new releases are made.

### Windows

**To be written** at some point.

### OSX

**To be written** at some point.

## Conda

If you don't have the ability to install these programmes at the system level an alternative is to use a Conda
environment.

```bash
conda create -n isoslam python==3.12
conda activate isoslam
conda install mamba
mamba install -c conda-forge -c bioconda bedtools
mamba install -c conda-forge -c bioconda cgat-apps
mamba install -c conda-forge -c bioconda samtools bcftools
mamba install -c conda-forge -c bioconda subread
mamba install -c conda-forge -c bioconda varscan
```

### HPC

**To be written** at some point.

[arch]: https://archlinux.org/
[aur]: https://aur.archlinux.org/
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[bioconda]: https://bioconda.github.io/
[cgat]: https://cgat-developers.github.io/cgat-core/
[conda]: https://docs.conda.io/en/latest/
[gh_fork]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo
[htslib]: https://www.htslib.org/
[isoslam]: https://github.com/sudlab/IsoSLAM
[pip]: https://pip.pypa.io/en/stable/
[pypi]: https://pypi.org/
[ruffus]: http://www.ruffus.org.uk/
[subread]: https://github.com/ShiLab-Bioinformatics/subread
[subread_docs]: https://subread.sourceforge.net/
[sudlab]: https://github.com/sudlab/
[varscan]: https://github.com/dkoboldt/varscan
[varscan_docs]: https://dkoboldt.github.io/varscan/
[wigtobigwig]: https://github.com/ucscGenomeBrowser/kent
[wigtobigwig-min]: https://github.com/ucscGenomeBrowser/kent-core
