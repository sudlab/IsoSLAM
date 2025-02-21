# Usage

## Stand alone usage

Once installed the `isolsam` command is should be available and you can view the options and sub-commands with

```shell
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

Global options that set the configuration file to use, the base and output directory and log-level can be
specified. Users then have a number of programs from IsoSLAM.

- `process`
- `create-config`
- `summary-counts`

### Processing

Help is available on using the `isoslam process` command and can be viewed with the `--help` option.

```shell
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

It is important to remember that the `.bam` file requires pre-processing and instructions on how to do that can be found
below.

Three additional input files are required.

- `.gtf` : DESCRIPTION
- `.bed` : DESCRIPTION
- `.vcf` : DESCRIPTION

By default output is written to the `./output` directory which will be created if it does not exist and the default
output file is `results.parquet` which is in [Parquet][parquet] format. The output directory is configurable, as is the
output file name. If you would rather your output is saved as ASCII delimited file it is easy to do so by specifying the
`--delim ,` (for comma-delimited) or `--delim \t` (for tab-delimited) and using the appropriate file extension in the
`--output-file`. For example the following outputs to `new_results` directory in a file `output_20250213.csv`

```shell
❱ isoslam --output-dir new_results process \
         --bam-file d0_no4sU_EKRN230046545-1A_HFWGNDSX7_L4.star.bam \
         --gtf-file test_wash1.gtf \
         --bed-file test_coding_introns.bed \
         --vcf-file d0.vcf.gz \
         --delim "," \
         --output-file output_20250213.csv
```

### Configuration File

It is possible to place all configuration options in a [YAML][yaml] configuration file. A sample configuration file can
be generated using the `isoslam create-config` command and there are options to specify the output filename, the default
is `config.yaml` in the current directory (see `isoslam create-config --help` for all options).

```shell
❱ isoslam create-config --help
usage: isoslam create-config [-h] [-f FILENAME] [-o OUTPUT_DIR]

Create a configuration file using the defaults.

options:
  -h, --help            show this help message and exit
  -f, --filename FILENAME
                        Name of YAML file to save configuration to (default 'config.yaml').
  -o, --output-dir OUTPUT_DIR
                        Path to where the YAML file should be saved (default './' the current directory).
```

Once created you can edit the `config.yaml` to specify all fields and parameters and then invoke processing using this
file.

```shell
❱ isoslam --config-file config.yaml process
```

Any command line options given such as the `--upper-pairs-limit` or `--first-matched-limit` will override those details
in the custom `config.yaml`. This means it is particularly useful for usage in the full pipeline as `.gtf`, `.bed` and
`.vcf` files are often common and it is the `.bam` input file that changes between runs. Thus you can specify the common
files in your configuration file and then use the `--bam-file <file_path>` to override the field in the configuration
file.

## Preparing files

Assuming there are two `.bam` input files, a `.gtf` and Fasta file with the names shown below the following steps must
be performed on the files prior to processing.

- `d0_no4sU_file1.bam`
- `d0_0hr1_file2.bam`
- `test_wash1.gtf`
- `hg38_noalt.fa`

```shell
# Sort both files via read name (insert _sorted before .bam in file name)
samtools sort -n d0_no4sU_file1.bam -o d0_no4sU_file1_sorted.bam
samtools sort -n d0_0hr1_file2.bam -o d0_ohr1_file2_sorted.bam

# Add XT tag to indicate trnascript reads map to and XS tag for which strand the gene is transcribed from
# The -p flag explicitly says to assume libraries contain paired-end reads
mkdir counts
featureCounts -p -a test_wash1.gtf -o counts/d0_no4sU_counts.txt -T2 -R BAM d0_no4sU_file1_sorted.bam
mv counts/d0_no4sU_file1_sorted.bam.featureCounts.bam d0_no4sU_file1_sorted_assigned.bam
featureCounts -p -a test_wash1.gtf -o counts/d0_0hr1_counts.txt -T2 -R BAM d0_0hr1_file2_sorted.bam
mv counts/d0_0hr1_file2_sorted.bam.featureCounts.bam d0_0hr1_file2_sorted_assigned.bam

# Resort the files
samtools sort -n d0_no4sU_file1_sorted_assigned.bam -o d0_no4sU_file1_sorted_assigned_sorted.bam
samtools sort -n d0_0hr1_file2_sorted_assigned.bam -o d0_0hr1_file2_sorted_assigned_sorted.bam

# Take just the negative control (no4sU) and create a SNP VCF file using Varscan
samtools mpileup -B -A -f hg38_noalt.fa d0_no4sU_file_sorted_assigned_sorted.bam |
  varscan mpileup2snp --variants 1 --output-vcf 1 > snp.vcf &&
  bcftools view snp.vcf -Oz -o snp.vcf.gz

# Index the VCF file for reading by PySam
tabix -p vcf snp.vcf.gz

# Run isoslam on the two files
isoslam process --bam-file d0_no4sU_file1_sorted_assigned_sorted.bam --bed test_coding_introns.bed --gtf test_wash1.gtf \
  --vcf-file snp.vcf.gz
isoslam process --bam-file d0_0hr1_file2_sorted_assigned_sorted.bam --bed test_coding_introns.bed --gtf test_wash1.gtf \
  --vcf-file snp.vcf.gz
```

## Ruffus/CGAT pipeline

Typically `isoslam` is part of a workflow pipeline that processes a large number of files. The example below uses the
[ruffus][ruffus] package to control the pipeline but workflow software could be used such as [Nextflow][nextflow].

Ruffus uses [decorators][python_decorators] that control the input and output of functions at each step. Several are
used in this example pipeline.

- [`@ruffus.follows()`][ruffus_follows] : ensures that the task, or list of tasks specified are run before the current
  function is run. Includes the ability to `mkdir("dir_name")` to create directories if they don't already exist.
- [`@ruffus.transforms()`][ruffus_transforms] : Takes input from a previous task/function, applies a filter which is
  typically a regular expression of some description to identify one or more files, and produces corresponding output
  for each matched file.

This allows multiple files to be processed in parallel on HPC systems using the [Slurm][slurm] workload manager without
having to explicitly specify any Slurm jobs or calls.

EXPAND THIS SECTION WITH WORKED EXAMPLES BASED ON EXISTING SCRIPTS.

[nextflow]: https://www.nextflow.io/docs/latest/index.html
[parquet]: https://parquet.apache.org/docs/file-format/
[python_decorators]: https://realpython.com/primer-on-python-decorators/
[ruffus]: http://www.ruffus.org.uk/
[ruffus_follows]: http://www.ruffus.org.uk/decorators/follows.html#follows
[ruffus_transforms]: http://www.ruffus.org.uk/decorators/transform_ex.html
[slurm]: https://slurm.schedmd.com/documentation.html
[yaml]: https://yaml.org/
