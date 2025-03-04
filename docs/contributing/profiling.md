# Profiling

The speed of execution of `isoslam process` can be investigated using Profiling which analyses which functions and
methods take up the most processing time. In order to undertake profiling the [`cProfile`][cprofile] standard library
can be used. Visualisation of the results can be aided using the [SnakeViz][snakeviz] package.

## Performing Profiling

You need a set of sample files to process in order to undertake profiling. Here we use the sample files that are
included as part of the test suite used in development that can be found under `tests/resources`. If you do not yet have
these locally you should `git clone` the [repository][isoslam] and install it in a clean virtual environment with the
development dependencies.

```shell
git clone git@github.com:sudlab/IsoSLAM.git
cd IsoSLAM
pip install -e .[dev]
```

We now make a `tmp/test-YYYYMMDD` directory and copy the necessary files here.

```shell
mkdir -p tmp/test-$(date +%Y%m%d)                            # Uses the current date
cp -r tests/resources/{bam,gtf,bed,vcf} tmp/test-20250221    # Modify to reflect the current date
```

We can run profiling on these samples using the following which writes the profiling to `isoslam-YYYYMMDD.prof`.

```shell
cd tmp/test-20250221
python -m cProfile -o isoslam-$(date +%Y%m%d).prof $(isoslam process \
   --bam-file bam/sorted_assigned/d0_no4sU_filtered_remapped_sorted.sorted.assigned.bam \
   --gtf-file gtf/test_wash1.gtf \
   --bed-file bed/test_coding_introns.bed \
   --vcf vcf/d0.vcf.gz)
```

You can verify output has been produced using [`parquet-tools`][parquettools] which is part of the `dev` dependencies.

```shell
parquet-tools show output/results.parquet
```

The profiling data should have been written to `isoslam.prof`.

```shell
head isoslam.prof
```

## Visualisation of Profiling

To visualise the results of profiling you can invoke `snakeviz` with the `.prof` file that has been generated.

```shell
snakeviz isoslam-<YYYYMMDD>.prof
```

This should launch a new browser tab with the [icicle][snakevizicicle] where the amount of time spent within a function
is proportional to the size of the bar.

[cprofile]: https://docs.python.org/3/library/profile.html#module-cProfile
[isoslam]: https://github.com/sudlab/IsoSLAM
[parquettools]: https://github.com/ktrueda/parquet-tools
[snakeviz]: https://jiffyclub.github.io/snakeviz/
[snakevizicicle]: https://jiffyclub.github.io/snakeviz/#interpreting-results
