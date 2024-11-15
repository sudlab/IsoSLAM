"""
This pipeline takes in mapped RNA seq reads in `.bam`, sorts and
annotates them, then passes them detects and quantifies the number
of 4sU-linked nucleotide conversions (T>C or A>G), exclusing SNPs.
"""

import os
import sys

import ruffus
from cgatcore import pipeline as P

# Read in pipeline.yml
PARAMS = P.get_parameters(
    [
        "%s/pipeline.yml" % os.path.splitext(__file__)[0],
        "../pipeline.yml",
        "pipeline.yml",
    ]
)


@ruffus.follows(ruffus.mkdir("name_sorted_bams"))
@ruffus.transform(
    "bam/*.bam",
    ruffus.regex("bam/(.+).bam"),
    r"name_sorted_bams/\1.sorted.bam",
)
def name_sort_bams(infile, outfile):
    """
    Takes `.bam` file and sorts it in name order, therefore reads from the same pair appear sequentially
    """
    statement = f"samtools sort -n {infile} -o {outfile}"
    P.run(statement, job_memory="8G")


@ruffus.follows(ruffus.mkdir("read_assignments"))
@ruffus.transform(
    name_sort_bams,
    ruffus.regex("name_sorted_bams/(.+).sorted.bam"),
    r"read_assignments/\1/\1.sorted.assigned.bam",
)
def feature_count_read_assignments(infile, outfile):
    """
    Takes `.bam` file and adds 2 tags to it:
        - XT tag tells us the transcript the reads map to
        - XS tag tells us which strand the gene is transcribed from

    NB The `-p` flag is required for SubRead v2.0.6, may need to make this flexible.
    """
    annotation = PARAMS["transcript_gtf"]
    outdir = os.path.dirname(outfile)
    rename_outfile = outfile.replace(".assigned.bam", "")
    statement = (
        "featureCounts"
        "-p "
        f"-a {annotation} "
        f"-o {outdir}/counts.txt "
        "-T 2 "
        "-R BAM"
        f"{infile} &&"
        "mv %(rename_outfile)s.bam.featureCounts.bam %(outfile)s"
    )
    P.run(statement, job_threads=2, job_memory="8G")


@ruffus.transform(
    "input_bams.dir/*no4sU*.bam",
    ruffus.regex("input_bams.dir/(.+)_no4sU_(.+).bam"),
    r"snp_vcf/\1.vcf.gz",
)
def VCF_from_no4sU_control(infile, outfile):
    """
    Takes the negative control `.bam` file (via the no4sU regex) and creates a SNP VCF file
    using Varscan. It looks for differences between the `.bam` file and the genome `.fasta`
    file to detect differences, where these are always present at the same position (above
    a treshold, see Varscan docs) it is classed as a SNP and included in output VCF file.
    """
    outfile_uncompressed = outfile.replace(".vcf.gz", ".vcf")
    genome_fasta = PARAMS["genome_fasta"]
    statement = (
        f"samtools mpileup -B -A -f {genome_fasta} {infile} | "
        f"varscan mpileup2snp --variants 1 --output-vcf 1 > {outfile_uncompressed} && "
        # f"sort -k1,1 -k2,2n {outfile_uncompressed} -o {outfile_uncompressed} &&"
        f"bcftools view {outfile_uncompressed} -Oz -o {outfile}"
    )
    P.run(statement, job_threads=12, job_memory="1G")


@ruffus.transform(VCF_from_no4sU_control, ruffus.suffix(".vcf.gz"), ".vcf.gz.tbi")
def index_VCF(infile, outfile):
    """
    Indexes the VCF file so it can be read with Pysam, generates the .vcf.gz.tbi file
    """
    statement = f"tabix -p vcf {infile}"
    P.run(statement)


@ruffus.follows(index_VCF, ruffus.mkdir("all_introns_counts_and_info"))
@ruffus.transform(
    feature_count_read_assignments,
    ruffus.regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"),
    ruffus.add_inputs(r"snp_vcf/\2.vcf.gz"),
    r"all_introns_counts_and_info/\1.tsv",
)
def all_introns_counts_and_info(infiles, outfile):
    """
    Takes in the sorted and feature assigned `.bam` file from previous steps and passes them to a
    python script that uses pysam (python wrapper for htslib). This script iterates over the `.bam`
    file pair-by-pair (representing the sequencing insert), determines whether the read-pair shows
    evidence of intron splicing/retention and assigns these to specific events by referencing the
    `.gtf` and `.bed` files, and XT tag from feature_count_read_assignments. Next, the script uses the
    XS tag from featureCountsReadAssignment to assign each read in the pair as the forward or reverse
    read, relative to the direction of transcription. Finally, it looks for T>C in FW read, A>G in RV
    read, checks these are not present in the SNP VCF file, and outputs metadata on each read-pair
    about it's event assignment, number of conversions, coverage etc.
    """
    bamfile, vcffile = infiles
    annotation = PARAMS["transcript_gtf"]
    utron_bed = PARAMS["all_introns_bed6"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/all_introns_counts_and_info.py"

    statement = (
        f"python {script_path} -b {bamfile} "
        f"-g {annotation} "
        f"-o {outfile} "
        f"-vcf {vcffile} "
        f"-bed {utron_bed} -v5"
    )
    P.run(statement, job_memory="16G")


@ruffus.follows(ruffus.mkdir("all_introns_counts_and_info/summarized"))
@ruffus.collate(
    all_introns_counts_and_info,
    ruffus.regex("(.+)/(.+)_.+(hr|4sU).+"),
    r"\1/summarized/\2_summarized.tsv",
)
def summarize_all_introns_counts(infiles, outfile):
    """
    R script to summarize the number of converted counts from each gene in each sample relative to
    unconverted counts, therefore we calcualte the "percent_converted" for each gene in each sample.
    """
    input_folder = os.path.dirname(infiles[0])
    day_regex = os.path.basename(outfile).split("_")[0]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/summarize_counts.R"
    statement = f"Rscript {script_path} {input_folder} {day_regex} {outfile}"

    P.run(statement, job_memory="64G")


@ruffus.follows(
    name_sort_bams,
    feature_count_read_assignments,
    VCF_from_no4sU_control,
    index_VCF,
    all_introns_counts_and_info,
    summarize_all_introns_counts,
)
def full():
    """
    Utility function to run the whole pipeline via:
    python pipeline_slam_3UIs.py make full -v5
    """
    pass


### MISC ###
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
############
