"""
This pipeline takes in mapped RNA seq reads in `.bam`, sorts and 
annotates them, then passes them detects and quantifies the number
of 4sU-linked nucleotide conversions (T>C or A>G), exclusing SNPs.
"""

import sys
import os
from cgatcore import pipeline as P
from ruffus import *


# Read in pipeline.yml
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


@follows(mkdir("name_sorted_bams"))
@transform("input_bams.dir/*.bam",
           regex("input_bams.dir/(.+).star.bam"),
           r"name_sorted_bams/\1.sorted.bam")
def name_sort_bams(infile, outfile):
    """
    Takes `.bam` file and sorts it in name order, therefore reads from the same pair appear sequentially
    """

    statement = "samtools sort -n %(infile)s -o %(outfile)s"
    P.run(statement, job_memory="8G")


@follows(mkdir("read_assignments"))
@transform(name_sort_bams, regex("name_sorted_bams/(.+).sorted.bam"), 
           r"read_assignments/\1/\1.sorted.assigned.bam")
def featureCountsReadAssignments(infile, outfile):
    """
    Takes `.bam` file and adds 2 tags to it: 
        - XT tag tells us the transcript the reads map to
        - XS tag tells us which strand the gene is transcribed from
    """

    annotation = PARAMS["transcript_gtf"]
    outdir = os.path.dirname(outfile)
    rename_outfile = outfile.replace(".assigned.bam", "")
    statement = """featureCounts -a %(annotation)s
                                 -o %(outdir)s/counts.txt
                                 -T 2
                                 -R BAM
                                 %(infile)s && 
                    mv %(rename_outfile)s.bam.featureCounts.bam %(outfile)s"""
    P.run(statement,
          job_threads=2,
          job_memory="8G")


@transform("input_bams.dir/*no4sU*.bam", 
           regex("input_bams.dir/(.+)_no4sU_(.+).bam"), 
           r"snp_vcf/\1.vcf.gz")
def VCF_from_no4sU_control(infile, outfile):
    """
    Takes the negative control `.bam` file (via the no4sU regex) and creates a SNP VCF file
    using Varscan. It looks for differences between the `.bam` file and the genome `.fasta` 
    file to detect differences, where these are always present at the same position (above 
    a treshold, see Varscan docs) it is classed as a SNP and included in output VCF file. 
    """

    outfile_uncompressed = outfile.replace(".vcf.gz", ".vcf")
    genome_fasta = PARAMS["genome_fasta"]
    statement = """samtools mpileup -B -A -f %(genome_fasta)s %(infile)s | 
                    varscan mpileup2snp --variants 1 --output-vcf 1 > %(outfile_uncompressed)s &&
                    bcftools view %(outfile_uncompressed)s -Oz -o %(outfile)s"""
    #removed sort -k1,1 -k2,2n %(outfile_uncompressed)s -o %(outfile_uncompressed)s &&
    #before bcftools view
    P.run(statement, 
          job_threads=12,
          job_memory="1G")


@transform(VCF_from_no4sU_control, 
           suffix(".vcf.gz"), 
           ".vcf.gz.tbi")
def index_VCF(infile, outfile):
    """
    Indexes the VCF file so it can be read with Pysam, generates the .vcf.gz.tbi file
    """
    
    statement = "tabix -p vcf %(infile)s"
    P.run(statement)
    

@follows(index_VCF, mkdir("all_introns_counts_and_info"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
           r"all_introns_counts_and_info/\1.tsv")
def all_introns_counts_and_info(infiles, outfile):
    """
    Takes in the sorted and feature assigned `.bam` file from previous steps and passes them to a
    python script that uses pysam (python wrapper for htslib). This script iterates over the `.bam` 
    file pair-by-pair (representing the sequencing insert), determines whether the read-pair shows
    evidence of intron splicing/retention and assigns these to specific events by referencing the
    `.gtf` and `.bed` files, and XT tag from featureCountsReadAssignments. Next, the script uses the
    XS tag from featureCountsReadAssignment to assign each read in the pair as the forward or reverse
    read, relative to the direction of transcription. Finally, it looks for T>C in FW read, A>G in RV
    read, checks these are not present in the SNP VCF file, and outputs metadata on each read-pair
    about it's event assignment, number of conversions, coverage etc. 
    """
    
    bamfile, vcffile = infiles
    annotation = PARAMS["transcript_gtf"]
    utron_bed = PARAMS["all_introns_bed6"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/all_introns_counts_and_info.py"

    statement = """python %(script_path)s -b %(bamfile)s 
                                            -g %(annotation)s
                                            -o %(outfile)s
                                            -vcf %(vcffile)s 
                                            -bed %(utron_bed)s -v5"""
    
    P.run(statement,
          job_memory="16G")
    

@follows(mkdir("all_introns_counts_and_info/summarized"))
@collate(all_introns_counts_and_info, 
       regex("(.+)/(.+)_.+(hr|4sU).+"),
       r"\1/summarized/\2_summarized.tsv")
def summarize_all_introns_counts(infiles, outfile):
    """
    R script to summarize the number of converted counts from each gene in each sample relative to 
    unconverted counts, therefore we calcualte the "percent_converted" for each gene in each sample.
    """
    
    input_folder = os.path.dirname(infiles[0])
    day_regex = os.path.basename(outfile).split("_")[0]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/summarize_counts.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(day_regex)s %(outfile)s"

    P.run(statement, job_memory="64G")


@follows(name_sort_bams,
         featureCountsReadAssignments, 
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