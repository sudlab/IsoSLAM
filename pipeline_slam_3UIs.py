"""
This pipeline takes in pre-mapped alignments and puts them through slam dunk
"""

import sys
import os
import sqlite3
import csv
import re
import glob

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

#@follows(mkdir("map"))
#@collate("input_fastq.dir/*", regex(".+/(.+).fastq.(.+).gz"), r"map/\1.bam")
#def map_paired_reads(infiles, outfile):
#    read1 = infiles[0]
#    read2 = infiles[1]
#    genome_fasta = PARAMS["genome_fasta"]
#    star_genome_dir = PARAMS["star_genome_dir"]
#    star_options = PARAMS["star_options"]
#    n_threads = PARAMS["number_of_threads"]
#    job_memory = PARAMS["memory"]
#    statement = """STAR --runMode alignReads
#                        --runThreadN %(number_of_threads)s
#                        --genomeDir %(star_genome_dir)s"""
#    P.run(statement,
#         job_memory= job_memory,
#         job_threads=n_threads)

#@follows(mkdir("filter"))
#@transform(map_paired_reads,
#           regex("map/(.+).bam"),
#           r"filter/\1_filtered.bam")
#def slamdunk_filter(infile, outfile):
#    n_threads = 1
#    out_dir = os.path.dirname(outfile)
#    job_memory="8G"
#    utr_bed = PARAMS["utr_bed"]
#    statement = """slamdunk filter -o %(out_dir)s -nm 60 -t %(n_threads)s -b %(utr_bed)s %(infile)s"""
#    P.run(statement,
#          job_memory="8G",
#          job_threads=n_threads)

#@follows(mkdir("rates"))
#@transform(slamdunk_filter, regex("filter/(.+)_slamdunk_mapped_filtered.bam"), r"rates/\1/\1_slamdunk_mapped_filtered_overallrates.csv")
#def conversion_rates(infile, outfile):
#    genome_fasta = PARAMS["reference_genome"]
#    out_dir = os.path.dirname(outfile)
#    job_memory = "64G"
#    n_threads = 1
#    statement = "alleyoop rates -o %(out_dir)s -r %(genome_fasta)s -t %(n_threads)s %(infile)s"
#    P.run(statement,
#          job_memory=job_memory,
#          job_threads=n_threads)



@follows(mkdir("conversions_per_read"))
@transform("input_bams.dir/*.bam",
           regex("input_bams.dir/(.+).star.bam"),
           r"conversions_per_read/\1.txt" )
def conversions_per_read(infile, outfile):
    '''Place docstring here'''
    script_path = os.path.dirname(os.path.abspath(__file__)) +\
          "/pipeline_slam_3UIs/conversions_per_read.py"
    statement = "python %(script_path)s -b %(infile)s -o %(outfile)s -v5"
    P.run(statement)


@follows(mkdir("conversions_per_read/plots"))
@transform(conversions_per_read,
            regex("(.+)/(.+).txt"), 
            r"\1/plots/\2.png")
def plot_conversions_per_read(infile, outfile):
    '''Place docstring here'''
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "pipeline_slam_3UIs/plot_conversions_per_read.R")
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"
    P.run(statement,
          job_memory="8G")


@follows(mkdir("name_sorted_bams"))
@transform("input_bams.dir/*.bam",
           regex("input_bams.dir/(.+).star.bam"),
           r"name_sorted_bams/\1.sorted.bam")
def name_sort_bams(infile, outfile):
    statement = "samtools sort -n %(infile)s -o %(outfile)s"
    P.run(statement, job_memory="8G")


@follows(mkdir("conversions_per_pair"))
@transform(name_sort_bams, regex("name_sorted_bams/(.+).sorted.bam"), r"conversions_per_pair/\1.txt")
def conversions_per_pair(infile, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/conversions_per_pair.py"
    statement = "python %(script_path)s -b %(infile)s -o %(outfile)s -v5"
    P.run(statement)


@follows(mkdir("conversions_per_pair/plots"))
@transform(conversions_per_pair, regex("(.+)/(.+).txt"), r"\1/plots/\2.png")
def plot_conversions_per_pair(infile, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/plot_conversions_per_pair.R"
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"
    P.run(statement,
            job_memory="8G")


@follows(mkdir("read_assignments"))
@transform(name_sort_bams, regex("name_sorted_bams/(.+).sorted.bam"), r"read_assignments/\1/\1.sorted.assigned.bam")
def featureCountsReadAssignments(infile, outfile):
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
          job_memory="4G")


@follows(mkdir("stranded_conversions_per_pair"))
@transform(featureCountsReadAssignments, regex("read_assignments/(.+)/(.+).sorted.assigned.bam"), r"stranded_conversions_per_pair/\1.txt")
def stranded_conversions_per_pair(infile, outfile):
    annotation = PARAMS["transcript_gtf"]
    outdir_subdir = outfile.replace(".txt", "")
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/stranded_conversions_per_pair.py"

    statement = """mkdir -p %(outdir_subdir)s &&
                    python %(script_path)s -b %(infile)s 
                                            -g %(annotation)s
                                            -o %(outfile)s
                                            -om %(outdir_subdir)s -v5"""
    
    P.run(statement,
          job_memory="8G")
    
    
@follows(mkdir("stranded_conversions_per_pair/plots"))
@transform(stranded_conversions_per_pair, regex("(.+)/(.+).txt"), r"\1/plots/\2.png")
def plot_stranded_conversions_per_pair(infile, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/plot_conversions_per_pair.R"
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"
    P.run(statement,
            job_memory="8G")
   

@transform("input_bams.dir/*0uM*.bam", regex("input_bams.dir/(.+)_0uM_(.+).bam"), r"snp_vcf/\1.vcf.gz")
def VCF_from_no4sU_control(infile, outfile):

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


@transform(VCF_from_no4sU_control, suffix(".vcf.gz"), ".vcf.gz.tbi")
def index_VCF(infile, outfile):
    statement = "tabix -p vcf %(infile)s"
    P.run(statement)


@follows(index_VCF, mkdir("stranded_conversions_per_pair_no_snp"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+).sorted.assigned.bam"), 
           add_inputs(VCF_from_no4sU_control),
           r"stranded_conversions_per_pair_no_snp/\1.txt")
def stranded_conversions_per_pair_no_snp(infiles, outfile):
    bamfile, vcffile = infiles
    annotation = PARAMS["transcript_gtf"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/stranded_conversions_per_pair_no_snp.py"

    statement = """python %(script_path)s -b %(bamfile)s 
                                            -g %(annotation)s
                                            -o %(outfile)s
                                            -vcf %(vcffile)s -v5"""
    
    P.run(statement,
          job_memory="8G")


@follows(mkdir("stranded_conversions_per_pair_no_snp/plots"))
@transform(stranded_conversions_per_pair_no_snp, regex("(.+)/(.+).txt"), r"\1/plots/\2.png")
def plot_stranded_conversions_per_pair_no_snp(infile, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/plot_conversions_per_pair.R"
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"
    P.run(statement,
            job_memory="8G")


@follows(index_VCF, mkdir("pulled_3UI_reads"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+).sorted.assigned.bam"), 
           add_inputs(VCF_from_no4sU_control),
           r"pulled_3UI_reads/\1.txt")
def pull_3UI_reads(infiles, outfile):
    bamfile, vcffile = infiles
    annotation = PARAMS["transcript_gtf"]
    utron_bed = PARAMS["utrons_bedfile"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/retrieve_3UI_reads.py"

    statement = """python %(script_path)s -b %(bamfile)s 
                                            -g %(annotation)s
                                            -o %(outfile)s
                                            -vcf %(vcffile)s 
                                            -bed %(utrons_bedfile)s -v5"""
    
    P.run(statement,
          job_memory="128G")

@follows(plot_conversions_per_read, 
         plot_conversions_per_pair, 
         featureCountsReadAssignments, 
         stranded_conversions_per_pair,
         plot_stranded_conversions_per_pair,
         VCF_from_no4sU_control,
         index_VCF,
         stranded_conversions_per_pair_no_snp,
         plot_stranded_conversions_per_pair_no_snp,
         pull_3UI_reads)

def full():
    pass

### MISC ###

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

############
