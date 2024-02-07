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
    P.run(statement,
            job_memory="16G")


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
          job_memory="8G")


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
          job_memory="32G")
    
    
@follows(mkdir("stranded_conversions_per_pair/plots"))
@transform(stranded_conversions_per_pair, regex("(.+)/(.+).txt"), r"\1/plots/\2.png")
def plot_stranded_conversions_per_pair(infile, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/plot_conversions_per_pair.R"
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"
    P.run(statement,
            job_memory="8G")
   

@transform("input_bams.dir/*no4sU*.bam", regex("input_bams.dir/(.+)_no4sU_(.+).bam"), r"snp_vcf/\1.vcf.gz")
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
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
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
          job_memory="64G")


@follows(mkdir("stranded_conversions_per_pair_no_snp/plots"))
@transform(stranded_conversions_per_pair_no_snp, regex("(.+)/(.+).txt"), r"\1/plots/\2.png")
def plot_stranded_conversions_per_pair_no_snp(infile, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/plot_conversions_per_pair.R"
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"
    P.run(statement,
            job_memory="8G")


@follows(index_VCF, mkdir("pulled_3UI_reads"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
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
          job_memory="64G")
    
@follows(index_VCF, mkdir("3UI_counts_and_info"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
           r"3UI_counts_and_info/\1.tsv")
def spliced_counts_and_info(infiles, outfile):
    bamfile, vcffile = infiles
    annotation = PARAMS["transcript_gtf"]
    utron_bed = PARAMS["utrons_bed6"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/3UI_spliced_counts_and_info.py"

    statement = """python %(script_path)s -b %(bamfile)s 
                                            -g %(annotation)s
                                            -o %(outfile)s
                                            -vcf %(vcffile)s 
                                            -bed %(utrons_bedfile)s -v5"""
    
    P.run(statement,
          job_memory="16G")
    
@follows(index_VCF, mkdir("split_labelled_vs_unlabelled"))
@subdivide(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
           [r"split_labelled_vs_unlabelled/\1/\1.labelled.fastq.1.gz",
            r"split_labelled_vs_unlabelled/\1/\1.labelled.fastq.2.gz",
            r"split_labelled_vs_unlabelled/\1/\1.unlabelled.fastq.1.gz",
            r"split_labelled_vs_unlabelled/\1/\1.unlabelled.fastq.2.gz"])
def split_labelled_vs_unlabelled(infiles, outfiles):
    bamfile, vcffile = infiles
    lfq1, lfq2, ulfq1, ulfq2 = outfiles
    annotation = PARAMS["transcript_gtf"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/split_labelled_vs_unlabelled.py"

    statement = """python %(script_path)s -b %(bamfile)s 
                                            -g %(annotation)s
                                            -vcf %(vcffile)s
                                            -lfq1 %(lfq1)s
                                            -lfq2 %(lfq2)s
                                            -ulfq1 %(ulfq1)s
                                            -ulfq2 %(ulfq2)s -v5"""
    
    P.run(statement,
          job_memory="64G")
    
@follows(mkdir("salmon_index"), mkdir("salmon_index/DAT"))
@transform(PARAMS["transcript_gtf"],
           regex(".+/(.+).gtf.gz"),
           r"salmon_index/\1.salmon.index")
def makeSalmonIndex(infile,outfile):
    # Long transcripts cause indexing to use lots of memory?
    job_memory=PARAMS["salmon_index_memory"]
    job_threads=PARAMS["salmon_index_threads"]

    gtf_basename = P.snip(os.path.basename(infile), ".gtf.gz")
    transcript_fasta = "salmon_index/" + gtf_basename + "transcripts.fa"
    fastaref = PARAMS["genome_fasta"]
    index_options=PARAMS["salmon_indexoptions"]
    tmpfile = P.get_temp_filename()
    
    #statement since salmon >v1.0 Selective Alignment update
    #statement now generates decoy.txt from reference genome, and concats the genome.fa and transcriptome.fa into gentrome.fa
    #gentrome.fa is then passed through salmon alongside the decoys to create decoy-aware transcriptome (DAT)
    statement = '''
    gunzip -c %(infile)s > %(tmpfile)s &&
    gffread %(tmpfile)s -g %(fastaref)s -w %(transcript_fasta)s &&
    grep "^>" <%(fastaref)s | cut -d " " -f 1 > salmon_index/DAT/decoys.txt &&
    sed -i.bak -e 's/>//g' salmon_index/DAT/decoys.txt &&
    cat %(transcript_fasta)s %(fastaref)s > salmon_index/DAT/gentrome.fa.gz &&
    salmon index
      -p %(job_threads)s
      %(index_options)s
      -t salmon_index/DAT/gentrome.fa.gz
      -d salmon_index/DAT/decoys.txt
      -i %(outfile)s &&
    rm %(tmpfile)s
    '''
    P.run(statement)


@follows(mkdir("quant_labelled_vs_unlabelled"))
@collate(split_labelled_vs_unlabelled,
           regex(".+/(.+)\.(.+)\.fastq.(.+).gz"),
           add_inputs(makeSalmonIndex),
           r"quant_labelled_vs_unlabelled/\1/\2/quant.sf")
def quant_labelled_vs_unlabelled(infiles, outfile):
    '''Quantify labelled vs unlabelled reads'''
    job_threads=PARAMS["salmon_threads"]
    job_memory=PARAMS["salmon_memory"]
    #print("infiles=")
    #print(infiles)
    input1, input2 = infiles
    #print("input1 = ")
    #print(input1)
    #print("input2 = ")
    #print(input2)
    fastq1, salmon_index = input1
    #print("fastq1 =")
    #print(fastq1)
    #print("salmon_index = ")
    #print(salmon_index)
    fastq2, salmon_index = input2
    #print("fastq2 = ")
    #print(fastq2)
    #print("salmon index = ")
    #print(salmon_index)

    job_options = PARAMS["cluster_options"] 
    outdir = os.path.dirname(outfile)
    salmon_options=PARAMS["salmon_quantoptions"]

    statement = '''
    salmon quant -i %(salmon_index)s
        --libType IU
        -1 %(fastq1)s
        -2 %(fastq2)s
        -o %(outdir)s
        -p %(job_threads)s
        %(salmon_options)s
    '''

    P.run(statement)

@follows(plot_conversions_per_read, 
         plot_conversions_per_pair, 
         featureCountsReadAssignments, 
         stranded_conversions_per_pair,
         plot_stranded_conversions_per_pair,
         VCF_from_no4sU_control,
         index_VCF,
         stranded_conversions_per_pair_no_snp,
         plot_stranded_conversions_per_pair_no_snp,
         pull_3UI_reads,
         split_labelled_vs_unlabelled,
         makeSalmonIndex,
         quant_labelled_vs_unlabelled,
         spliced_counts_and_info)

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
