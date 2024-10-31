"""
This pipeline takes in pre-mapped alignments and puts them through slam dunk
test to check how branch switching works
"""

from heapq import merge
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
                                            -bed %(utron_bed)s -v5"""
    
    P.run(statement,
          job_memory="16G")
    
### FOR ALL SPLICED EVENTS ###################################################################################################

@follows(index_VCF, mkdir("all_introns_counts_and_info"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
           r"all_introns_counts_and_info/\1.tsv")
def all_introns_counts_and_info(infiles, outfile):
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
    input_folder = os.path.dirname(infiles[0])
    day_regex = os.path.basename(outfile).split("_")[0]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/summarize_counts.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(day_regex)s %(outfile)s"

    P.run(statement, job_memory="64G")

@follows(mkdir("all_introns_counts_and_info/pairs"))
@transform(summarize_all_introns_counts,
         regex("(.+)/summarized/(.+)_summarized.tsv"),
       r"\1/pairs/\2_pvalues.tsv")
def get_pair_pvalues_all_events(infiles, outfile):
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/get_pair_pvalues_all_events.R"
    infile = str(infiles)
    statement = "Rscript %(script_path)s %(infile)s %(outfile)s"

    P.run(statement, job_memory="64G")

#################################################################################################################################

@follows(index_VCF, mkdir("gene_level_counts_and_info"))
@transform(featureCountsReadAssignments, 
           regex("read_assignments/(.+)/(.+)_(.+)_(.+)_(.+)_(.+).sorted.assigned.bam"), 
           add_inputs(r"snp_vcf/\2.vcf.gz"),
           r"gene_level_counts_and_info/\1.tsv")
def gene_level_counts_and_info(infiles, outfile):
    bamfile, vcffile = infiles
    annotation = PARAMS["transcript_gtf"]
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/gene_level_counts_and_info.py"

    statement = """python %(script_path)s -b %(bamfile)s 
                                            -g %(annotation)s
                                            -o %(outfile)s
                                            -vcf %(vcffile)s -v5"""
    
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

@follows(mkdir("3UI_counts_and_info/pairs"))
@collate(spliced_counts_and_info,
         regex("(.+)/(.+)_.+(hr|4sU).+"),
       r"\1/pairs/\2_pvalues.tsv")
def get_pair_pvalues(infiles, outfile):
    day_regex = os.path.basename(outfile).split("_")[0]
    input_folder = os.path.dirname(infiles[0])
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/get_pair_pvalues.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(day_regex)s %(outfile)s"

    P.run(statement, job_memory="16G")

### Run half lives after running pvals, then we can omit those where the initial 
### half life estimate is negative or >24hr. 

@follows(get_pair_pvalues)
@collate(spliced_counts_and_info, 
       regex("(.+)/(.+)_.+(hr|4sU).+"),
       r"\1/pairs/\2_half_lives.tsv")
def get_pair_half_lives(infiles, outfile):
    day_regex = os.path.basename(outfile).split("_")[0]
    input_folder = os.path.dirname(infiles[0])
    num_bootstraps = 100
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/get_pair_halflives.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(day_regex)s %(num_bootstraps)s %(outfile)s"

    P.run(statement, job_memory="16G")

@follows(mkdir("3UI_counts_and_info/over_time"))
@collate(spliced_counts_and_info, 
       regex("(.+)/(.+)_.+(hr|4sU).+"),
       r"\1/over_time/over_time_pvalues.tsv")
def get_over_time_pvalues(infiles, outfile):
    input_folder = os.path.dirname(infiles[0])
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/get_over_time_pvalues.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(outfile)s"

    P.run(statement, job_memory="32G")

@follows(mkdir("3UI_counts_and_info/over_time"))
@collate(spliced_counts_and_info, 
       regex("(.+)/(.+)_.+(hr|4sU).+"),
       r"\1/over_time/interaction_over_time_pvalues.tsv")
def get_interaction_over_time_pvals(infiles, outfile):
    input_folder = os.path.dirname(infiles[0])
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/get_interaction_over_time_pvals.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(outfile)s"

    P.run(statement, job_memory="32G")

@follows(mkdir("gene_level_counts_and_info/over_time"))
@collate(gene_level_counts_and_info, 
       regex("(.+)/(.+)_.+(hr|4sU).+"),
       r"\1/over_time/over_time_pvalues.tsv")
def get_over_time_pvalues_gene_level(infiles, outfile):
    input_folder = os.path.dirname(infiles[0])
    script_path = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/get_over_time_pvalues_gene_level.R"
    statement = "Rscript %(script_path)s %(input_folder)s %(outfile)s"

    P.run(statement, job_memory="64G")

@follows(mkdir("dapars"))
@transform(PARAMS["transcript_gtf"],
            regex("(.+)/(.+).gtf.gz"),
            r"dapars/\2.bed12")
def gene_model_bed12(infile, outfile):
    statement = """cgat gff2bed 
                    --bed12-from-transcripts
                    -I %(infile)s
                    -S %(outfile)s
                    -L %(outfile)s.log"""
    P.run(statement, job_memory="8G")

@transform(gene_model_bed12, regex("(.+)/(.+).bed12"), r"\1/\2.sorted.bed12")
def sort_bed12(infile, outfile):
    statement = "sort -k 1,1 %(infile)s > %(outfile)s"
    P.run(statement)

@follows(mkdir("dapars/sorted_bams"))
@transform("input_bams.dir/*.bam", regex("(.+)/(.+).bam"), r"dapars/sorted_bams/\2.sorted.bam")
def pos_sort_bam(infile, outfile):
    statement = "samtools sort %(infile)s -o %(outfile)s"
    P.run(statement, job_memory="8G")

@follows(mkdir("dapars/bedGraphs"))
@transform(pos_sort_bam, regex("(.+)/(.+)/(.+).sorted.bam"), r"dapars/bedGraphs/\3.bedGraph")
def get_bedGraphs(infile, outfile):
    genome_file = PARAMS["genome_coords"]
    statement = "genomeCoverageBed -split -bg -ibam %(infile)s -g %(genome_file)s > %(outfile)s"
    P.run(statement, job_memory="8G")

@transform(PARAMS["transcript_gtf"],
           regex("(.+)/(.+).gtf.gz"),
           r"dapars/\2.tx2gene")
def get_dapars_tx2gene(infile, outfile):
    outlines = []
    for entry in GTF.iterator(IOTools.open_file(infile)):
        if not entry.feature == "transcript":
            continue  

        transcript_id = entry.transcript_id
        if(len(entry.attributes.split(";")) == 5 ):
            mstrg_gene = entry.attributes.split(";")[0].strip("gene_id \"").strip("\"").strip()
            ref_gene_id = entry.attributes.split(";")[3].strip("ref_gene_id \"").strip("\"").strip()

        outlines.append((transcript_id, ref_gene_id))

    outlines = list(set(outlines))

    with open(outfile, "wt") as out:
        for line in outlines:
            tx, gene = line
            out.write(f"{tx}\t{gene}\n")

@transform(gene_model_bed12, 
           regex("(.+)/(.+).bed12"),
           add_inputs(get_dapars_tx2gene),
           r"dapars/\2.dapars.anno")
def dapars_extract_anno(infiles, outfile):
    bed12, tx2gene = infiles
    dapars_scripts = "/users/mbp20jjr/dapars/src"
    statement = """python %(dapars_scripts)s/DaPars_Extract_Anno.py
                     -b %(bed12)s
                     -s %(tx2gene)s
                     -o %(outfile)s"""
    P.run(statement)

@follows(mkdir("dapars/output"))
@subdivide("dapars_design.tsv",
           formatter(),
           add_inputs(get_bedGraphs, dapars_extract_anno),
           ["dapars/output/%s_dapars_config.txt" % line.split()[0]
            for line in IOTools.open_file("dapars_design.tsv")])
def generate_dapars_config(infiles, outfiles):
    bedGraphs = infiles[1:-1]
    anno = infiles[-1]
    config_base = os.path.dirname(os.path.abspath(__file__)) + "/pipeline_slam_3UIs/dapars_default_config.txt"

    comparisons = [line.split() for line in IOTools.open_file("dapars_design.tsv")]

    import fnmatch

    for name, condition1, condition2 in comparisons:
        condition1  = "*" + condition1 + "*"
        condition2 = "*" + condition2 + "*"
        condition1_files = [f for f in bedGraphs if fnmatch.fnmatch(f, condition1)]
        condition1_filestring = ','.join(condition1_files)
        condition2_files = [f for f in bedGraphs if fnmatch.fnmatch(f, condition2)]
        condition2_filestring = ','.join(condition2_files)
        name=name
        name2 = "*" + name + "*"
        outfile = [o for o in outfiles if fnmatch.fnmatch(o, name2)]
        with open(outfile[0], "wt") as out:
            out.write(f"Annotated_3UTR={anno}\n")
            out.write(f"Group1_Tophat_aligned_Wig={condition1_filestring}\n")
            out.write(f"Group2_Tophat_aligned_Wig={condition2_filestring}\n")
            out.write(f"Output_directory=dapars/output/{name}\n")
            out.write(f"Output_result_file={name}.dapars.output\n")

            default_config = IOTools.open_file(config_base)
            for line in default_config:
                out.write(line)

@transform(generate_dapars_config,
           regex("dapars/output/(.+)_dapars_config.txt"),
           r"dapars/output/\1/\1.dapars.output")
def run_dapars_main(infile, outfile):
    dapars_scripts = "/users/mbp20jjr/dapars/src"
    statement = "python %(dapars_scripts)s/DaPars_main.py %(infile)s"
    P.run(statement, job_memory="32G")


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
         spliced_counts_and_info,
         gene_level_counts_and_info,
         all_introns_counts_and_info,
         get_pair_pvalues,
         get_pair_half_lives,
         get_over_time_pvalues,
         get_over_time_pvalues_gene_level,
         get_interaction_over_time_pvals,
         summarize_all_introns_counts,
         get_pair_pvalues_all_events
         )

def full():
    pass

@follows(gene_model_bed12,
         sort_bed12,
         pos_sort_bam,
         get_bedGraphs,
         get_dapars_tx2gene,
         dapars_extract_anno,
         run_dapars_main)
def run_dapars():
    pass

### MISC ###

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

############
