"""
retrive_3UI_reads.py
====================

This script......



Help:

Run script by 
"""



import pysam as pysam
from statistics import mean
from statistics import median
from matplotlib import pyplot as plt
import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import pandas as pd
import cgat.GTF as GTF
from collections import defaultdict
from collections import Counter
import gzip
import io

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-b", "--bam", dest="infile_bam", type=str,
                        help="Supply a path to the bam file that has undergone read assignment with featureCounts")
    
    parser.add_argument("-g", "--gtf", dest="gtf_path", type=str,
                        help="Supply a path to the transcript assembly gtf file")

    parser.add_argument("-o", "--out", dest="outfile_txt", type=str,
                        help="""Supply a path to the output file. This file will contain 
                        conversions per pair, accounting for stranding""")
    
    parser.add_argument("-vcf", "--vcf", dest="vcf_path", type=str,
                        help="""Supply a path to the VCF.gz file""")    
    
    parser.add_argument("-lfq1", "--labelled-fastq-1", dest="labelled_fastq_1_path", type=str,
                        help="""Supply a path to the gzipped output file. This file will contain
                        the read1s from labelled pairs""")   
    
    parser.add_argument("-lfq2", "--labelled-fastq-2", dest="labelled_fastq_2_path", type=str,
                        help="""Supply a path to the gzipped output file. This file will contain
                        the read2s from labelled pairs""") 

    parser.add_argument("-ulfq1", "--unlabelled-fastq-1", dest="unlabelled_fastq_1_path", type=str,
                        help="""Supply a path to the gzipped output file. This file will contain
                        the read1s from unlabelled pairs""") 
    
    parser.add_argument("-ulfq2", "--unlabelled-fastq-2", dest="unlabelled_fastq_2_path", type=str,
                        help="""Supply a path to the gzipped output file. This file will contain
                        the read2s from unlabelled pairs""")       


    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args.infile_bam)
    #bamfile = pysam.AlignmentFile("../STAR-custom/read_assignments/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3.sorted.assigned.bam")
    
    vcffile = pysam.VariantFile(args.vcf_path)
    #vcffile = pysam.VariantFile('../STAR-custom/snp_vcf/D0.vcf.gz')

    tx2gene = defaultdict(list)
    strand_dict = defaultdict(str)

    for entry in GTF.iterator(iotools.open_file(args.gtf_path)):
        if not entry.feature == "transcript":
            continue
        strand_dict[entry.gene_id] = entry.strand        
        tx2gene[entry.gene_id].append(entry.transcript_id)

    def fragment_iterator(read_iterator):

        read_list = list()
        last_read = None

        for read in read_iterator:
            if last_read is not None and last_read != read.query_name:
                yield read_list
                read_list = list()
                last_read = read.query_name
            last_read = read.query_name
            read_list.append(read)

        yield read_list
    
    i_progress = 0
    i_total_progress = 0

    with gzip.open(args.labelled_fastq_1_path, "wt") as labelled_fastq_1, \
          gzip.open(args.labelled_fastq_2_path, "wt") as labelled_fastq_2, \
          gzip.open(args.unlabelled_fastq_1_path, "wt") as unlabelled_fastq_1, \
          gzip.open(args.unlabelled_fastq_2_path, "wt") as unlabelled_fastq_2:
    
        for pair in fragment_iterator(bamfile.fetch(until_eof=True)):
        
            #if i_total_progress >= 200000: 
            #    break

            if len(pair)!=2:
                continue

            read1, read2 = pair

            if read1.is_unmapped or read2.is_unmapped:
                continue 

            i_total_progress+=1
            i_progress+=1

            if i_progress == 10000:
                E.debug(str(i_total_progress) + " pairs processed")
                i_progress = 0


            read1_status = read1.get_tag("XS")
            read2_status = read2.get_tag("XS")

            status_list = list()
            status_list.append(read1_status)
            status_list.append(read2_status)
        
            if not (any(status in ["Assigned", "+", "-"] for status in status_list)):
                continue
            
            # pass if either is assigned
            if(all(status in ["Assigned", "+", "-"] for status in status_list)):
                # pass if both are assigned
                assignment1 = read1.get_tag("XT")
                assignment2 = read2.get_tag("XT")
                strand1 = strand_dict[assignment1]
                strand2 = strand_dict[assignment2]
                if(strand1 == strand2):
                    # pass if both on same strand
                    strand = strand1
                else:
                    # if not on same strand bin the pair
                    continue
            else:
                # pass if only 1 is assigned
                if(read1_status=="Assigned"):
                    # if read 1 is the 1 assigned
                    assignment = read1.get_tag("XT")
                    strand = strand_dict[assignment]
                else:
                    # if read 2 is the 1 assigned]
                    assignment = read2.get_tag("XT")
                    strand = strand_dict[assignment]

            # assigned a "forward" and "reverse" read relative to the genome
            if(read1.is_reverse and not read2.is_reverse):
                reverse_read = read1
                forward_read = read2
            elif(read2.is_reverse and not read1.is_reverse):
                reverse_read = read2
                forward_read = read1
            else:
                #Not proper pair
                continue

            # if we are mapped to a +ve stranded transcript, then count T>C in 
            # the forward read and A>G in the reverse read. if we are mapped to
            # a -ve stranded transcript, count T>C in the reverse read and A>G
            # in the forward read.
            if strand == "+":
                # pass if mapped to +ve transcript
                conversions = 0
                for base in forward_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if(None in base):
                        continue

                    read_seq = forward_read.query_sequence[read_pos]
                                        
                    if read_seq == "C" and genome_seq == "t":
                        variants_at_position = list(vcffile.fetch(forward_read.reference_name, genome_pos, genome_pos+1)) 
                        if variants_at_position:                            
                            if(any(variant_at_pos.alts[0]=="C" for variant_at_pos in variants_at_position)):
                                pass
                            else:
                                conversions += 1                  
                        else:
                            conversions += 1 

                for base in reverse_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if(None in base):
                        continue

                    read_seq = reverse_read.query_sequence[read_pos]

                    if read_seq == "G" and genome_seq == "a":
                        variants_at_position = list(vcffile.fetch(reverse_read.reference_name, genome_pos, genome_pos+1)) 
                        if variants_at_position:
                            if(any(variant_at_pos.alts[0]=="G" for variant_at_pos in variants_at_position)):
                                pass
                            else:
                                conversions += 1                          

                        else:
                            conversions += 1 
                
                if(conversions >= 1): #if there is a conversion, pulldown
                    ### ADD TO LABELLED FASTQ
                    labelled_fastq_1.write(f"@{read1.query_name}\n{read1.query_sequence}\n+\n{read1.qual}\n")
                    labelled_fastq_2.write(f"@{read2.query_name}\n{read2.query_sequence}\n+\n{read2.qual}\n")
                else:
                    ### ADD TO UNLABELLED FASTQ
                    unlabelled_fastq_1.write(f"@{read1.query_name}\n{read1.query_sequence}\n+\n{read1.qual}\n")
                    unlabelled_fastq_2.write(f"@{read2.query_name}\n{read2.query_sequence}\n+\n{read2.qual}\n")
                                
            elif strand == "-":
                # pass if mapped to -ve transcript
                conversions = 0
                for base in forward_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if(None in base):
                        continue

                    read_seq = forward_read.query_sequence[read_pos]

                    if read_seq == "G" and genome_seq == "a":
                        variants_at_position = list(vcffile.fetch(forward_read.reference_name, genome_pos, genome_pos+1)) 
                        if variants_at_position:
                            if(any(variant_at_pos.alts[0]=="G" for variant_at_pos in variants_at_position)):
                                pass
                            else:
                                conversions += 1                           

                        else:
                            conversions += 1 

                for base in reverse_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if(None in base):
                        continue

                    read_seq = reverse_read.query_sequence[read_pos]

                    if read_seq == "C" and genome_seq == "t":
                        variants_at_position = list(vcffile.fetch(reverse_read.reference_name, genome_pos, genome_pos+1)) 
                        if variants_at_position:
                            if(any(variant_at_pos.alts[0]=="C" for variant_at_pos in variants_at_position)):
                                pass
                            else:
                                conversions += 1                           

                        else:
                            conversions += 1 
                if(conversions >= 1): #if there is a conversion, pulldown
                    ### ADD TO LABELLED FASTQ
                    labelled_fastq_1.write(f"@{read1.query_name}\n{read1.query_sequence}\n+\n{read1.qual}\n")
                    labelled_fastq_2.write(f"@{read2.query_name}\n{read2.query_sequence}\n+\n{read2.qual}\n")
                else:
                    ### ADD TO UNLABELLED FASTQ
                    unlabelled_fastq_1.write(f"@{read1.query_name}\n{read1.query_sequence}\n+\n{read1.qual}\n")
                    unlabelled_fastq_2.write(f"@{read2.query_name}\n{read2.query_sequence}\n+\n{read2.qual}\n")
            else:
                # should not be possible - but just in case
                pass 

    # write footer and output benchmark information.
    E.stop()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
