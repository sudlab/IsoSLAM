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
    
    parser.add_argument("-bed", dest="utron_bed", type=str,
                        help="Supply a path to the utron bed file")

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

    utron_coords = []
    
    
    utron_coords = defaultdict(list)

    #bed_path = "/shared/sudlab1/General/projects/stem_utrons/existing_hPSC_data/HipSci/HIPSCI-REANNOTATE/utron_beds.dir/agg-agg-agg.all_utrons.bed.gz"
    
    with iotools.open_file(args.utron_bed) as bedfile:
        for line in bedfile:
            contents = line.strip().split("\t")
            chromosome, start, end, transcript_id, bedstrand = contents[0], int(contents[1]), int(contents[2]), contents[3], contents[5]
            bed_tuple = (chromosome, start, end, transcript_id, bedstrand)
            utron_coords[transcript_id].append(bed_tuple)
            #utron_coords.append(bed_tuple)

    print(utron_coords)

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
    
    i = 0
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
                E.debug(str(i) + " 3UI spliced/retained pairs proccessed")
                i_progress = 0

            read1_start = read1.reference_start
            #print("read1 start: " + str(read1_start))
            read1_end = read1.reference_end
            #print("read1 end: " + str(read1_end))
            read2_start = read2.reference_start
            #print("read2 start: " + str(read2_start))
            read2_end = read2.reference_end
            #print("read2 end: " + str(read2_end))

            if read1.get_tag("XS") == "Assigned":
                transcripts_read1 = tx2gene[read1.get_tag("XT")]
                utrons1 = [utron_coords[tx] for tx in transcripts_read1]
                utrons1 = sum(utrons1, [])
            else:
                utrons1 = list()
                
            if read2.get_tag("XS") == "Assigned":
                transcripts_read2 = tx2gene[read2.get_tag("XT")]
                utrons2 = [utron_coords[tx] for tx in transcripts_read2]
                utrons2 = sum(utrons2, [])
            else:
                utrons2 = list()

            # Dictionary, keyed on transcript id, of introns coords where the read overlaps
            # the intron.
            read1_within_intron = set(transcript_id + "_retained"
                                    for _, start, end, transcript_id, _ in utrons1
                                    if (start <= read1_start <= end or 
                                        start <= read1_end <= end))
            
            read2_within_intron = set(transcript_id + "_retained"
                                    for chromosome, start, end, transcript_id, bedstrand in utrons2
                                    if start <= read2_start <= end or
                                        start <= read2_end <= end)

            # Dictionary, keyed on transcript id of introns where the read splices at identical
            # locations to a utron. 
            block_starts, block_ends = zip(*read1.get_blocks())
            read1_spliced_3UI = {transcript_id: (start, end)
                                for chromosome, start, end, transcript_id, bedstrand in utrons1
                                if start in block_ends and 
                                    end in block_starts}

            block_starts, block_ends = zip(*read2.get_blocks())
            read2_spliced_3UI = {transcript_id: (start, end)
                                for chromosome, start, end, transcript_id, bedstrand in utrons2
                                if start in block_ends and 
                                    end in block_starts}

            # Check one read in the pair has caused the creation of at least one entry.
            all_dicts = [read1_within_intron,
                        read2_within_intron,
                        read1_spliced_3UI,
                        read2_spliced_3UI]
            all_empty = all(not contents for contents in all_dicts)

            if all_empty:
                continue

            i+=1

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
