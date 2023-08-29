import pysam as pysam
from statistics import mean
from statistics import median
from matplotlib import pyplot as plt
import sys
import cgatcore.experiment as E
import pandas as pd
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


    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args.infile_bam)
    #bamfile = pysam.AlignmentFile("../STAR-custom/read_assignments/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3.sorted.assigned.bam")
    
    vcffile = pysam.VariantFile(args.vcf_path)
    #vcffile = pysam.VariantFile('../STAR-custom/snp_vcf/D0.vcf.gz')

    utron_coords = []
    with gzip.open(args.utron_bed, "rt") as bedfile:
        for line in bedfile:
            contents = line.strip().split("\t")
            chromosome, start, end, transcript_id, bedstrand = contents[0], int(contents[1]), int(contents[2]), contents[3], contents[5]
            bed_tuple = (chromosome, start, end, transcript_id, bedstrand)
            utron_coords.append(bed_tuple)


    gtf_cols = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf_df = pd.read_csv(args.gtf_path, sep="\t", comment="#", header=None, names=gtf_cols)
    #gtf_df = pd.read_csv("/shared/sudlab1/General/projects/stem_utrons/existing_hPSC_data/HipSci/HIPSCI-REANNOTATE/filtered_genesets.dir/agg-agg-agg.filtered.gtf.gz", sep="\t", comment="#", header=None, names=gtf_cols)
    #filter the dataframe to only deal with transcripts (as they have stand info)
    gtf_df = gtf_df[gtf_df["feature"]=="transcript"]

    #define a function so we can parse the gene_id/gene_name from the messy attributes column
    def attibute_to_dict(attribute_string):
        attributes = {}
        for attribute in attribute_string.split(";"):
            attribute_not_empty = attribute.split()
            #sometimes attribute is empty therefore splitting it to key, value causes error
            if attribute_not_empty:
                key, value = attribute.split(" ", 1)
                attributes[key] = value.strip("\"")
        return attributes
    
    # convert the attribute column from string to dict
    gtf_df["attribute"] = gtf_df["attribute"].apply(attibute_to_dict)

    # generate a dict of gene_ids to strands
    gene_id = gtf_df["attribute"].apply(lambda x: x.get("gene_id"))
    strands = gtf_df["strand"]
    strand_dict = dict(zip(gene_id, strands))

    conversion_counts = list()
    conversion_dict = {}

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

    status_list_1 = list()
    status_list_2 = list()

    for pair in fragment_iterator(bamfile.fetch(until_eof=True)):
    
        #if i_total_progress >= 100000: 
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

        read1_within_intron = {transcript_id: [(start <= read1_start <= end and 
                                                start <= read1_end <= end and 
                                                chromosome == read1.reference_name), start, end]
                                    for chromosome, start, end, transcript_id, bedstrand in utron_coords
                                    if start <= read1_start <= end and 
                                        start <= read1_end <= end and 
                                        chromosome == read1.reference_name}
        
        read2_within_intron = {transcript_id: [(start <= read2_start <= end and 
                                                start <= read2_end <= end and 
                                                chromosome == read2.reference_name), start, end]
                                    for chromosome, start, end, transcript_id, bedstrand in utron_coords
                                    if start <= read2_start <= end and
                                        start <= read2_end <= end and 
                                        chromosome == read2.reference_name}

        """ if not read1_within_intron:
            pass
            #print("Read1 is not within any 3'UTR introns")
        else:
            print("\n\nRead1 is within at least 1 3'UTR intron")
            print(read1_within_intron)
            print(read1)
        if not read2_within_intron:
            pass
            #print("Read2 is not within any 3'UTR introns")
        else:
            print("\n\nRead2 is within at least 1 3'UTR intron")
            print(read2_within_intron)
            print(read2) """

        read1_spliced_3UI = {transcript_id: [(read1_start < start and 
                                              read1_end > end and 
                                              "N" in read1.cigarstring and 
                                              chromosome == read1.reference_name), start, end]
                             for chromosome, start, end, transcript_id, bedstrand in utron_coords
                             if read1_start < start and 
                                read1_end > end and 
                                "N" in read1.cigarstring and
                                chromosome == read1.reference_name}


        read2_spliced_3UI = {transcript_id: [(read2_start < start and 
                                              read2_end > end and 
                                              "N" in read2.cigarstring and
                                              chromosome == read2.reference_name), start, end]
                             for chromosome, start, end, transcript_id, bedstrand in utron_coords
                             if read2_start < start and 
                                read2_end > end and 
                                "N" in read2.cigarstring and
                                chromosome == read2.reference_name}
        
        """ if not read1_spliced_3UI:
            pass
            #print("Read1 is not a 3UI spliced read")
        else:
            print("\n\nRead1 is a 3UI spliced read")
            print(read1_spliced_3UI)
            print(read1)
        if not read2_spliced_3UI:
            pass
            #print("Read2 is not a 3UI spliced read")
        else:
            print("\n\nRead2 is a 3UI spliced read")
            print(read2_spliced_3UI)
            print(read2) """

        all_dicts = [read1_within_intron, read2_within_intron, read1_spliced_3UI, read2_spliced_3UI]

        all_empty = all(not contents for contents in all_dicts)

        if all_empty:
            continue

        assign_conversions_to_retained = []
        for d in [read1_within_intron, read2_within_intron]:
            assign_conversions_to_retained.extend(list(d.keys()))
        assign_conversions_to_retained = list(set(assign_conversions_to_retained))
        assign_conversions_to_retained = [element + "_retained" for element in assign_conversions_to_retained]

        assign_conversions_to_spliced = []
        for d in [read1_spliced_3UI, read2_spliced_3UI]:
            assign_conversions_to_spliced.extend(list(d.keys()))
        assign_conversions_to_spliced = list(set(assign_conversions_to_spliced))
        assign_conversions_to_spliced = [element + "_spliced" for element in assign_conversions_to_spliced]

        assign_conversions_to = assign_conversions_to_retained + assign_conversions_to_spliced


        i+=1

        read1_status = read1.get_tag("XS")
        read2_status = read2.get_tag("XS")

        ############################
        ### For metadata outputs ###
        status_list_1.append(read1_status)
        status_list_2.append(read2_status)  
        status_list = list()
        status_list.append(read1_status)
        status_list.append(read2_status)
        ############################
    
        if(any(status in ["Assigned", "+", "-"] for status in status_list)):
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
            if(read1.is_reverse):
                reverse_read = read1
                forward_read = read2
            else:
                reverse_read = read2
                forward_read = read1

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
                    for transcript in assign_conversions_to:
                        current_val = conversion_dict.get(transcript, 0)
                        new_val = current_val + 1
                        conversion_dict.update({transcript: new_val})             
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
                    for transcript in assign_conversions_to:
                        current_val = conversion_dict.get(transcript, 0)
                        new_val = current_val + 1
                        conversion_dict.update({transcript: new_val})      
            else:
                # should not be possible - but just in case
                pass 

    with open(args.outfile_txt, "w") as outfile2:
        for key, value in conversion_dict.items():
            outfile2.write(str(key) + "\t" + str(value) + "\n")

    # write footer and output benchmark information.
    E.stop()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
