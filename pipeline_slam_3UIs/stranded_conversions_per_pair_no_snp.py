import pysam as pysam
from statistics import mean
from statistics import median
from matplotlib import pyplot as plt
import sys
import cgatcore.experiment as E
import pandas as pd
from collections import Counter

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


    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args.infile_bam)
    #bamfile = pysam.AlignmentFile("../STAR-custom/read_assignments/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3.sorted.assigned.bam")
    
    #vcffile = pysam.VariantFile('../STAR-custom/snp_vcf/D0.vcf.gz')
    vcffile = pysam.VariantFile(args.vcf_path)

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

    status_list_1 = list()
    status_list_2 = list()

    for pair in fragment_iterator(bamfile.fetch(until_eof=True)):
    
        #if i >= 10000: 
        #    break

        if len(pair)!=2:
            continue

        read1, read2 = pair

        if read1.is_unmapped or read2.is_unmapped:
            continue 

        i+=1
        i_progress+=1

        if i_progress == 100000:
            E.debug(str(i) + " pairs processed")
            i_progress = 0

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
                print("both assigned")
                print(strand1)
                print(strand2)
                if(strand1 == strand2):
                    # pass if both on same strand
                    print("both same")
                    strand = strand1
                else:
                    # if not on same strand bin the pair
                    print("different, bin")
                    continue
            else:
                # pass if only 1 is assigned
                print("1 assigned")
                if(read1_status=="Assigned"):
                    # if read 1 is the 1 assigned
                    assignment = read1.get_tag("XT")
                    strand = strand_dict[assignment]
                    print(strand)
                else:
                    # if read 2 is the 1 assigned]
                    assignment = read2.get_tag("XT")
                    strand = strand_dict[assignment]
                    print(strand)

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
                            print(forward_read)
                            print(f"Read: {forward_read.query_name}, Position: {read_pos}, Base: {base}, Variants: {variants_at_position}")
                            print(forward_read.reference_name)
                            print(read_pos)
                            print(forward_read.reference_start)
                            
                            if(any(variant_at_pos.alts[0]=="C" for variant_at_pos in variants_at_position)):
                                print("false T>C +fw")
                            else:
                                print("true T>C +fw wrongvar")
                                conversions += 1                  
                        else:
                            print("true T>C +fw novar")
                            conversions += 1 

                for base in reverse_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if(None in base):
                        continue

                    read_seq = reverse_read.query_sequence[read_pos]

                    if read_seq == "G" and genome_seq == "a":
                        variants_at_position = list(vcffile.fetch(reverse_read.reference_name, genome_pos, genome_pos+1)) 
                        if variants_at_position:
                            print(reverse_read)
                            print(f"Read: {reverse_read.query_name}, Position: {read_pos}, Base: {base}, Variants: {variants_at_position}")
                            print(reverse_read.reference_name)
                            print(read_pos)
                            print(reverse_read.reference_start)
                            if(any(variant_at_pos.alts[0]=="G" for variant_at_pos in variants_at_position)):
                                print("false A>G +rv")
                            else:
                                print("true A>G +rv wrongvar")
                                conversions += 1                          

                        else:
                            print("true A>G +rv novar")
                            conversions += 1 

                conversion_counts.append(conversions)        
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
                            print(forward_read)
                            print(f"Read: {forward_read.query_name}, Position: {read_pos}, Base: {base}, Variants: {variants_at_position}")
                            print(forward_read.reference_name)
                            print(read_pos)
                            print(forward_read.reference_start)
                            if(any(variant_at_pos.alts[0]=="G" for variant_at_pos in variants_at_position)):
                                print("false A>G -fw")
                            else:
                                print("true A>G -fw wrongvar")
                                conversions += 1                           

                        else:
                            print("true A>G -fw novar")
                            conversions += 1 

                for base in reverse_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if(None in base):
                        continue

                    read_seq = reverse_read.query_sequence[read_pos]

                    if read_seq == "C" and genome_seq == "t":
                        variants_at_position = list(vcffile.fetch(reverse_read.reference_name, genome_pos, genome_pos+1)) 
                        if variants_at_position:
                            print(reverse_read)
                            print(f"Read: {reverse_read.query_name}, Position: {read_pos}, Base: {base}, Variants: {variants_at_position}")
                            print(reverse_read.reference_name)
                            print(read_pos)
                            print(reverse_read.reference_start)
                            if(any(variant_at_pos.alts[0]=="C" for variant_at_pos in variants_at_position)):
                                print("false T>C -rv")
                            else:
                                print("true T>C -rv wrongvar")
                                conversions += 1                           

                        else:
                            print("true T>C -rv novar")
                            conversions += 1 
                conversion_counts.append(conversions)   
            else:
                # should not be possible - but just in case
                pass 
            print("new pair")
                
    #print(conversion_counts)

    with open(args.outfile_txt, "w") as outfile:
        for count in conversion_counts:
            outfile.write(str(count)+"\n")
        print("written out to " + args.outfile_txt)

    # write footer and output benchmark information.
    E.stop()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
