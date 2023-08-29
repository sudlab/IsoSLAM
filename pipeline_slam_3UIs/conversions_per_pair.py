import pysam as pysam
from statistics import mean
from statistics import median
from matplotlib import pyplot as plt
import sys
import cgatcore.experiment as E

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-b", "--bam", dest="infile_bam", type=str,
                        help="Supply a path to the bam file you want the conversions of")

    parser.add_argument("-o", "--out", dest="outfile_txt", type=str,
                        help="Supply a path to the output file")

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args.infile_bam)
    #bamfile = pysam.AlignmentFile("../rates_nm60/name_sorted_bams/D2_65uM_EKRN230032564-1A_HGK2CDSX7_L3_slamdunk_mapped_sorted.bam")
    
    i = 0
    i_progress = 0

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
    

    for pair in fragment_iterator(bamfile.fetch(until_eof=True)):
    
        #if i > 100: 
        #    break

        if len(pair)!=2:
            continue

        read1_tc = 0
        read1_ag = 0

        read2_tc = 0
        read2_ag = 0

        read1, read2 = pair

        if read1.is_unmapped or read2.is_unmapped:
            continue 
        
        for base in read1.get_aligned_pairs(with_seq=True):
            read_pos, genome_pos, genome_seq = base
            if(None in base):
                continue

            read_seq = read1.query_sequence[read_pos]

            if read_seq == "C" and genome_seq == "t":
                read1_tc += 1 
                
            
            if read_seq == "G" and genome_seq == "a":
                read1_ag += 1 

        for base in read2.get_aligned_pairs(with_seq=True):
            read_pos, genome_pos, genome_seq = base
            if(None in base):
                continue

            read_seq = read2.query_sequence[read_pos]

            if read_seq == "C" and genome_seq == "t":
                read2_tc += 1 
                
            
            if read_seq == "G" and genome_seq == "a":
                read2_ag += 1 

        fw_tc = read1_tc + read2_ag
        rv_tc = read1_ag + read2_tc


        conversion_counts.append(max(fw_tc, rv_tc))
        i+=1
        i_progress += 1

        if i_progress == 100000:
            E.debug(str(i) + " pairs processed")
            i_progress = 0

    with open(args.outfile_txt, "w") as outfile:
        for count in conversion_counts:
            outfile.write(str(count)+"\n")
        print("written out to " + args.outfile_txt)

    # write footer and output benchmark information.
    E.stop()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
