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

    i = 0
    i_progress = 0

    conversion_counts = list()

    for read in bamfile.fetch(until_eof=True):
        #if i > 100000: 
        #    break

        local_tc = 0
        local_ag = 0

        if read.is_unmapped:
            continue 
        
        for base in read.get_aligned_pairs(with_seq=True):
            read_pos, genome_pos, genome_seq = base
            if(None in base):
                continue

            read_seq = read.query_sequence[read_pos]

            if read_seq == "C" and genome_seq == "t":
                local_tc += 1 
                
            
            if read_seq == "G" and genome_seq == "a":
                local_ag += 1 


        conversion_counts.append(max(local_tc, local_ag))

        i+=1
        i_progress+=1

        if i_progress == 100000:
            E.debug(str(i) + " reads processed")
            i_progress = 0

    with open(args.outfile_txt, "w") as outfile:
        for count in conversion_counts:
            outfile.write(str(count)+"\n")
        print("written out to " + args.outfile_txt)

    # write footer and output benchmark information.
    E.stop()



if __name__ == "__main__":
    sys.exit(main(sys.argv))
