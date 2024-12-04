"""
all_introns_counts_and_info.py
====================
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

# mypy: ignore-errors

from argparse import ArgumentParser
import sys

import pandas as pd

from isoslam import isoslam, io


def main(argv=None):
    """Script main.
    parses command line options in sys.argv, unless *argv* is given.
    """
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "-b",
        "--bam",
        dest="infile_bam",
        type=str,
        help="Supply a path to the bam file that has undergone read assignment with featureCounts",
    )
    parser.add_argument(
        "-g", "--gtf", dest="gtf_path", type=str, help="Supply a path to the transcript assembly gtf file"
    )
    parser.add_argument(
        "-bed", dest="utron_bed", type=str, help="Supply a path to the utron bed file. Must be bed6 format"
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="outfile_tsv",
        type=str,
        help="""Supply a path to the output file. This file will contain
                        conversions per pair, accounting for stranding""",
    )
    parser.add_argument("-vcf", "--vcf", dest="vcf_path", type=str, help="""Supply a path to the VCF.gz file""")

    argv_as_dictionary = vars(argv)

    # Load files...
    bamfile = io.load_file(argv_as_dictionary["infile_bam"])
    vcffile = io.load_file(argv_as_dictionary["vcf_path"])
    # .bed file
    utron_coords = isoslam.extract_transcripts(argv_as_dictionary["utron_bed"])
    # .gtf file
    strand_dict, tx2gene = isoslam.extract_strand_transcript(argv_as_dictionary["gtf_path"])

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
    first_matched = 0
    i_output = 0

    with open(argv_as_dictionary["outfile_tsv"], "w") as outfile:
        # Add column headers
        outfile.write(
            "Read_UID\tTranscript_id\tStart\tEnd\tChr\tStrand\tAssignment\tConversions\tConvertible\tCoverage\n"
        )
        results = pd.DataFrame()

        for pair in fragment_iterator(bamfile.fetch(until_eof=True)):
            # if i_total_progress >= 2000000:
            ##    print("length break")
            #    break

            # if first_matched >1000:
            #    break

            if len(pair) != 2:
                continue

            read1, read2 = pair

            if read1.is_unmapped or read2.is_unmapped:
                continue

            i_total_progress += 1
            i_progress += 1

            if i_progress == 10000:
                i_progress = 0

            read1_start = read1.reference_start
            read1_end = read1.reference_end
            read2_start = read2.reference_start
            read2_end = read2.reference_end
            read1_length = read1.query_length
            read2_length = read2.query_length

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

            ## Get blocks to see how the reads are spliced
            ## Mainly use this to assign spliced reads
            ## But also use to check retained reads aren't spliced
            ## ALSO: There are some examples where a read would be called
            ## as retained because it overlaps an intron by 1 or 2 bases
            ## but these bases within the intron are the same as those after
            ## i.e. its ambiguous, but we can use evidence from the blocks of
            ## the partner read if it overlaps is spliced there.
            ## If not, we should cull assignment to those events
            ## We are comparing ret vs spliced, so instead of adding 1 to both
            ## we will just ignore those few
            block_starts1, block_ends1 = zip(*read1.get_blocks())
            block_starts2, block_ends2 = zip(*read2.get_blocks())
            # RETAINED
            read1_within_intron = {}

            for chr, start, end, transcript_id, strand in utrons1:
                # if starts/ends in intron/(s)
                if (
                    (start <= read1_start <= end or start <= read1_end <= end)
                    and start not in block_ends1
                    and end not in block_starts1
                    and start not in block_ends2
                    and end not in block_starts2
                ):
                    if transcript_id not in read1_within_intron:
                        read1_within_intron[transcript_id] = []
                    read1_within_intron[transcript_id].append((start, end, chr, strand))

                # if spans an entire intron
                intron_length = end - start
                if (
                    read1_start < start
                    and read1_end > end
                    and intron_length < read1_length
                    and start not in block_ends1
                    and end not in block_starts1
                    and start not in block_ends2
                    and end not in block_starts2
                ):
                    if transcript_id not in read1_within_intron:
                        read1_within_intron[transcript_id] = []
                    read1_within_intron[transcript_id].append((start, end, chr, strand))

            # SPLICED
            read1_spliced_3UI = {}

            for chr, start, end, transcript_id, strand in utrons1:
                if start in block_ends1 and end in block_starts1:
                    if transcript_id not in read1_spliced_3UI:
                        read1_spliced_3UI[transcript_id] = []
                    read1_spliced_3UI[transcript_id].append((start, end, chr, strand))

            ## READ 2
            # RETAINED
            read2_within_intron = {}

            for chr, start, end, transcript_id, strand in utrons2:
                # if starts/ends in intron/(s)
                if (
                    (start <= read2_start <= end or start <= read2_end <= end)
                    and start not in block_ends2
                    and end not in block_starts2
                    and start not in block_ends1
                    and end not in block_starts1
                ):
                    if transcript_id not in read2_within_intron:
                        read2_within_intron[transcript_id] = []
                    read2_within_intron[transcript_id].append((start, end, chr, strand))

                # if spans an entire intron
                intron_length = end - start
                if (
                    read2_start < start
                    and read2_end > end
                    and intron_length < read2_length
                    and start not in block_ends2
                    and end not in block_starts2
                    and start not in block_ends1
                    and end not in block_starts1
                ):
                    if transcript_id not in read2_within_intron:
                        read2_within_intron[transcript_id] = []
                    read2_within_intron[transcript_id].append((start, end, chr, strand))

            # SPLICED
            read2_spliced_3UI = {}

            for chr, start, end, transcript_id, strand in utrons2:
                if start in block_ends2 and end in block_starts2:
                    if transcript_id not in read2_spliced_3UI:
                        read2_spliced_3UI[transcript_id] = []
                    read2_spliced_3UI[transcript_id].append((start, end, chr, strand))

            all_dicts = [read1_within_intron, read2_within_intron, read1_spliced_3UI, read2_spliced_3UI]
            all_empty = all(not contents for contents in all_dicts)

            if all_empty:
                continue
            first_matched += 1

            # Create a set of tuples: (tx_id,(start,end))
            # Retained
            assign_conversions_to_retained = []

            for transcript_id, start_end_list in read1_within_intron.items():
                for start_end in start_end_list:
                    assign_conversions_to_retained.append((transcript_id, start_end))

            for transcript_id, start_end_list in read2_within_intron.items():
                for start_end in start_end_list:
                    assign_conversions_to_retained.append((transcript_id, start_end))

            # Only need each unique element once
            assign_conversions_to_retained = set(assign_conversions_to_retained)

            # Spliced
            assign_conversions_to_spliced = []

            for transcript_id, start_end_list in read1_spliced_3UI.items():
                for start_end in start_end_list:
                    assign_conversions_to_spliced.append((transcript_id, start_end))

            for transcript_id, start_end_list in read2_spliced_3UI.items():
                for start_end in start_end_list:
                    assign_conversions_to_spliced.append((transcript_id, start_end))

            # Only need each unique element once
            assign_conversions_to_spliced = set(assign_conversions_to_spliced)

            ## If there are any events in both we want to remove them - this should be rare
            in_both = assign_conversions_to_retained.intersection(assign_conversions_to_spliced)
            assign_conversions_to_spliced -= in_both
            assign_conversions_to_retained -= in_both

            i += 1

            read1_status = read1.get_tag("XS")
            read2_status = read2.get_tag("XS")

            status_list = list()
            status_list.append(read1_status)
            status_list.append(read2_status)

            if not (any(status in ["Assigned", "+", "-"] for status in status_list)):
                continue

            # pass if either is assigned
            if all(status in ["Assigned", "+", "-"] for status in status_list):
                # pass if both are assigned
                assignment1 = read1.get_tag("XT")
                assignment2 = read2.get_tag("XT")
                strand1 = strand_dict[assignment1]
                strand2 = strand_dict[assignment2]
                if strand1 == strand2:
                    # pass if both on same strand
                    strand = strand1
                else:
                    # if not on same strand bin the pair
                    continue
            else:
                # pass if only 1 is assigned
                if read1_status == "Assigned":
                    # if read 1 is the 1 assigned
                    assignment = read1.get_tag("XT")
                    strand = strand_dict[assignment]
                else:
                    # if read 2 is the 1 assigned]
                    assignment = read2.get_tag("XT")
                    strand = strand_dict[assignment]

            # assigned a "forward" and "reverse" read relative to the genome
            if read1.is_reverse and not read2.is_reverse:
                reverse_read = read1
                forward_read = read2
            elif read2.is_reverse and not read1.is_reverse:
                reverse_read = read2
                forward_read = read1
            else:
                # Not proper pair
                continue

            # if we are mapped to a +ve stranded transcript, then count T>C in
            # the forward read and A>G in the reverse read. if we are mapped to
            # a -ve stranded transcript, count T>C in the reverse read and A>G
            # in the forward read.
            if strand == "+":
                # pass if mapped to +ve transcript
                convertible = set()
                # create a set (list that only allows unique values to be added)
                # we will add the genome_pos at each point for both reads
                # len(coverage) will be the # of uniquely covered positions
                coverage = set()

                # instead of counting conversions as an int +=1, just add the
                # position to a set, and len(set) will be the number of
                # unique conversions. ACCOUNTS FOR OVERLAP BETWEEN READ1+2.
                converted_position = set()

                for base in forward_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if None in base:
                        continue

                    coverage.add(genome_pos)

                    read_seq = forward_read.query_sequence[read_pos]

                    if genome_seq.upper() == "T":
                        convertible.add(genome_pos)

                    if read_seq == "C" and genome_seq == "t":
                        variants_at_position = list(
                            vcffile.fetch(forward_read.reference_name, genome_pos, genome_pos + 1)
                        )
                        if variants_at_position:
                            if any(variant_at_pos.alts[0] == "C" for variant_at_pos in variants_at_position):
                                pass
                            else:
                                converted_position.add(genome_pos)
                        else:
                            converted_position.add(genome_pos)

                for base in reverse_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if None in base:
                        continue

                    coverage.add(genome_pos)

                    read_seq = reverse_read.query_sequence[read_pos]

                    if genome_seq.upper() == "A":
                        convertible.add(genome_pos)

                    if read_seq == "G" and genome_seq == "a":
                        variants_at_position = list(
                            vcffile.fetch(reverse_read.reference_name, genome_pos, genome_pos + 1)
                        )
                        if variants_at_position:
                            if any(variant_at_pos.alts[0] == "G" for variant_at_pos in variants_at_position):
                                pass
                            else:
                                converted_position.add(genome_pos)

                        else:
                            converted_position.add(genome_pos)

            elif strand == "-":
                # pass if mapped to -ve transcript
                convertible = set()
                coverage = set()
                converted_position = set()
                for base in forward_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if None in base:
                        continue

                    coverage.add(genome_pos)

                    read_seq = forward_read.query_sequence[read_pos]

                    if genome_seq.upper() == "A":
                        convertible.add(genome_pos)

                    if read_seq == "G" and genome_seq == "a":
                        variants_at_position = list(
                            vcffile.fetch(forward_read.reference_name, genome_pos, genome_pos + 1)
                        )
                        if variants_at_position:
                            if any(variant_at_pos.alts[0] == "G" for variant_at_pos in variants_at_position):
                                pass
                            else:
                                converted_position.add(genome_pos)

                        else:
                            converted_position.add(genome_pos)

                for base in reverse_read.get_aligned_pairs(with_seq=True):
                    read_pos, genome_pos, genome_seq = base
                    if None in base:
                        continue

                    coverage.add(genome_pos)

                    read_seq = reverse_read.query_sequence[read_pos]

                    if genome_seq.upper() == "T":
                        convertible.add(genome_pos)

                    if read_seq == "C" and genome_seq == "t":
                        variants_at_position = list(
                            vcffile.fetch(reverse_read.reference_name, genome_pos, genome_pos + 1)
                        )
                        if variants_at_position:
                            if any(variant_at_pos.alts[0] == "C" for variant_at_pos in variants_at_position):
                                pass
                            else:
                                converted_position.add(genome_pos)

                        else:
                            converted_position.add(genome_pos)
            else:
                # should not be possible - but just in case
                pass

            i_output += 1

            # Stream output as a tsv
            # Format: read_uid, transcript_id, start, end, ret/spl, conversions, convertible, coverage
            # A read pair will cover multiple lines if it matches multiple events (but metadata will be same)
            # ns-rse : Add in building Pandas dataframe so the function can return something that is testable
            for transcript_id, position in assign_conversions_to_retained:
                start, end, chr, strand = position
                outfile.write(
                    f"{i_output}\t{transcript_id}\t"
                    f"{start}\t{end}\t{chr}\t{strand}\tRet\t{len(converted_position)}\t"
                    f"{len(convertible)}\t{len(coverage)}\n"
                )
                row = pd.DataFrame(
                    [
                        {
                            "read_uid": i_output,
                            "transcript_id": transcript_id,
                            "start": start,
                            "end": end,
                            "chr": chr,
                            "strand": strand,
                            "assignment": "Ret",
                            "conversions": len(converted_position),
                            "convertible": len(convertible),
                            "coverage": len(coverage),
                        }
                    ]
                )
                results = pd.concat([results, row])

            for transcript_id, position in assign_conversions_to_spliced:
                start, end, chr, strand = position
                outfile.write(
                    f"{i_output}\t{transcript_id}\t"
                    f"{start}\t{end}\t{chr}\t{strand}\tSpl\t{len(converted_position)}\t"
                    f"{len(convertible)}\t{len(coverage)}\n"
                )
                row = pd.DataFrame(
                    [
                        {
                            "read_uid": i_output,
                            "transcript_id": transcript_id,
                            "start": start,
                            "end": end,
                            "chr": chr,
                            "strand": strand,
                            "assignment": "Spl",
                            "conversions": len(converted_position),
                            "convertible": len(convertible),
                            "coverage": len(coverage),
                        }
                    ]
                )
                results = pd.concat([results, row])

    return results.sort_values(by=["read_uid", "transcript_id", "chr", "start", "end"])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
