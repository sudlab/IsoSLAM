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
# pylint: disable=fixme

import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import polars as pl

from isoslam import io, isoslam


def main(argv=None):
    """Script main.
    parses command line options in sys.argv, unless *argv* is given.
    """
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", dest="config_file", type=Path, help="Path to configuration file to load.")
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
    # TODO: Add option for locating vcf file
    parser.add_argument(
        "-o",
        "--out",
        dest="outfile_tsv",
        type=str,
        help="""Supply a path to the output file. This file will contain
                        conversions per pair, accounting for stranding""",
    )
    parser.add_argument("-vcf", "--vcf", dest="vcf_path", type=str, help="Supply a path to the VCF.gz file")
    parser.add_argument("--delim", dest="delim", type=str, help="Delimiter to use in output.")

    argv_as_dictionary = vars(argv)
    print(f"{argv_as_dictionary=}")

    # Load files...
    # bamfile = io.load_file(argv_as_dictionary["infile_bam"])
    vcffile = io.load_file(argv_as_dictionary["vcf_path"])
    # .bed file
    utron_coords = isoslam.extract_transcripts(argv_as_dictionary["utron_bed"])
    # .gtf file
    strand_dict, tx2gene = isoslam.extract_strand_transcript(argv_as_dictionary["gtf_path"])

    # Load configuration
    config = io.load_and_update_config(argv)

    def fragment_iterator(read_iterator):
        read_list = []
        last_read = None

        for read in read_iterator:
            if last_read is not None and last_read != read.query_name:
                yield read_list
                read_list = []
                last_read = read.query_name
            last_read = read.query_name
            read_list.append(read)

        yield read_list

    i_progress = 0
    i_total_progress = 0
    first_matched = 0
    i_output = 0
    delim = argv_as_dictionary["delim"]
    # Create Polars dataframe to hold results
    schema = {
        "read_uid": int,
        "transcript_id": str,
        "start": int,
        "end": int,
        "chr": str,
        "strand": str,
        "assignment": str,
        "conversions": int,
        "convertible": int,
        "coverage": int,
    }
    results = pl.DataFrame(schema=schema)
    with open(argv_as_dictionary["outfile_tsv"], "w", encoding="utf-8") as outfile:
        # Add column headers
        outfile.write(
            # TODO: Lowercase and 'chromosome' over 'Chr'
            f"Read_UID{delim}Transcript_id{delim}Start{delim}End{delim}Chr{delim}Strand{delim}"
            f"Assignment{delim}Conversions{delim}Convertible{delim}Coverage\n"
        )
        pd_results = pd.DataFrame()

        for pair in isoslam.extract_segment_pairs(argv_as_dictionary["infile_bam"]):
            if i_total_progress >= 2000000:
                break
            if first_matched > 1000:
                break

            if len(pair) != 2:
                continue

            read1, read2 = pair

            if read1.is_unmapped or read2.is_unmapped:
                continue

            i_total_progress += 1
            i_progress += 1

            if i_progress == 10000:
                i_progress = 0

            # Extract features -
            # TODO - It should be possible to modify extract_features_from_pair() to determine which is forward and
            #        which reverse in this initial call so that subsequently we have pair_features["forward"] and
            #        pair_features["reverse"] that we can work with.
            pair_features = isoslam.extract_features_from_pair(pair)

            # Use the intersection of sets, this skips if either read1 or read2 aren't assigned, - or +
            if not {"Assigned", "+", "-"} & {pair_features["read1"]["status"], pair_features["read2"]["status"]}:
                continue

            # pass if either is assigned
            status_list = [pair_features["read1"]["status"], pair_features["read2"]["status"]]
            # If both reads are assigned/+/-...
            if all(status in ["Assigned", "+", "-"] for status in status_list):
                # ...check if strands are equal if so they are assign one to strand
                if (
                    strand_dict[pair_features["read1"]["transcript"]]
                    == strand_dict[pair_features["read2"]["transcript"]]
                ):
                    strand = strand_dict[pair_features["read1"]["transcript"]]
                # ...if not we pass this pair
                else:
                    continue
            else:
                # If read1 is Assigned we get the ??? for this transcript, otherwise it comes from read2
                strand = (
                    strand_dict[pair_features["read1"]["transcript"]]
                    if pair_features["read1"]["status"] == "Assigned"
                    else strand_dict[pair_features["read2"]["transcript"]]
                )
            # Set forward and reverse reads
            # TODO: Extract is_reverse to the pair_features["read1"]["reverse"] in refactored code
            if read1.is_reverse and not read2.is_reverse:
                reverse_read = read1
                forward_read = read2
            elif read2.is_reverse and not read1.is_reverse:
                reverse_read = read2
                forward_read = read1
            else:
                # Not proper pair
                continue

            # Lists (utrons1 and utrons2) are iterated over further down
            pair_features["read1"]["utron"] = isoslam.extract_utron(
                features=pair_features["read1"], gene_transcript=tx2gene, coordinates=utron_coords
            )
            pair_features["read2"]["utron"] = isoslam.extract_utron(
                features=pair_features["read2"], gene_transcript=tx2gene, coordinates=utron_coords
            )

            ## @ns-rse 2023-12-20
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
            block_starts1, block_ends1 = isoslam.zip_blocks(read1)
            block_starts2, block_ends2 = isoslam.zip_blocks(read2)

            # Build a dictionary as its cleaner to work with
            blocks = {
                "read1": {"starts": block_starts1, "ends": block_ends1},
                "read2": {"starts": block_starts2, "ends": block_ends2},
            }

            def print_info(chromosome, start, end, transcript_id, strand, blocks, text, count=i_total_progress) -> None:
                """Print information on progress to aid debugging."""
                print(f"\n{count} =========================== {text=}\n")
                print(f"{transcript_id=}")
                print(f"{chromosome=}")
                print(f"{start=}")
                print(f"{end=}")
                print(f"{strand=}")
                print(f"{blocks=}")

            # Retain within introns
            read1_within_intron = isoslam.filter_within_introns(pair_features, blocks, read="read1")
            read2_within_intron = isoslam.filter_within_introns(pair_features, blocks, read="read2")
            # Retain spliced
            read1_spliced_3UI = isoslam.filter_spliced_utrons(pair_features, blocks, read="read1")
            read2_spliced_3UI = isoslam.filter_spliced_utrons(pair_features, blocks, read="read2")

            all_dicts = [read1_within_intron, read2_within_intron, read1_spliced_3UI, read2_spliced_3UI]

            # Check that there are some retained regions (all dictionaries would be empty if there are none retained and
            # not {} evaluates to True, hence wrapping in all())
            if all(not contents for contents in all_dicts):
                continue
            first_matched += 1

            assign_conversions_to_retained = isoslam.unique_conversions(read1_within_intron, read2_within_intron)
            assign_conversions_to_spliced = isoslam.unique_conversions(read1_spliced_3UI, read2_spliced_3UI)

            ## If there are any events in both we want to remove them - this should be rare
            assign_conversions_to_retained, assign_conversions_to_spliced = isoslam.remove_common_reads(
                assign_conversions_to_retained, assign_conversions_to_spliced
            )

            # if we are mapped to a +ve stranded transcript, then count T>C in
            # the forward read and A>G in the reverse read.
            # TODO - Once we have assigned `read1/read2` to `forward/reverse` much earlier we can update this
            if strand == "+":
                coverage_counts = isoslam.count_conversions_across_pairs(
                    forward_read=forward_read,
                    reverse_read=reverse_read,
                    vcf_file=vcffile,
                    forward_conversion=config["forward_reads"],
                    reverse_conversion=config["reverse_reads"],
                )
            elif strand == "-":
                # If we are mapped to a -ve stranded transcript, count T>C in the reverse read and A>G
                # in the forward read.
                # can either flip the reads or the conversions
                coverage_counts = isoslam.count_conversions_across_pairs(
                    forward_read=reverse_read,
                    reverse_read=forward_read,
                    vcf_file=vcffile,
                    forward_conversion=config["forward_reads"],
                    reverse_conversion=config["reverse_reads"],
                )
            else:
                # should not be possible - but just in case
                pass

            i_output += 1

            # Stream output as a tsv
            # Format: read_uid, transcript_id, start, end, ret/spl, conversions, convertible, coverage
            io.write_assigned_conversions(
                assigned_conversions=assign_conversions_to_retained,
                coverage_counts=coverage_counts,
                read_uid=i_output,
                assignment="Ret",
                outfile=outfile,
                delim=argv_as_dictionary["delim"],
            )
            io.write_assigned_conversions(
                assigned_conversions=assign_conversions_to_spliced,
                coverage_counts=coverage_counts,
                read_uid=i_output,
                assignment="Spl",
                outfile=outfile,
                delim=argv_as_dictionary["delim"],
            )
            # A read pair will cover multiple lines if it matches multiple events (but metadata will be same)
            # ns-rse : Add in building Pandas dataframe so the function can return something that is testable
            results = isoslam.append_data(
                assigned_conversions=assign_conversions_to_retained,
                coverage_counts=coverage_counts,
                read_uid=i_output,
                assignment="Ret",
                results=results,
                schema=schema,
            )
            results = isoslam.append_data(
                assigned_conversions=assign_conversions_to_spliced,
                coverage_counts=coverage_counts,
                read_uid=i_output,
                assignment="Spl",
                results=results,
                schema=schema,
            )

    return results.sort(by=["read_uid", "transcript_id", "chr", "start", "end"])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
