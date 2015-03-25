import os
import argparse
import taggd.io.barcode_utils as bu
import taggd.core.demultiplex_core_functions as core
import taggd.core.demultiplex_search_functions as srch


def main(argv=None):
    """Main application."""

    # Create a parser
    parser = argparse.ArgumentParser(description=__doc__)

    # Needed parameters
    parser.add_argument('barcodes_infile', metavar='barcodes-infile', help="The file with true barcode IDs and other properties.")
    parser.add_argument('reads_infile', metavar='reads-infile', help="The FASTQ, FASTA, SAM or BAM file with reads.")
    parser.add_argument('outfile_prefix', metavar='outfile-prefix', help="The output files prefix.")

    # Optional arguments.
    parser.add_argument('--start-position', type=int, help='The start position for barcodes in reads (default: %(default)d)', default=0, metavar="[int]")
    parser.add_argument('--k', type=int, help='The kmer length (default: %(default)d)', default=7, metavar="[int]")
    parser.add_argument('--max-edit-distance', type=int, help='The max edit distance for allowing hits (default: %(default)d)', default=5, metavar="[int]")
    parser.add_argument('--metric', help= "Distance metric: Subglobal, Levenshtein or Hamming (default: %(default)s)", default="Subglobal", metavar="[string]")
    parser.add_argument('--slider-increment', type=int, help="Space between kmer searches, 0 yields kmer length (default: %(default)d)", default=0, metavar="[int]")
    parser.add_argument('--overhang', type=int, help="Additional flanking bases around read barcode to allow for insertions (default: %(default)d)", default=2, metavar="[int]")
    parser.add_argument('--only-output-matched', help="Suppresses writing of output to only file with matched reads.", default=False, action='store_true')
    parser.add_argument('--seed', help="Random number generator seed for shuffling ambiguous hits (default: %(default)s)", default=None, metavar="[string]")
    parser.add_argument('--no-multiprocessing', help="If set, turns off multiprocessing of reads", default=False, action='store_true')
    parser.add_argument('--estimate-min-edit-distance', type=int, help="If set, estimates the min edit distance among true barcodes by comparing the specified number of pairs, 0 means no estimation (default: %(default)d)", default=0, metavar="[int]")
    parser.add_argument('--max-chunk-size', type=int, help="Chunk of maximum number of simultaneously processed reads (default: %(default)d)", default=500000, metavar="[int]")
    parser.add_argument('--no-offset-speedup', help="Turns off an offset speedup routine, increasing time but may yield more hits.", default=False, action='store_true')

    # Parse
    if argv == None:
        options = parser.parse_args()
    else:
        options = parser.parse_args(argv)

    # Validate all options.
    if not os.path.isfile(options.barcodes_infile) :
        raise ValueError("Invalid true barcodes input file path.")
    if not os.path.isfile(options.reads_infile) :
        raise ValueError("Invalid reads input file path.")
    if not (options.reads_infile.upper().endswith(".FASTQ") or \
                    options.reads_infile.upper().endswith(".FQ") or \
                    options.reads_infile.upper().endswith(".SAM") or \
                    options.reads_infile.upper().endswith(".FASTA") or \
                    options.reads_infile.upper().endswith(".FA") or \
                    options.reads_infile.upper().endswith(".BAM")):
        raise ValueError("Invalid reads input file format: must be FASTQ, FASTA, SAM or BAM format and file end with .fq, fastq, .fa, .fasta, .sam or .bam")
    if options.outfile_prefix is None or options.outfile_prefix == "":
        raise ValueError("Invalid output file prefix.")
    if options.k <= 0:
        raise ValueError("Invalid kmer length. Must be > 0.")
    if options.max_edit_distance < 0:
        raise ValueError("Invalid max edit distance. Must be >= 0.")
    if options.metric not in ("Subglobal", "Levenshtein", "Hamming"):
        raise ValueError("Invalid metric. Must be Subglobal, Levenshtein or Hamming.")
    if options.slider_increment < 0:
        raise ValueError("Invalid slider increment. Must be >= 0.")
    if options.slider_increment == 0:
        options.slider_increment = int(options.k)
    if options.start_position < 0:
        raise ValueError("Invalid start position. Must be >= 0.")
    if options.overhang < 0:
        raise ValueError("Invalid overhang. Must be >= 0.")
    if options.metric == "Hamming" and options.overhang > 0:
        raise ValueError("Invalid overhang. Must be 0 for Hamming metric.")
    if options.max_chunk_size <= 0:
        raise ValueError("Invalid max chunk size. Must be > 0.")

    # Read barcodes file
    true_barcodes = bu.read_barcode_file(options.barcodes_infile)

    # Check barcodes file.
    if options.estimate_min_edit_distance > 0:
        min_dist = bu.estimate_min_edit_distance(true_barcodes, options.estimate_min_edit_distance)
        if min_dist <= options.max_edit_distance:
            raise ValueError("Invalid max edit distance: exceeds or equal to estimated minimum edit distance among true barcodes.")
        print "# Minimum edit distance between true barcodes was estimated to (may be less) " + str(min_dist)


    # Initialize.
    core.init(true_barcodes, options.reads_infile, os.path.abspath(options.outfile_prefix), options.max_edit_distance, \
              options.start_position, min(options.start_position, options.overhang), options.overhang, options.seed, \
              options.no_multiprocessing, options.only_output_matched, options.max_chunk_size)
    srch.init(true_barcodes, options.k, options.max_edit_distance, options.metric, options.slider_increment, \
              min(options.start_position, options.overhang), options.overhang, options.no_offset_speedup)

    # Demultiplex
    core.print_pre_stats()
    print "# Starting demultiplexing with the following options..." + str(options)
    core.demultiplex()
    print "# ...finished demultiplexing"
    core.print_post_stats()