""" 
Runs taggd, a tool to demultiplex (link molecular barcodes back to) a file of genetic reads,
typically obtained by sequencing. For matched reads, the barcode and its related properties
are added to the read. Unmatched reads, ambiguously matched reads, and stats are by default
produced as output files as well.

The input ID file should be tab-delimited with the following format:
<barcode>   <prop1>   <prop2>   ...
<barcode>   <prop1>   <prop2>   ...

The input files can be in FASTA, FASTQ, SAM or BAM format. Matched files will be appended
with the barcode and properties like this:

B0:Z:<barcode> B1:Z:<prop1> B2:Z:<prop3> ...

Source:          https://github.com/SpatialTranscriptomicsResearch/taggd
Python package:  https://pypi.python.org/pypi/taggd
Contact:         joel.sjostrand@gmail.com;jose.fernandez.navarro@scilifelab.se
"""
import os
import time
import multiprocessing as mp
import argparse
import taggd.io.barcode_utils as bu
import taggd.core.demultiplex_core_functions as core
import taggd.core.demultiplex_sub_functions as sub
import taggd.core.demultiplex_search_functions as srch

def main(argv=None):
    """
    Main application.
    Starts a timer, create parameter parsers, parsers parameters
    and run all the steps for the demultiplexing.
    """

    start_time = time.time()

    # Create a parser
    parser = argparse.ArgumentParser(description=__doc__, 
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Needed parameters
    parser.add_argument('barcodes_infile', 
                        metavar='barcodes-infile', 
                        help="The file with true barcode IDs and other properties.")
    parser.add_argument('reads_infile', 
                        metavar='reads-infile', 
                        help="The FASTQ, FASTA, SAM or BAM file with reads.")
    parser.add_argument('outfile_prefix', 
                        metavar='outfile-prefix', help="The output files prefix.")

    # Optional arguments.
    parser.add_argument('--no-matched-output',
                         help='Do not output matched reads',
                         default=False, action='store_true')
    parser.add_argument('--no-ambiguous-output',
                         help='Do not output ambiguous reads',
                         default=False, action='store_true')
    parser.add_argument('--no-unmatched-output',
                         help='Do not output unmatched reads',
                         default=False, action='store_true')
    parser.add_argument('--no-results-output',
                         help='Do not output a tab-separated results file with stats on the reads',
                         default=False, action='store_true')
    parser.add_argument('--start-position', 
                        type=int,
                         help='The start position for barcodes in reads (default: %(default)d)', 
                         default=0, metavar="[int]")
    parser.add_argument('--k', 
                        type=int, 
                        help='The kmer length (default: %(default)d)', 
                        default=6, metavar="[int]")
    parser.add_argument('--max-edit-distance', 
                        type=int, 
                        help='The max edit distance for allowing hits (default: %(default)d)', 
                        default=2, metavar="[int]")
    parser.add_argument('--metric', 
                        help= "Distance metric: Subglobal, Levenshtein or Hamming (default: %(default)s)", 
                        default="Subglobal", metavar="[string]")
    parser.add_argument('--ambiguity-factor',
                        type=float,
                        help='Top matches within this factor from the best match are considered ambiguous,\n'
                             'for instance with factor=1.5, having one match with distance 2 and two matches\n'
                             'with distance 4 yields all three matches as ambiguous hits. Perfect hits are always\n'
                             'considered non-ambiguous, irrespective of factor. (default: %(default)d)',
                        default=1.0, metavar="[int]")
    parser.add_argument('--slider-increment',
                        type=int, help="Space between kmer searches, " \
                        "0 yields kmer length (default: %(default)d)", 
                        default=0, metavar="[int]")
    parser.add_argument('--overhang', 
                        type=int, 
                        help="Additional flanking bases around read barcode\n" \
                        "to allow for insertions when matching (default: %(default)d)",
                        default=2, metavar="[int]")
    parser.add_argument('--seed', 
                        help="Random number generator seed for shuffling ambiguous hits (default: %(default)s)", 
                        default=None, metavar="[string]")
    parser.add_argument('--homopolymer-filter',
                        type=int,
                        help="If set, excludes reads where the barcode part contains\n" \
                        "a homopolymer of the given length,\n" \
                        "0 means no filter (default: %(default)d)",
                        default=8, metavar="[int]")
    parser.add_argument('--subprocesses',
                        type=int,
                        help="Number of subprocesses started (default: 0, yielding number of machine cores - 1)",
                        default=0, metavar="[int]")
    parser.add_argument('--estimate-min-edit-distance',
                        type=int, 
                        help="If set, estimates the min edit distance among true\n" \
                        "barcodes by comparing the specified number of pairs,\n" \
                        "0 means no estimation (default: %(default)d)", 
                        default=0, metavar="[int]")
    parser.add_argument('--no-offset-speedup', 
                        help="Turns off an offset speedup routine.\n" \
                        "Increases running time but may yield more hits.",
                        default=False, action='store_true')
    parser.add_argument('--multiple-hits-keep-one', 
                        help="When multiple kmer hits are found for a record\n" \
                        "keep one as unambiguous and the rest as ambiguous",
                        default=False, action='store_true')
    parser.add_argument('--trim-sequences', nargs='+', type=int, default=None, 
                        help="Trims from the barcodes in the input file\n" \
                        "The bases given in the list of tuples as START END START END .. where\n" \
                        "START is the integer position of the first base (0 based) and END is the integer\n" \
                        "position of the last base.\nTrimmng sequences can be given several times.")
    parser.add_argument('--version', action='version', version='%(prog)s ' + "0.3.1")

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
        raise ValueError("Invalid reads input file format: must be FASTQ, " \
                         "FASTA, SAM or BAM format and file end with .fq, fastq, .fa, .fasta, .sam or .bam")
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
    if options.subprocesses < 0:
        raise ValueError("Invalid no. of subprocesses. Must be >= 0.")
    if options.ambiguity_factor < 1.0:
        raise ValueError("Invalid ambiguity factor. Must be >= 1.")
    # Check the the trimming sequences given are valid
    if options.trim_sequences is not None \
    and (len(options.trim_sequences) % 2 != 0 or min(options.trim_sequences)) < 0:
        raise ValueError("Invalid trimming sequences given " \
                         "The number of positions given must be even and they must fit into the barcode length.")
        
    # Read barcodes file
    true_barcodes = bu.read_barcode_file(options.barcodes_infile)

    # Paths
    frmt = options.reads_infile.split(".")[-1]
    fn_bc = os.path.abspath(options.barcodes_infile)
    fn_reads = os.path.abspath(options.reads_infile)
    fn_prefix = os.path.abspath(options.outfile_prefix)
    fn_matched = None if options.no_matched_output else fn_prefix + "_matched." + frmt
    fn_ambig = None if options.no_ambiguous_output else fn_prefix + "_ambiguous." + frmt
    fn_unmatched = None if options.no_unmatched_output else fn_prefix + "_unmatched." + frmt
    fn_results = None if options.no_results_output else fn_prefix + "_results.tsv"

    # Subprocesses
    if options.subprocesses == 0:
        options.subprocesses = mp.cpu_count() - 1

    print "# Options: " + str(options).split("Namespace")[-1]
    print "# Barcodes input file: " + str(fn_bc)
    print "# Reads input file: " + str(fn_reads)
    print "# Matched output file: " + str(fn_matched)
    print "# Ambiguous output file: " + str(fn_ambig)
    print "# Unmatched output file: " + str(fn_unmatched)
    print "# Results output file: " + str(fn_results)
    print "# Number of barcodes in input: " + str(len(true_barcodes))
    lngth = len(true_barcodes.keys()[0])
    print "# Barcode length: " + str(lngth)
    print "# Barcode length when overhang added: " + \
    str(lngth + min(options.start_position, options.overhang) + options.overhang)
        
    # Check barcodes file.
    if options.estimate_min_edit_distance > 0:
        min_dist = estimate_min_edit_distance(true_barcodes, options.estimate_min_edit_distance)
        if min_dist <= options.max_edit_distance:
            raise ValueError("Invalid max edit distance: exceeds or equal " \
                             "to estimated minimum edit distance among true barcodes.")
        print "# Estimate of minimum edit distance between true barcodes (may be less): " + str(min_dist)
    else:
        print "# Estimate of minimum edit distance between true barcodes (may be less): Not estimated"

    # Make the input trim coordinates a list of tuples
    trim_sequences = None
    if options.trim_sequences is not None:
        trim_sequences = list()
        for i in xrange(len(options.trim_sequences) - 1):
            if i % 2 == 0:
                trim_sequences.append((options.trim_sequences[i], 
                                       options.trim_sequences[i+1]))
    
    # Initialize main components
    sub.init(true_barcodes,
             options.start_position,
             min(options.start_position, options.overhang),
             options.overhang,
             options.max_edit_distance,
             options.homopolymer_filter,
             options.seed,
             options.multiple_hits_keep_one,
             trim_sequences)

    srch.init(true_barcodes,
              options.k,
              options.max_edit_distance,
              options.metric,
              options.slider_increment, 
              min(options.start_position, options.overhang), 
              options.overhang,
              options.ambiguity_factor,
              options.no_offset_speedup)

    # Demultiplex
    print "# Starting demultiplexing..."
    stats = core.demultiplex(fn_reads,
                             fn_matched,
                             fn_ambig,
                             fn_unmatched,
                             fn_results,
                             options.subprocesses)
    print "# ...finished demultiplexing"
    print "# Wall time in secs: " + str(time.time() - start_time)
    print str(stats)
