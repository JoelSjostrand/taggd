""" 
This class contains the main interface of taggd
Parameters will be captured, parsed and validated. 
The main workflow will be executed and output files
will be generated
"""
import os
import sys
from subprocess import Popen
import argparse
import taggd.io.barcode_utils as bu
import taggd.core.demultiplex_core_functions as core
import taggd.core.demultiplex_record_functions as rec
import taggd.core.demultiplex_search_functions as srch
import tempfile
import shutil
import pysam

def main(argv=None):
    """Main application."""

    # Create a parser
    parser = argparse.ArgumentParser(description=__doc__)

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
    parser.add_argument('--start-position', 
                        type=int,
                         help='The start position for barcodes in reads (default: %(default)d)', 
                         default=0, metavar="[int]")
    parser.add_argument('--k', 
                        type=int, 
                        help='The kmer length (default: %(default)d)', 
                        default=7, metavar="[int]")
    parser.add_argument('--max-edit-distance', 
                        type=int, 
                        help='The max edit distance for allowing hits (default: %(default)d)', 
                        default=5, metavar="[int]")
    parser.add_argument('--metric', 
                        help= "Distance metric: Subglobal, Levenshtein or Hamming (default: %(default)s)", 
                        default="Subglobal", metavar="[string]")
    parser.add_argument('--slider-increment', 
                        type=int, help="Space between kmer searches, " \
                        "0 yields kmer length (default: %(default)d)", 
                        default=0, metavar="[int]")
    parser.add_argument('--overhang', 
                        type=int, 
                        help="Additional flanking bases around read barcode " \
                        "to allow for insertions (default: %(default)d)", 
                        default=2, metavar="[int]")
    parser.add_argument('--only-output-matched', 
                        help="Suppresses writing of output to only file with matched reads.", 
                        default=False, action='store_true')
    parser.add_argument('--seed', 
                        help="Random number generator seed for shuffling ambiguous hits (default: %(default)s)", 
                        default=None, metavar="[string]")
    parser.add_argument('--no-multiprocessing',
                         help="If set, turns off multiprocessing of reads", 
                         default=False, action='store_true')
    parser.add_argument('--number-of-child-processes',
                        type=int, help="Number of child processes to process the input. " \
                        "This option must be used together with the option ----no-multiprocessing (default: %(default)d)",
                        default=0, metavar="[int]")
    parser.add_argument('--read-every-nth-entry-modulo',
                        type=int, help="Parse every nth entry of the input, " \
                        "The value 1, means handle all entries (default: %(default)d)",
                        default=1, metavar="[int]")
    parser.add_argument('--read-every-nth-entry-index',
                        type=int, help="Skip the first n entries, " \
                        "n should be less than --read-every-nth-entry-modulo (default: %(default)d)",
                        default=0, metavar="[int]")
    parser.add_argument('--homopolymer-filter',
                        type=int,
                        help="If set, excludes reads where the barcode part contains " \
                        "a homopolymer of the given length, " \
                        "0 means no filter (default: %(default)d)",
                        default=6, metavar="[int]")
    parser.add_argument('--estimate-min-edit-distance', 
                        type=int, 
                        help="If set, estimates the min edit distance among true " \
                        "barcodes by comparing the specified number of pairs, " \
                        "0 means no estimation (default: %(default)d)", 
                        default=0, metavar="[int]")
    parser.add_argument('--chunk-size',
                        type=int, 
                        help="Chunk of maximum number of simultaneously " \
                        "processed reads (default: %(default)d)", 
                        default=50000, metavar="[int]")
    parser.add_argument('--mp-chunk-size',
                        type=int,
                        help="Chunk of maximum number of " \
                        "processed reads sent to each thread (default: %(default)d)",
                        default=500, metavar="[int]")
    parser.add_argument('--no-offset-speedup', 
                        help="Turns off an offset speedup routine, " \
                        "increasing time but may yield more hits.", 
                        default=False, action='store_true')

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
    if options.chunk_size <= 0:
        raise ValueError("Invalid max chunk size. Must be > 0.")
    if options.mp_chunk_size <= 0:
        raise ValueError("Invalid multiprocessing max chunk size. Must be > 0.")
    if options.read_every_nth_entry_index >= options.read_every_nth_entry_modulo:
        raise ValueError("The value from --read-every-nth-entry-modulo must be bigger than the value from --read-every-nth-entry-index")
    if options.read_every_nth_entry_index < 0:
        raise ValueError("--read-every-nth-entry-index must not be negative")
    if not options.read_every_nth_entry_modulo > 0:
        raise ValueError("--read-every-nth-entry-modulo must be bigger than 0")
    if options.number_of_child_processes != 0 and not options.no_multiprocessing:
        raise ValueError("If --number-of-child-processes is not zero, the --no-multiprocessing must also be given")
    if options.number_of_child_processes != 0 and (options.read_every_nth_entry_modulo != 1):
        raise ValueError("--number-of-child-processes can't be used together with --read-every-nth-entry-modulo")

    # Read barcodes file
    true_barcodes = bu.read_barcode_file(options.barcodes_infile)

    # Check barcodes file.
    if options.estimate_min_edit_distance > 0:
        min_dist = bu.estimate_min_edit_distance(true_barcodes, options.estimate_min_edit_distance)
        if min_dist <= options.max_edit_distance:
            raise ValueError("Invalid max edit distance: exceeds or equal " \
                             "to estimated minimum edit distance among true barcodes.")
        print "# Minimum edit distance between true barcodes was estimated to (may be less) " + str(min_dist)

    if options.number_of_child_processes != 0:
      dirpath = tempfile.mkdtemp()
      processes = []
      barcodes_infile, reads_infile, outfile_prefix = sys.argv[-3:]
      arg_vector = sys.argv[:-3]
      num_children = options.number_of_child_processes
      arg_vector.append("--read-every-nth-entry-modulo")
      arg_vector.append(str(num_children))

      if arg_vector.count("--number-of-child-processes") > 1:
        raise ValueError("More than one --number-of-child-processes were given")
      assert arg_vector.count("--number-of-child-processes") != 0
      option_index = arg_vector.index("--number-of-child-processes")
      assert option_index < (len(arg_vector) -1)
      assert option_index > 0
      arg_vector_without_number_of_child_processes = arg_vector[:option_index] + arg_vector[option_index+2:]
      for i in range(num_children):
        arg_vector2 = arg_vector_without_number_of_child_processes[:]
        arg_vector2.append("--read-every-nth-entry-index")
        arg_vector2.append(str(i))
        arg_vector2.append(barcodes_infile)
        arg_vector2.append(reads_infile)
        arg_vector2.append(dirpath + "/" + str(i))
        print "arg_vector2 = " + ', '.join(arg_vector2)
        processes.append(Popen(arg_vector2))
      exit_codes = [p.wait() for p in processes]
      exit_codes_set = set(exit_codes)
      if (len(exit_codes_set) != 1 or 0 not in exit_codes_set):
        raise ValueError("At least one of the child processes failed")
      suffix = options.reads_infile.lower().split(".")[-1]
      for filename_ending in ['_matched.' + suffix, '_unmatched.' + suffix, '_ambiguous.' + suffix]:
        infile_template = pysam.AlignmentFile(dirpath  + "/" + str(0) + filename_ending, "r")
        outfile = pysam.AlignmentFile(outfile_prefix + filename_ending, "w", template=infile_template)
        for i in range(num_children):
          infile = pysam.AlignmentFile(dirpath  + "/" + str(i) + filename_ending, "r")
          for s in infile:
            outfile.write(s)

      filename_ending = "_results.tsv"
      with open(outfile_prefix + filename_ending,'wb') as outfile:
        for i in range(num_children):
          with open(dirpath  + "/" + str(i) + filename_ending, 'rb') as infile:
            shutil.copyfileobj(infile, outfile, 1024*1024*10)
      shutil.rmtree(dirpath)
      sys.exit(0)

    # Initialize main components
    core.init(true_barcodes,
              options.reads_infile,
              os.path.abspath(options.outfile_prefix),
              options.start_position,
              options.max_edit_distance,
              options.no_multiprocessing,
              options.read_every_nth_entry_modulo,
              options.read_every_nth_entry_index,
              options.only_output_matched,
              options.chunk_size,
              options.mp_chunk_size)

    rec.init(true_barcodes,
             options.start_position,
             min(options.start_position, options.overhang),
             options.overhang,
             options.max_edit_distance,
             options.homopolymer_filter,
             options.seed)

    srch.init(true_barcodes,
              options.k,
              options.max_edit_distance,
              options.metric,
              options.slider_increment, 
              min(options.start_position, options.overhang), 
              options.overhang,
              options.no_offset_speedup)

    # Demultiplex
    core.print_pre_stats()
    rec.print_pre_stats()
    print "# Starting demultiplexing with the following options..." + str(options)
    core.demultiplex()
    print "# ...finished demultiplexing"
    core.print_post_stats()
