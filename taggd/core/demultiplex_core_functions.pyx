"""
Main functions for carrying out the demultiplexing with
multithreading.
"""
import os
import time
import multiprocessing as mp
import taggd.core.match as match
cimport taggd.core.match as match
import taggd.core.match_type as match_type
cimport taggd.core.match_type as match_type
import taggd.core.demultiplex_record_functions as demulti
import taggd.io.reads_reader_writer as rw
from taggd.misc.counter import *
from cpython cimport bool

# Global variables for user options
cdef dict true_barcodes
cdef str reads_infile
cdef str outfile_prefix
cdef int max_edit_distance
cdef bool no_multiprocessing
cdef bool only_output_matched
cdef int chunk_size
cdef int start_position
cdef int barcode_length

# Reader-writer
cdef object re_wr
cdef object f_match = None
cdef object f_res = None
cdef object f_ambig = None
cdef object f_unmatch = None


# Global variables for stats.
cdef object stats_total_reads
cdef object stats_total_reads_wr
cdef object stats_unmatched
cdef object stats_perfect_matches
cdef object stats_imperfect_unambiguous_matches
cdef object stats_imperfect_ambiguous_matches  # Non-unique
cdef list stats_edit_distance_counts
cdef int start_time

# Threading objects
cdef object manager


def init(dict true_barcodes_,
         str reads_infile_,
         str outfile_prefix_,
         int start_position_,
         int max_edit_distance_,
         bool no_multiprocessing_,
         bool only_output_matched_,
         int chunk_size_,
         int mp_chunk_size_):
    """
    Initializes settings (global variables).
    :param true_barcodes_: true barcodes.
    :param reads_infile_: input file with read sequences (containing read barcodes).
    :param outfile_prefix_: output file prefix.
    :param max_edit_distance_: Max edit distance.
    :param no_multiprocessing_: No multiprocessing.
    :param only_output_matched_: Only matched file out.
    """

    global true_barcodes
    true_barcodes = true_barcodes_
    global reads_infile
    reads_infile = reads_infile_
    global outfile_prefix
    outfile_prefix = outfile_prefix_
    global start_position
    start_position = start_position_
    global barcode_length
    barcode_length = len(true_barcodes.keys()[0])
    global max_edit_distance
    max_edit_distance = max_edit_distance_
    global no_multiprocessing
    no_multiprocessing = no_multiprocessing_
    global only_output_matched
    only_output_matched = only_output_matched_
    global chunk_size
    chunk_size = chunk_size_
    global mp_chunk_size
    mp_chunk_size = mp_chunk_size_

    # Reader writer
    global re_wr
    re_wr = rw.ReadsReaderWriter(reads_infile)

    # Stats
    global stats_total_reads
    stats_total_reads = Counter(0)
    global stats_total_reads_wr
    stats_total_reads_wr = Counter(0)
    global stats_unmatched
    stats_unmatched = Counter(0)
    global stats_perfect_matches
    stats_perfect_matches = Counter(0)
    global stats_imperfect_unambiguous_matches
    stats_imperfect_unambiguous_matches = Counter(0)
    global stats_imperfect_ambiguous_matches
    stats_imperfect_ambiguous_matches = Counter(0)
    global stats_edit_distance_counts
    stats_edit_distance_counts = []
    cdef int i
    for i in xrange(max_edit_distance+1):
        stats_edit_distance_counts.append(Counter(0))
    global start_time
    start_time = time.time()

    # Threading.
    global manager
    manager = mp.Manager()


def demultiplex():
    """
    Demultiplexes the contents of a reads file.
    It tries to find all the matches using the reads and barcodes
    present in global variables and it will write the results
    to different files.
    NOTE that this functions expects init() to have been called before.
    """

    # Open files.
    __open_files()

    # Demultiplex.
    cdef list chunk = list()
    cdef object rec

    # Process the reads in chunks either parallel or in single thread mode
    # Reference to function append to avoid overhead.
    append = chunk.append
    for rec in re_wr.reader_open():
        append(rec)
        stats_total_reads.increment()
        if stats_total_reads.value() % chunk_size == 0:
            if no_multiprocessing:
                __demultiplex_linearly_chunk(chunk)
            else:
                __demultiplex_mp_chunk(chunk)
            del chunk[:]
    # Process the rest if applies
    if no_multiprocessing:
        __demultiplex_linearly_chunk(chunk)
    else:
        __demultiplex_mp_chunk(chunk)
    del chunk[:]

    # Close all files.
    __close_files()


def __open_files():
    """
    Opens files for writing.
    """
    global f_match
    f_match = re_wr.get_writer(outfile_prefix + "_matched." + re_wr.get_format())
    if not only_output_matched:
        f_res_name = outfile_prefix + "_results.tsv"
        if os.path.exists(f_res_name):
            os.remove(f_res_name)
        global f_res
        f_res = open(f_res_name, 'w')
        f_res.write(match.get_match_header() + "\n")
        f_res.flush()
        global f_ambig
        f_ambig = re_wr.get_writer(outfile_prefix + "_ambiguous." + re_wr.get_format())
        global f_unmatch
        f_unmatch = re_wr.get_writer(outfile_prefix + "_unmatched." + re_wr.get_format())


def __close_files():
    """
    Closes output files.
    """
    re_wr.reader_close()
    f_match.close()
    if not only_output_matched:
        f_res.close()
        f_ambig.close()
        f_unmatch.close()


def __demultiplex_linearly_chunk(list chunk):
    """
    Demultiplexes a list of reads in a single-threaded approach
    """
    cdef object q = manager.Queue()

    cdef object rec = None
    # this loop will put messages in the queue containing
    # records to write
    demulti.demultiplex_record_wrapper(q, chunk)
    # write matches
    __write_matches(q)


def __demultiplex_mp_chunk(list chunk):
    """
    Demultiplexes a list of reads in a multi-threaded approach
    """

    # Must use Manager queue here, or will not work
    cdef object q = manager.Queue()
    cdef object pool = mp.Pool(mp.cpu_count() - 1)

    # Fire off workers
    cdef list jobs = []
    cdef object job = None
    cdef list recs = list()
    cdef object rec = None
    cdef str read_barcode = None
    cdef object mtch = None
    cdef int i = 0
    append = jobs.append
    for rec in chunk:
        recs.append(rec)
        i += 1
        if i % mp_chunk_size == 0:
            job = pool.apply_async(demulti.demultiplex_record_wrapper, (q, recs,))
            append(job)
            recs = list()
    job = pool.apply_async(demulti.demultiplex_record_wrapper, (q, recs,))
    append(job)

    # Collect results from the workers through the pool result queue.
    for job in jobs:
        job.get()

    pool.close()
    pool.join()

    # write matches
    __write_matches(q)


def __write_matches(object q):
    """
    Processes matches in queue to write results.
    """

    cdef list mtchs = None
    cdef object mtch = None
    cdef object rec = None
    cdef str bcseq = None
    cdef object bc = None
    cdef list props = None
    cdef list tags = None
    cdef int i = 0

    while not q.empty():

        # Extract records
        mtchs = q.get()

        for mtch in mtchs:

            # Write to info file.
            if not only_output_matched:
                f_res.write(str(mtch) + "\n")

            # No match.
            if mtch.match_type == match_type.UNMATCHED:
                if not only_output_matched:
                    re_wr.write_record(f_unmatch, mtch.record)
                    stats_total_reads_wr.increment()
                stats_unmatched.increment()
                continue

            # Append record with properties. B0:Z:Barcode, B1:Z:Prop1, B2:Z:prop3 ...
            bc = true_barcodes[mtch.barcode]
            tags = list()
            tags.append(("B0:Z", mtch.barcode))
            for i in xrange(len(bc.attributes)):
                tags.append(("B" + str(i+1) + ":Z", bc.attributes[i]))
            mtch.record.add_tags(tags)

            # Write to file.
            if mtch.match_type == match_type.MATCHED_PERFECTLY:
                re_wr.write_record(f_match, mtch.record)
                stats_perfect_matches.increment()
                stats_edit_distance_counts[0].increment()
                stats_total_reads_wr.increment()
            elif mtch.match_type == match_type.MATCHED_UNAMBIGUOUSLY:
                re_wr.write_record(f_match, mtch.record)
                stats_imperfect_unambiguous_matches.increment()
                stats_edit_distance_counts[mtch.edit_distance].increment()
                stats_total_reads_wr.increment()
            elif mtch.match_type == match_type.MATCHED_AMBIGUOUSLY:
                if not only_output_matched:
                    re_wr.write_record(f_ambig, mtch.record)
                    stats_total_reads_wr.increment()
                stats_imperfect_ambiguous_matches.increment()


def print_pre_stats():
    """
    Prints pre stats
    """
    print "# Absolute output prefix path: " + outfile_prefix
    print "# Only writing matched reads: " + str(only_output_matched)
    print "# Reads format: " + re_wr.get_format()


def print_post_stats():
    """
    Prints post stats
    """
    print "# Total reads processed from infile: " \
        + str(stats_total_reads.value())
    print "# Total reads written (including multiple ambiguities): " \
        + str(stats_total_reads_wr.value())
    cdef int matched_unam = stats_perfect_matches.value() + stats_imperfect_unambiguous_matches.value()
    cdef float tot = float(stats_total_reads.value())
    print "# Reads matched unambiguously: " \
        + str(matched_unam) + " [" + str(matched_unam / tot * 100) + "%]"
    print "#    - Reads matched perfectly: " \
        + str(stats_perfect_matches.value()) + " [" + str(stats_perfect_matches.value() / tot * 100) + "%]"
    print "#    - Reads matched imperfectly: " \
        + str(stats_imperfect_unambiguous_matches.value()) \
        + " [" + str(stats_imperfect_unambiguous_matches.value() / tot * 100) + "%]"
    cdef int uniq_amb = stats_total_reads.value() - stats_perfect_matches.value() - \
                        stats_imperfect_unambiguous_matches.value() - stats_unmatched.value()
    print "# Reads matched ambiguously: " + str(uniq_amb) + " unique (making " \
        + str(stats_imperfect_ambiguous_matches.value()) + " overall) " \
        + " [" + str(uniq_amb / tot * 100) + "%]"
    print "# Reads unmatched: " \
        + str(stats_unmatched.value()) + " [" + str(stats_unmatched.value() / tot * 100) + "%]"
    cdef list distr = []
    cdef int i = 0
    for i in xrange(max_edit_distance + 1):
        distr.append(str(stats_edit_distance_counts[i].value()))
    print("# Edit distance counts for [0,...," + str(max_edit_distance) + "]: [" + ", ".join(distr) + "]")
    cdef int wall_time = time.time() - start_time
    cdef int days = wall_time / (24*60*60)
    cdef int hours = (wall_time - (24*60*60*days)) / (60*60)
    cdef int mins = (wall_time - (24*60*60*days) - (60*60*hours)) / 60
    cdef int secs = (wall_time - (24*60*60*days) - (60*60*hours) - (60*mins))
    print("# Wall time: %d:%d:%d:%d" % (days, hours, mins, secs))

