"""
Main functions for carrying out the demultiplexing with
multithreading.
"""

import multiprocessing as mp
import random
import taggd.core.demultiplex_search_functions as srch
cimport taggd.core.demultiplex_search_functions as srch
from taggd.io.reads_reader_writer import *
from taggd.misc.counter import *
from cpython cimport bool


# Global variables for user options
cdef dict true_barcodes
cdef str reads_infile
cdef str outfile_prefix
cdef unsigned int max_edit_distance
cdef unsigned int start_position
cdef unsigned int barcode_length
cdef unsigned int pre_overhang
cdef unsigned int post_overhang
cdef str seed
cdef bool no_mp


# Reader-writer
cdef object reader_writer


# Global variables for stats
cdef object stats_total_reads
cdef object stats_unmatched
cdef object stats_perfect_matches
cdef object stats_imperfect_unambiguous_matches
cdef object stats_imperfect_ambiguous_matches
cdef list stats_edit_distance_counts


# Global variables for multithreading.
cdef object manager
cdef object q
cdef object pool
cdef object watcher


# Global variables for match types.
cdef int KILL, UNMATCHED, MATCHED_PERFECTLY, MATCHED_UNAMBIGUOUSLY, MATCHED_AMBIGUOUSLY


def init(\
        dict true_barcodes_, \
        str reads_infile_, \
        str outfile_prefix_, \
        unsigned int max_edit_distance_, \
        unsigned int start_position_, \
        unsigned int pre_overhang_, \
        unsigned int post_overhang_, \
        str seed_, \
        bool no_multiprocessing_):
    """
    Initializes settings (global variables).
    :param true_barcodes_: barcodes dict to attributes.
    :param reads_infile_: input file with read sequences (containing read barcodes).
    :param outfile_prefix_: output file prefix.
    :param max_edit_distance_: max allowed distance.
    :param start_position_: the start position.
    :param pre_overhang_: how many flanking bases pre for indels. Not applicable for Hamming distance.
    :param post_overhang_: how many flanking bases post for indels. Not applicable for Hamming distance.
    :param seed_: PRNG seed.
    :param no_multiprocessing: No multiprocessing.
    """

    global true_barcodes
    true_barcodes = true_barcodes_
    global reads_infile
    reads_infile = reads_infile_
    global outfile_prefix
    outfile_prefix = outfile_prefix_
    global max_edit_distance
    max_edit_distance = max_edit_distance_
    global start_position
    start_position = start_position_
    global pre_overhang
    pre_overhang = pre_overhang_
    global post_overhang
    post_overhang = post_overhang_
    global seed
    seed = seed_
    random.seed(seed)
    global no_multiprocessing
    no_multiprocessing = no_multiprocessing_

    # Read chip file
    global barcode_length
    barcode_length = len(true_barcodes.keys()[0])

    # Reader writer
    global reader_writer
    reader_writer = ReadsReaderWriter(reads_infile)

    # Stats
    global stats_total_reads
    stats_total_reads = Counter(0)
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
    cdef unsigned int i
    for i in xrange(max_edit_distance+1):
        stats_edit_distance_counts.append(Counter(0))

    # Match types
    global KILL
    global UNMATCHED
    global MATCHED_PERFECTLY
    global MATCHED_UNAMBIGUOUSLY
    global MATCHED_AMBIGUOUSLY
    KILL, UNMATCHED, MATCHED_PERFECTLY, MATCHED_UNAMBIGUOUSLY, MATCHED_AMBIGUOUSLY = range(-1,4)




def demultiplex():
    """Demultiplexes the contents of a reads file."""
    if no_multiprocessing:
        __demultiplex_linearly()
    else:
        __demultiplex_mp()



def __demultiplex_linearly():
    """Demultiplexes the contents of a reads file in a linear manner."""
    global manager
    manager = mp.Manager()
    global q
    q = manager.Queue()

    cdef object rec
    for rec in reader_writer.reader_open():
        __demultiplex_record_wrapper(rec)
    reader_writer.reader_close()
    q.put((None, KILL, None, None))
    __listener()



def __demultiplex_mp():
    """
    Demultiplexes the contents of a reads file in a multithreaded manner.
    """

    # Must use Manager queue here, or will not work
    global manager
    manager = mp.Manager()
    global q
    q = manager.Queue()
    global pool
    pool = mp.Pool(mp.cpu_count())

    # Put listener to work first
    global watcher
    watcher = pool.apply_async(__listener, ())

    # Fire off workers
    cdef list jobs = []
    cdef object job
    cdef object rec
    for rec in reader_writer.reader_open():
        job = pool.apply_async(__demultiplex_record_wrapper, (rec,))
        jobs.append(job)
    reader_writer.reader_close()

    # Collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    # Now we are done, kill the listener
    q.put((None, KILL, None, None))
    pool.close()



def __listener():
    '''Opens files, listens for messages on the q, writes output to files. '''

    # Open files.
    cdef object f_res, f_match, f_ambig, f_unmatch
    f_res = open(outfile_prefix + ".result", 'w')
    f_res.write("#Annotation\tMatch_result\tBarcode\tEdit_distance\tAmbiguous_top_hits\tQualified_candidates\tRaw_candidates\tLast_position\tApprox_insertions\tApprox_deletions\n")
    f_match = reader_writer.get_writer(outfile_prefix + ".matched." + reader_writer.get_format())
    f_ambig = reader_writer.get_writer(outfile_prefix + ".ambiguous." + reader_writer.get_format())
    f_unmatch = reader_writer.get_writer(outfile_prefix + ".unmatched." + reader_writer.get_format())

    cdef object rec
    cdef int result
    cdef str bcseq
    cdef object bc
    cdef list props
    cdef list tags
    cdef unsigned int i
    while True:
        # Extract record
        (rec, result, bcseq, props) = q.get()

        # End of line.
        if result == KILL:
            break

        # Write to info file.
        f_res.write("%s\t%s\t%s\t%s\n" % (rec.annotation, __match_type_to_str(result), bcseq, "\t".join(props)))

        # No match.
        if result == UNMATCHED:
            reader_writer.write_record(f_unmatch, rec)
            continue

        # Append record with properties. B0:Z:Barcode, B1:Z:Prop1, B2:Z:prop3 ...
        bc = true_barcodes[bcseq]
        tags = list()
        tags.append(("B0:Z", bcseq))
        for i in xrange(len(bc.attributes)):
            tags.append(("B" + str(i+1) + ":Z", bc.attributes[i]))
        rec.add_tags(tags)

        # Write to file.
        if result == MATCHED_PERFECTLY:
            reader_writer.write_record(f_match, rec)
        elif result == MATCHED_UNAMBIGUOUSLY:
            reader_writer.write_record(f_match, rec)
        elif result == MATCHED_AMBIGUOUSLY:
            reader_writer.write_record(f_ambig, rec)

    # Close all files.
    f_res.close()
    f_match.close()
    f_ambig.close()
    f_unmatch.close()


def __demultiplex_record_wrapper(object rec):
    """Wrapper for cdef:ed multithreading function."""
    return __demultiplex_record(rec)




cdef bool __demultiplex_record(object rec):
    """
    Demultiplexes a record.
    :return: true if mapped.
    """
    stats_total_reads.increment()

    cdef str read_barcode = rec.sequence[int(start_position):int(start_position+barcode_length)]
    cdef dict candidates = None
    cdef list qual_hits = None
    cdef list top_hits = None
    cdef str bcseq = None
    cdef unsigned int dist = 0
    cdef str last_pos = None
    cdef str a = None
    cdef str b = None
    if read_barcode in true_barcodes:
        # PERFECT MATCH
        # Record   Match type   Barcode   [Edit distance, Ambiguous top hits, Qualified candidates, Raw candidates, Last_pos, Insertions read, Insertions true]
        #print "PER_MATCH" + str(rec)
        q.put((rec, MATCHED_PERFECTLY, read_barcode, ["0", "1", "1", "-", str(barcode_length-1), "0", "0"]))
        stats_perfect_matches.increment()
        stats_edit_distance_counts[0].increment()
        return True
    else:
        # Include overhang.
        read_barcode = rec.sequence[int(start_position-pre_overhang):min(len(rec.sequence), int(start_position+barcode_length+post_overhang))]

        # Narrow down hits.
        candidates = srch.get_candidates(read_barcode)
        qual_hits = srch.get_distances(read_barcode, candidates)
        top_hits = srch.get_top_hits(qual_hits)

        if not top_hits:
            # NO MATCH.
            # Record   Match type   Barcode   [Edit distance, Ambiguous top hits, Qualified candidates, Raw candidates, Last_pos, Insertions read, Insertions true]
            q.put((rec, UNMATCHED, "-", ["-", "0", str(len(qual_hits)), str(len(candidates)), "-", "-", "-"]))
            stats_unmatched.increment()
            return False
        else:
            if len(top_hits) == 1:
                # SINGLE MATCH.
                bcseq, dist, last_pos, a, b = top_hits[0]
                # Record   Match type   Barcode   [Edit distance, Ambiguous top hits, Qualified candidates, Raw candidates, Last_pos, Insertions read, Insertions true]
                q.put((rec, MATCHED_UNAMBIGUOUSLY, bcseq, [str(dist), str(len(top_hits)), str(len(qual_hits)), str(len(candidates)), last_pos, a, b]))
                stats_imperfect_unambiguous_matches.increment()
                stats_edit_distance_counts[dist].increment()
                return True
            else:
                # AMBIGUOUS MATCHES.
                # Multiple best hits. Shuffle results to avoid biases.
                random.shuffle(top_hits)
                for (bcseq, dist, last_pos, a, b) in top_hits:
                    # Record   Match type   Barcode   [Edit distance, Ambiguous top hits, Qualified candidates, Raw candidates, Last_pos, Insertions read, Insertions true]
                    q.put((rec, MATCHED_AMBIGUOUSLY, bcseq, [str(dist), str(len(top_hits)), str(len(qual_hits)), str(len(candidates)), last_pos, a, b]))
                stats_imperfect_ambiguous_matches.increment()
                stats_edit_distance_counts[dist].increment()
                return True



def print_pre_stats():
    """Prints pre stats"""
    print "# Reads format: " + reader_writer.get_format()
    print "# True barcodes length: " + str(barcode_length)
    print "# Read barcodes length when overhang added: " + str(barcode_length + pre_overhang + post_overhang)



def print_post_stats():
    """Prints post stats"""
    print "# Total reads processed: " + str(stats_total_reads.value())
    cdef int matched_unam = stats_total_reads.value() - stats_unmatched.value() - stats_imperfect_ambiguous_matches.value()
    cdef float tot = float(stats_total_reads.value())
    print "# Reads matched unambiguously: " + str(matched_unam) + " [" + str(100 * matched_unam / tot) + "%]"
    print "#    - Reads matched perfectly: " + str(stats_perfect_matches.value()) + " [" + str(100* stats_perfect_matches.value() / tot) + "%]"
    print "#    - Reads matched imperfectly: " + str(stats_imperfect_unambiguous_matches.value()) + " [" + str(100* stats_imperfect_unambiguous_matches.value() / tot) + "%]"
    print "# Reads matched ambiguously: " + str(stats_imperfect_ambiguous_matches.value()) + " [" + str(100* stats_imperfect_ambiguous_matches.value() / tot) + "%]"
    print "# Reads unmatched: " + str(stats_unmatched.value()) + " [" + str(100* stats_unmatched.value() / tot) + "%]"
    cdef list distr = []
    cdef unsigned int i = 0
    for i in xrange(max_edit_distance+1):
        distr.append(str(stats_edit_distance_counts[i].value()))
    print("# Edit distance counts for [0,...," + str(max_edit_distance) + "]: [" + ", ".join(distr) + "]")



cdef str __match_type_to_str(int match_type):
    '''Simple converter.'''
    if match_type == UNMATCHED:
        return "UNMATCHED"
    elif match_type == MATCHED_PERFECTLY:
        return "MATCHED_PERFECTLY"
    elif match_type == MATCHED_UNAMBIGUOUSLY:
        return "MATCHED_UNAMBIGUOUSLY"
    elif match_type == MATCHED_AMBIGUOUSLY:
        return "MATCHED_AMBIGUOUSLY"
    return "KILL"