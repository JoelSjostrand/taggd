"""
Contains functions for demultiplexing a record.
"""

import taggd.core.match as match
cimport taggd.core.match as match
import taggd.core.match_type as match_type
cimport taggd.core.match_type as match_type
import taggd.core.demultiplex_search_functions as srch
cimport taggd.core.demultiplex_search_functions as srch
import random
from cpython cimport bool

# User settings.
cdef dict true_barcodes
cdef int start_position
cdef int barcode_length
cdef int pre_overhang
cdef int post_overhang
cdef int max_edit_distance
cdef int homopolymer_filter
cdef list homopolymers
cdef str seed


def init(dict true_barcodes_,
         int start_position_,
         int pre_overhang_,
         int post_overhang_,
         int max_edit_distance_,
         int homopolymer_filter_,
         str seed_):
    """
    Initializes variables for matching.
    :param true_barcodes_: barcodes dict to attributes.
    :param start_position_: the start position.
    :param pre_overhang_: how many flanking bases pre for indels. Not applicable for Hamming distance.
    :param post_overhang_: how many flanking bases post for indels. Not applicable for Hamming distance.
    :param max_edit_distance_: max allowed distance.
    :param homopolymer_filter_: if containing homopolymer of this length, exclude. 0 means no filter.
    :param seed_: PRNG seed.
    """

    global true_barcodes
    true_barcodes = true_barcodes_

    # User settings.
    global start_position
    start_position = start_position_
    global pre_overhang
    pre_overhang = pre_overhang_
    global post_overhang
    post_overhang = post_overhang_
    global barcode_length
    barcode_length = len(true_barcodes.keys()[0])
    global max_edit_distance
    max_edit_distance = max_edit_distance_
    global homopolymer_filter
    homopolymer_filter = max(homopolymer_filter_, 0)
    global homopolymers
    homopolymers = list()
    for c in "ACGT":
        homopolymers.append(c * homopolymer_filter)
    global seed
    seed = seed_
    random.seed(seed)


def demultiplex_record_wrapper(object q, object rec):
    """
    Non cdef wrapper for cdef:ed multithreading function.
    """
    return demultiplex_record(q, rec)



cdef bool demultiplex_record(object q, object recs):
    """
    Demultiplexes a record. Init must have been called before.
    :param q: the queue where results are stored.
    :param recs: the records.
    :return: true if matched.
    """

    cdef str read_barcode = None
    cdef dict candidates = None
    cdef list qual_hits = None
    cdef list top_hits = None
    cdef str bcseq = None
    cdef int dist = 0
    cdef int last_pos = -1
    cdef int a = -1
    cdef int b = -1
    cdef object mtch = None
    cdef list mtchs = list()
    cdef bool homo = False

    # Iterate over records.
    for rec in recs:

        # Try perfect hit first.
        read_barcode = rec.sequence[start_position:(start_position+barcode_length)]
        if read_barcode in true_barcodes:
            mtch = match.Match(rec, match_type.MATCHED_PERFECTLY, read_barcode, 0, 1, 1, -1, barcode_length-1, 0, 0)
            mtchs.append(mtch)
            continue

        # Homopolymer filter.
        homo = False
        for filter in homopolymers:
            if filter in read_barcode:
                homo = True
                break
        if homo:
            mtch = match.Match(rec, match_type.UNMATCHED, "-", -1, 0, -1, -1, -1, -1, -1)
            mtchs.append(mtch)
            continue

        # Include overhang.
        read_barcode = rec.sequence[(start_position - pre_overhang):min(len(rec.sequence), \
                            (start_position + barcode_length + post_overhang))]

        # Narrow down hits.
        candidates = srch.get_candidates(read_barcode)
        qual_hits = srch.get_distances(read_barcode, candidates)
        top_hits = srch.get_top_hits(qual_hits)

        if not top_hits:
            # NO MATCH.
            # Record   Match type   Barcode
            # [Edit distance, Ambiguous top hits, Qualified candidates,
            # Raw candidates, Last_pos, Insertions read, Insertions true]
            mtch = match.Match(rec, match_type.UNMATCHED, "-", -1, 0, len(qual_hits), len(candidates), -1, -1, -1)
            mtchs.append(mtch)
        else:
            if len(top_hits) == 1:
                # SINGLE MATCH.
                bcseq, dist, last_pos, a, b = top_hits[0]
                # Record   Match type   Barcode
                # [Edit distance, Ambiguous top hits, Qualified candidates,
                # Raw candidates, Last_pos, Insertions read, Insertions true]
                mtch = match.Match(rec, match_type.MATCHED_UNAMBIGUOUSLY, bcseq, dist, len(top_hits), len(qual_hits),
                                   len(candidates), last_pos, a, b)
                mtchs.append(mtch)
            else:
                # AMBIGUOUS MATCHES.
                # Multiple best hits. Shuffle results to avoid biases.
                random.shuffle(top_hits)
                for (bcseq, dist, last_pos, a, b) in top_hits:
                    # Record   Match type   Barcode
                    # [Edit distance, Ambiguous top hits, Qualified candidates,
                    # Raw candidates, Last_pos, Insertions read, Insertions true]
                    mtch = match.Match(rec, match_type.MATCHED_AMBIGUOUSLY, bcseq, dist, len(top_hits),
                            len(qual_hits), len(candidates), last_pos, a, b)
                    mtchs.append(mtch)
    q.put(mtchs)
    return True



def print_pre_stats():
    """
    Prints pre stats
    """
    print "# True barcodes length: " + str(barcode_length)
    print "# Read barcodes length when overhang added: " + str(barcode_length + pre_overhang + post_overhang)