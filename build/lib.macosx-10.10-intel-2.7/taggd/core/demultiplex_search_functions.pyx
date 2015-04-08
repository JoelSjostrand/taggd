"""
Main functions for barcodes, like search and find matches
"""
cimport taggd.misc.kmer_utils as ku
cimport taggd.misc.distance_metrics as dm
from cpython cimport bool

# User options
cdef dict true_barcodes
cdef int k
cdef int max_edit_distance
cdef str metric
cdef int slider_increment
cdef int pre_overhang
cdef int post_overhang
cdef bool no_offset_speedup

# Dictionary
cdef dict kmer2seq = None

# Metrics
cdef int SUBGLOBAL
cdef int LEVENSHTEIN
cdef int HAMMING
cdef int metric_choice


def init(dict true_barcodes_,
        int k_,
        int max_edit_distance_,
        str metric_,
        int slider_increment_,
        int pre_overhang_,
        int post_overhang_,
        bool no_offset_speedup_):
    """
    Initializes settings (global variables).
    :param true_barcodes_: true barcodes dict to attributes.
    :param k_: the k-mer length used in the search algorithm.
    :param max_edit_distance_: max allowed distance.
    :param metric_: distance metric to use, Subglobal, Levenshtein or Hamming.
    :param slider_increment_: the number of increments for kmer search. The default, 0 sets it to k.
    :param pre_overhang: pre-read flanking bases.
    :param post_overhang: post-read flanking bases.
    """

    global true_barcodes
    true_barcodes = true_barcodes_
    global k
    k = k_
    global max_edit_distance
    max_edit_distance = max_edit_distance_
    global metric
    metric = metric_
    global slider_increment
    slider_increment = slider_increment_
    global pre_overhang
    pre_overhang = pre_overhang_
    global post_overhang
    post_overhang = post_overhang_
    global no_offset_speedup
    no_offset_speedup = no_offset_speedup_

    # Create k-mer mappings with ALL kmers
    global kmer2seq
    seq2kmer, kmer2seq = ku.get_kmers_dicts(true_barcodes.keys(), k, False, 1)

    # Metrics
    global SUBGLOBAL
    global LEVENSHTEIN
    global HAMMING
    SUBGLOBAL, LEVENSHTEIN, HAMMING = range(0, 3)
    global metric_choice
    if metric == "Subglobal":
        metric_choice = SUBGLOBAL
    elif metric == "Levenshtein":
        metric_choice = LEVENSHTEIN
    elif metric == "Hamming":
        metric_choice = HAMMING
    else:
        raise ValueError("Invalid distance metric specified")


cdef dict get_candidates(str read_barcode):
    """
    Returns candidate barcodes for a read barcode.
    :param read_barcode: read barcode.
    :return: the candidates as a dict with <barcode, read kmer hits>.
    """
    cdef dict candidates = dict()
    cdef list kmers_offsets = ku.get_kmers(read_barcode, k, False, slider_increment)
    cdef str kmer
    cdef int offset
    cdef dict hits
    cdef str hit
    cdef list hit_offsets
    cdef int hit_offset
    cdef int penalty = 0
    cdef int min_penalty = 0

    # NON-OPTIMIZED CASE
    # For each kmer in read (typically incremented by k positions at a time).
    if no_offset_speedup:
        for kmer, offset in kmers_offsets:
            try:
                hits = kmer2seq[kmer]
            except KeyError:
                continue
            # For each true barcode containing read's kmer.
            for hit, hit_offsets in hits.iteritems():
                candidates[hit] = 0
        return candidates

    # OPTIMIZED CASE
    # For each kmer in read (typically incremented by k positions at a time).
    for kmer, offset in kmers_offsets:
        # Obtain all the barcodes for a given kmer
        try:
            hits = kmer2seq[kmer]
        except KeyError:
            continue

        # For each true barcode containing read's kmer.
        for hit, hit_offsets in hits.iteritems():
            min_penalty = 100000000
            # For each position that kmer occurred in the true barcode.
            for hit_offset in hit_offsets:
                # Kmer may be shifted overhang positions without penalty, due to subglobal alignment.
                penalty = max(0, abs(offset - hit_offset) - pre_overhang - post_overhang)
                if penalty < min_penalty:
                    min_penalty = penalty
            #catching KeyError is faster than makin look-ups
            try:
                candidates[hit] = max(min_penalty, candidates[hit])
            except KeyError:
                candidates[hit] = min_penalty

    # Assure we remove the kmers
    del kmers_offsets
    
    # Clear out all candidates with a forced offset penalty greater than the max edit distance:
    for hit in candidates.keys():
        if candidates[hit] > max_edit_distance:
            del candidates[hit]
    return candidates


cdef list get_distances(str read_barcode, dict candidates):
    """
    Returns all qualified hits ordered by distance.
    :param read_barcode: the read barcode.
    :param candidates: the candidate dictionary.
    :return: a list of tuples, (barcode,distance,lastpos,ins,del).
    """
    cdef list qual_hits = []
    cdef str candidate = None
    cdef int penalty = 0
    cdef int dist = 0
    cdef int max_lim = 0
    cdef int read_last_pos = -1
    cdef int a = -1
    cdef int b = -1

    append = qual_hits.append
    for candidate, penalty in candidates.iteritems():
        if metric_choice == SUBGLOBAL:
            dist, read_last_pos, a, b = dm.subglobal_distance(read_barcode, candidate)
        elif metric_choice == LEVENSHTEIN:
            # Account for the added overhang!
            # Note: This does NOT equate subglobal and may miss cases subglobal would catch!
            max_lim = max_edit_distance + max(0, len(read_barcode) - len(candidate))
            dist = dm.levenshtein_distance(read_barcode, candidate, max_lim)
            read_last_pos = -1
            a = -1
            b = -1
        elif metric_choice == HAMMING:
            dist = dm.hamming_distance(read_barcode, candidate, max_edit_distance)
            read_last_pos = -1
            a = -1
            b = -1
        else:
            raise ValueError("Invalid distance metric specified")
        
        if dist <= max_edit_distance:
            append((candidate, dist, read_last_pos, a, b))
                
    return qual_hits


cdef list get_top_hits(list qual_hits):
    """
    Returns the top hits.
    :param qual_hits: the list of qualified hits.
    :return: the top hits.
    """
    if len(qual_hits) == 0:
        return None
    
    #TODO this can be optimized
    
    # Find smallest.
    cdef int mini = 10000000
    cdef str cand = None
    cdef int dist = 0
    cdef int last_pos = -1
    cdef int a = -1
    cdef int b = -1
    for cand, dist, last_pos, a, b in qual_hits:
        if dist < mini:
            mini = dist
    
    # Filter out the smallest
    cdef list top_hits = []
    append = top_hits.append
    for cand, dist, last_pos, a, b in qual_hits:
        if dist == mini:
            append((cand, dist, last_pos, a, b))
    return top_hits
