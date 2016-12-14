"""
Main functions for barcodes, like search and find matches
"""
cimport taggd.misc.kmer_utils as ku
cimport taggd.misc.distance_metrics as dm
from cpython cimport bool
from collections import defaultdict

# User options
cdef int k
cdef int max_edit_distance
cdef str metric
cdef int slider_increment
cdef int pre_overhang
cdef int post_overhang
cdef float ambiguity_factor
cdef bool no_offset_speedup

# Dictionary of kmer to sequences
cdef object kmer2seq = None

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
        float ambiguity_factor_,
        bool no_offset_speedup_):
    """
    Initializes settings (global variables).
    """

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
    global ambiguity_factor
    ambiguity_factor = ambiguity_factor_
    global no_offset_speedup
    no_offset_speedup = no_offset_speedup_

    # Create k-mer mappings with ALL kmers
    global kmer2seq
    kmer2seq = ku.get_kmers_dicts(true_barcodes_.keys(), k, False)

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

cdef list get_candidates(str read_barcode):
    """
    Returns candidate barcodes for a read barcode
    as a list of barcodes
    First, split the barcode in kmers and then
    finds all the candidate barcodes of the kmers
    :param read_barcode the barcode from which to get candidates
    :return a list of candidates barcodes
    """
    # NOTE probably faster to keep kmer_offsets in memory as we will call
    #      this function several times with the same barcode but we get a penalty in memory use
    cdef object candidates = defaultdict(int)
    cdef list kmers_offsets = ku.get_kmers(read_barcode, k, False, slider_increment)
    cdef int penalty = 0
    cdef int min_penalty = 0
    cdef str kmer
    cdef int offset
    cdef object hits
    cdef int hit_offset
    cdef str hit
    cdef list hit_offsets
    
    # Iterate all the kmer-offset combinations found in the input barcode
    for kmer, offset in kmers_offsets:
        # Obtain all the barcodes that matched for the current kmer
        try:
            hits = kmer2seq[kmer]
        except KeyError:
            continue
        # For each true barcode containing read's kmer.
        # Hit refers to barcode and hit_offsets to where the kmer was in the barcode
        for hit, hit_offsets in hits.iteritems():
            if no_offset_speedup:
                # NON-OPTIMIZED CASE
                # For each kmer in read (typically incremented by k positions at a time).
                candidates[hit] = 0
                continue
            # OPTIMIZED CASE
            # For each kmer in read (typically incremented by k positions at a time).
            min_penalty = 100000000
            # For each position that kmer occurred in the true barcode.
            # Get the minimum penalty
            for hit_offset in hit_offsets:
                # Kmer may be shifted overhang positions without penalty, due to subglobal alignment.
                penalty = max(0, abs(offset - hit_offset) - pre_overhang - post_overhang)
                if penalty < min_penalty: min_penalty = penalty
            # Assign the min penalty to the candidate (if exists already take max)
            # TODO if there are several equal barcode candidates for different kmers,
            #      why keep the max penalty and not an average?
            candidates[hit] = max(min_penalty, candidates[hit])
            
    # Clear out all candidates with a forced offset penalty greater than the max edit distance and return
    return [hit for hit,penal in candidates.iteritems() if penal <= max_edit_distance]
    
cdef list get_distances(str read_barcode, list candidates):
    """
    Returns all qualified hits ordered by distance as
    a list of tuples, (barcode,distance).
    :param read_barcode the original barcode
    :param candidates a list of possible candidates
    :return a list of (candidate, distance score)
    """
    cdef list qual_hits = []
    cdef int dist = 0
    cdef int max_lim = 0
    cdef int read_last_pos = -1
    cdef int a = -1
    cdef int b = -1
    cdef str candidate
    cdef int length_barcode = len(read_barcode)
    # Iterate candidates, compute distance with the original barcode
    # and create a list of candidate hits
    for candidate in candidates:
        if metric_choice == SUBGLOBAL:
            dist = dm.subglobal_distance(read_barcode, candidate)
        elif metric_choice == LEVENSHTEIN:
            # Account for the added overhang!
            # Note: This does NOT equate subglobal and may miss cases subglobal would catch!
            max_lim = max_edit_distance + max(0, length_barcode - len(candidate))
            dist = dm.levenshtein_distance(read_barcode, candidate, max_lim)
        else:
            dist = dm.hamming_distance(read_barcode, candidate, max_edit_distance)
        # Only add if distance is good
        if dist <= max_edit_distance: qual_hits.append((candidate, dist))
    return qual_hits

cdef list get_top_hits(list qual_hits):
    """
    Returns the top hits candidates filtering by minimum distance
    :param qual_hits the list of possible candidate tuples (barcode,distance)
    :return the filtered list
    """
    if len(qual_hits) == 0: return None
 
    # Find smallest distance
    from operator import itemgetter
    sorted_qual_hits = sorted(qual_hits, key=itemgetter(1), reverse=False)
    cdef int mini = round(sorted_qual_hits[0][1] * ambiguity_factor)
    # Filter out elements with similar or lower distance than the min*ambiguity_factor
    return [candidate for candidate in qual_hits if candidate[1] <= mini]
