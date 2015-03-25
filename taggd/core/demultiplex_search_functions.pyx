cimport taggd.misc.kmer_utils as ku
cimport taggd.misc.distance_metrics as dm
from cpython cimport bool


# User options
cdef dict true_barcodes
cdef unsigned int k
cdef unsigned int max_edit_distance
cdef str metric
cdef unsigned int slider_increment
cdef unsigned int pre_overhang
cdef unsigned int post_overhang
cdef bool no_offset_speedup


# Dictionary
cdef dict kmer2seq = None


# Metrics
cdef int SUBGLOBAL
cdef int LEVENSHTEIN
cdef int HAMMING
cdef int metric_choice



def init(\
        dict true_barcodes_, \
        unsigned int k_,
        unsigned int max_edit_distance_, \
        str metric_, \
        unsigned int slider_increment_, \
        unsigned int pre_overhang_, \
        unsigned int post_overhang_,\
        bool no_offset_speedup_\
    ):
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
    cdef unsigned int offset
    cdef dict hits
    cdef str hit
    cdef list hit_offsets
    cdef unsigned int hit_offset
    cdef unsigned int penalty = 0
    cdef unsigned int min_penalty = 0

    # NON-OPTIMIZED CASE
    # For each kmer in read (typically incremented by k positions at a time).
    if no_offset_speedup:
        for kmer, offset in kmers_offsets:
            if kmer in kmer2seq:
                hits = kmer2seq[kmer]
            else:
                continue
            # For each true barcode containing read's kmer.
            for hit, hit_offsets in hits.iteritems():
                candidates[hit] = 0
        return candidates

    # OPTIMIZED CASE
    # For each kmer in read (typically incremented by k positions at a time).
    for kmer, offset in kmers_offsets:
        if kmer in kmer2seq:
            hits = kmer2seq[kmer]
        else:
            continue

        # For each true barcode containing read's kmer.
        for hit, hit_offsets in hits.iteritems():
            min_penalty = 100000000
            # For each position that kmer occurred in the true barcode.
            for hit_offset in hit_offsets:
                # Kmer may be shifted overhang positions without penalty, due to subglobal alignment.
                penalty = max(0, abs(int(offset) - int(hit_offset)) - int(pre_overhang) - int(post_overhang))
                if penalty < min_penalty:
                    min_penalty = penalty
            if hit in candidates:
                candidates[hit] = max(min_penalty, candidates[hit])
            else:
                candidates[hit] = min_penalty

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
    cdef unsigned int penalty = 0
    cdef unsigned int dist = 0
    cdef unsigned int max_lim = 0
    cdef unsigned int read_last_pos = 0
    cdef unsigned int a = 0
    cdef unsigned int b = 0

    for candidate, penalty in candidates.iteritems():
        if metric_choice == SUBGLOBAL:
            dist, read_last_pos, a, b = dm.subglobal_distance(read_barcode, candidate)
            if dist <= max_edit_distance:
                qual_hits.append((candidate, dist, str(read_last_pos),  str(a), str(b)))
        elif metric_choice == LEVENSHTEIN:
            # Account for the added overhang!
            # Note: This does NOT equate subglobal and may miss cases subglobal would catch!
            max_lim = max_edit_distance + max(0, len(read_barcode) - len(candidate))
            dist = dm.levenshtein_distance(read_barcode, candidate, max_lim)
            if dist <= max_lim:
                qual_hits.append((candidate, min(dist, max_edit_distance), "-", "-", "-"))
        elif metric_choice == HAMMING:
            dist = dm.hamming_distance(read_barcode, candidate, max_edit_distance)
            if dist <= max_edit_distance:
                qual_hits.append((candidate, dist, "-", "-", "-"))
        else:
            raise ValueError("Invalid distance metric specified")
    return qual_hits




cdef list get_top_hits(list qual_hits):
    """
    Returns the top hits.
    :param qual_hits: the list of qualified hits.
    :return: the top hits.
    """
    if len(qual_hits) == 0:
        return None
    # Find smallest.
    cdef unsigned int mini = 10000000
    cdef str cand = None
    cdef unsigned int dist = 0
    cdef str a = None
    cdef str b = None
    for cand, dist, last_pos, a, b in qual_hits:
        if dist < mini:
            mini = dist
    # Filter out the smallest
    cdef list top_hits = []
    for cand, dist, last_pos, a, b in qual_hits:
        if dist == mini:
            top_hits.append((cand, dist, last_pos, a, b))
    return top_hits
