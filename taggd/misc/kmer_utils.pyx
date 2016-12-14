# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
"""
Contains utilities for working with k-mer chunks of the barcodes or other sequences.
"""
cimport cython 
from cpython cimport bool
from collections import defaultdict

cdef object get_kmers_dicts(list seqs, int k, bool round_robin=False, int slider_increment=1):
    """
    Returns dictionaries for kmers of a list of sequences.
    The last kmer is always included, irrespective of slider increment.
    :param seqs: the sequences.
    :param k: the k-mer length.
    :param round_robin: if to treat the sequence as circular.
    :param slider_increment determines how the kmers are obtained
    :returns: a dictionary: kmer -> dictionary of seqs holding kmer -> 
    list of offsets of kmer in seq
    """
    cdef object kmer2seq = defaultdict(lambda : defaultdict(list))
    cdef str seq
    cdef str seqq
    cdef str kmer
    cdef int i
    for seq in seqs:
        # Adjust barcode if round robin
        seqq = seq + seq[0:(k-1)] if round_robin else seq
        # Create the Kmers of length k with slider increment
        for i in xrange(0, len(seqq)-k+1, slider_increment):
            kmer = seqq[i:i+k]
            kmer2seq[kmer][seq].append(i)
        # Special treatment of last in case we would skip it during incrementation
        if len(seqq) % slider_increment != 0:
            i = len(seqq)-k
            kmer = seqq[i:len(seqq)]
            kmer2seq[kmer][seq].add(i)       
    # Important to be able to generate KeyError
    kmer2seq.default_factory = None
    return kmer2seq

cdef list get_kmers(str seq, int k, bool round_robin=False, int slider_increment=0):
    """
    Returns the kmers of a sequence as a list of kmer-offset tuples. 
    The last kmer will always be included irrespective of the slider increment.
    :param seq: sequence.
    :param k: kmer length.
    :param round_robin: if to treat the sequence as circular.
    :param slider_increment determines how the kmers are obtained
    :return: the kmers list (kmer, offset).
    """
    cdef list kmer_list = list()
    cdef str seqq = seq + seq[0:(k-1)] if round_robin else seq
    cdef str kmer
    cdef int i
    # Simply compute kmers for the sequence
    #TODO this function could be used in get_kmers_dictst to avoid code duplication
    for i in xrange(0, len(seqq)-k+1, slider_increment):
        kmer = seqq[i:i+k]
        kmer_list.append((kmer, i))
    # Special treatment of last
    if len(seqq) % slider_increment != 0:
        i = len(seqq)-k
        kmer = seqq[i:len(seqq)]
        kmer_list.append((kmer,i))
    return kmer_list