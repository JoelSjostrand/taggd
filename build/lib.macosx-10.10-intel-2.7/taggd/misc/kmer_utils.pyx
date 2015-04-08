"""
Contains utilities for working with k-mer chunks of the barcodes or other sequences.
"""
from cpython cimport bool

cdef list get_kmers_dicts(list seqs, int k, bool round_robin=False, int slider_increment=1):
    """
    Returns dictionaries for kmers of a list of sequences.
    The last kmer is always included, irrespective of slider increment.
    :param seqs: the sequences.
    :param k: the k-mer length.
    :param round_robin: if to treat the sequence as circular.
    :return: 1) a dictionary: seq -> set of kmer's of seq, and
             2) a dictionary: kmer -> dictionary of seqs holding kmer -> list of offsets of kmer in seq
    """
    cdef dict seq2kmer = {}
    cdef dict kmer2seq = {}
    cdef str seq
    cdef str seqq
    cdef set l
    cdef str kmer
    cdef int i
    for seq in seqs:
        l = set()
        seq2kmer[seq] = l
        seqq = seq
        if round_robin:
            seqq = seq + seq[0:(k-1)]
        for i in xrange(0, len(seqq)-k+1, slider_increment):
            kmer = seqq[i:i+k]
            l.add(kmer)
            if not kmer in kmer2seq:
                kmer2seq[kmer] = dict()
            if not seq in kmer2seq[kmer]:
                kmer2seq[kmer][seq] = list()
            kmer2seq[kmer][seq].append(i)
        # Special treatment of last in case we would skip it during incrementation
        if len(seqq) % slider_increment != 0:
            i = len(seqq)-k
            kmer = seqq[i:len(seqq)]
            l.add(kmer)
            if not kmer in kmer2seq:
                kmer2seq[kmer] = dict()
            if not seq in kmer2seq[kmer]:
                kmer2seq[kmer][seq] = list()
            kmer2seq[kmer][seq].add(i)
    return [seq2kmer, kmer2seq]


cdef list get_kmers(str seq, int k, bool round_robin=False, int slider_increment=1):
    """
    Returns the kmers of a sequence as a list of kmer-offset tuples. The last kmer will always be included
    irrespective of the slider increment.
    :param seq: sequence.
    :param k: kmer length.
    :param round_robin: if to treat the sequence as circular.
    :return: the kmers.
    """
    cdef list l = list()
    cdef str seqq = seq
    if round_robin:
        seqq = seq + seq[0:(k-1)]
    cdef str kmer
    cdef int i
    for i in xrange(0, len(seqq)-k+1, slider_increment):
        kmer = seqq[i:i+k]
        l.append((kmer, i))
    # Special treatment of last
    if len(seqq) % slider_increment != 0:
        i = len(seqq)-k
        kmer = seqq[i:len(seqq)]
        l.append((kmer,i))
    return l