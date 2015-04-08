"""
This package contains utilities for the Phred quality measurement.
By definition, Phred quality score Q is related to base-calling error probability P
as follows:  Q = -10 log10 P, or alternatively P = 10^{-Q/10}.
"""

import math

cpdef double probability_to_phred(double p):
    """
    Returns the Phred score of a probability.
    """
    if p < 0 or p > 1:
        raise ValueError("Probability must be in range [0,1]")
    return (-10 * math.log10(p))

cpdef double phred_to_probability(double q):
    """
    Returns the probability of a Phred score.
    """
    if q < 0:
        raise ValueError("Phred score must be in [0,inf).")
    return (math.pow(10, (-q/10)))


cpdef str sanger_phred_to_ascii(list phred):
    """
    Converts a list of sanger phred sequences to an ascii string
    :param phred: the list.
    :return: the string.
    """
    cdef list s = []
    cdef int p
    cdef int q
    for p in phred:
        q = p + 33
        if q < 33 or q > 126:
            raise ValueError("Invalid Sanger ASCII value: " + q + ". Should be in range [33,126].")
        s.append(chr(q))
    return ''.join(s)


cpdef list sanger_ascii_to_phred(str ascii_seq):
    """
    Converts an ASCII string from a SAM file or FASTQ file with Sanger Phred format
    to a list of phred scores
    """
    # Ascii range [33,126], phred range [0,93].
    cdef list phred = []
    cdef char c
    cdef int v
    for c in ascii_seq:
        v = ord(c)
        if v < 33 or v > 126:
            raise ValueError("Invalid Sanger ASCII value: " + v + ". Should be in range [33,126].")
        phred.append(v - 33)
    return phred

cpdef list sanger_ascii_to_probability(str ascii_seq):
    """
    Converts an ASCII string from a SAM file or FASTQ file with Sanger Phred format
    to a list of probabilities
    """
    cdef list probs = []
    cdef int phred
    for phred in sanger_ascii_to_phred(ascii_seq):
        probs.append(phred)
    return probs

cpdef list illumina_1_3_ascii_to_phred(str ascii_seq):
    """
    Converts an ASCII string from a Illumina 1.3-1.7 Phred format
    to a list of phred scores.
    """
    # Ascii range [64,126], phred range [0,62].
    cdef list phred = []
    cdef char c
    cdef int v
    for c in ascii_seq:
        v = ord(c)
        if v < 64 or v > 126:
            raise ValueError("Invalid Illumina 1.3 ASCII value: " + v + ". Should be in range [64,126].")
        phred.append(v - 64)
    return phred

cpdef list illumina_1_3_ascii_to_probability(str ascii_seq):
    """
    Converts an ASCII string from a Illumina 1.3-1.7 Phred format
    to a list of probabilities
    """
    cdef list probs = []
    cdef int phred
    for phred in illumina_1_3_ascii_to_phred(ascii_seq):
        probs.append(phred)
    return probs

cpdef list illumina_1_8_casava_ascii_to_phred(str ascii_seq):
    """
    Converts an ASCII string from a Illumina 1.8-x.x Phred format
    to a list of phred scores.
    """
    return sanger_ascii_to_phred(ascii_seq)

cpdef list illumina_1_8_casava_ascii_to_probability(str ascii_seq):
    """
    Converts an ASCII string from a Illumina 1.8-x.x Phred format
    to a list of probabilities
    """
    return sanger_ascii_to_probability(ascii_seq)

cpdef list solid_ascii_to_phred(str ascii_seq):
    """
    Converts an ASCII string from a SOLiD Phred format
    to a list of phred scores.
    """
    return sanger_ascii_to_phred(ascii_seq)

cpdef list solid_ascii_to_probability(str ascii_seq):
    """
    Converts an ASCII string from a SOLiD Phred format
    to a list of probabilities
    """
    return sanger_ascii_to_probability(ascii_seq)


