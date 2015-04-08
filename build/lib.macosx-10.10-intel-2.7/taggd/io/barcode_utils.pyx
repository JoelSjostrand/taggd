"""
Contains utilities for handling barcodes.
"""

from taggd.misc.distance_metrics cimport *
import random

cdef class Barcode:
    """
    Holds a Spatial Transcriptomics barcode sequence and feature coordinates.
    """

    def __cinit__(self, str seq_, list attributes_):
        self.sequence = seq_
        self.attributes = attributes_


cpdef dict read_barcode_file(str infile_path):
    """
    Reads a barcode ("chip") file and returns the barcodes.
    :param infile_path: chip file.
    :return: a dictionary with barcode sequences as keys and Barcode instances as values.
    """
    cdef dict res_dict = {}
    cdef str line
    cdef list tmp
    cdef str seq
    cdef int length = 0
    with open(infile_path) as infile:
        for line in infile:
            tmp = line.strip().split()
            seq = tmp[0]
            if length == 0:
                length = len(seq)
            if len(seq) != length:
                raise ValueError("Barcode file incorrect, varying lengths among barcodes")
            if seq in res_dict:
                raise ValueError("Barcode file incorrect, duplicate barcode: " + seq)
            res_dict[tmp[0]] = Barcode(tmp[0], tmp[1:])
    return res_dict

cpdef int estimate_min_edit_distance(dict true_barcodes, int max_iters):
    """
    Reads a barcodes dict and estimates the minimum edit distance
    by comparing a certain number of pairs.
    """

    cdef int length = 0
    cdef int min_dist = 1000000000
    cdef str ln
    cdef list tmp
    cdef str seq

    # Compute minimum edit distance.
    cdef list seqslist = true_barcodes.keys()
    random.shuffle(seqslist)
    cdef int i
    cdef int j
    cdef int dist
    cdef int iter = 0
    for i in xrange(len(seqslist)):
        for j in xrange(int(i+1), len(seqslist)):
            dist = hamming_distance(seqslist[i], seqslist[j], min_dist)
            if dist < min_dist:
                min_dist = dist
            iter += 1
            if iter >= max_iters:
                return min_dist
            #print str(i) + " vs " + str(j) + ": " + str(dist)
    return min_dist