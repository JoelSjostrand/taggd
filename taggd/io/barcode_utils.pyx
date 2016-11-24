# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
"""
Contains utilities for handling barcodes.
"""

from taggd.misc.distance_metrics cimport *
import random

cdef class Barcode:
    """
    Holds a Spatial Transcriptomics barcode sequence and 
    attributes such as feature coordinates.
    """
    def __cinit__(self, str seq, list attributes):
        self.sequence = seq
        self.attributes = attributes

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
    cdef object infile
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
    :param a dict of barcode -> Barcode
    :param max_iters the max number of iterations to try
    """
    # Compute minimum edit distance.
    # Get the barcodes and shuffle them
    cdef list seqslist = true_barcodes.keys()
    random.shuffle(seqslist)
    # Iterate the barcodes to see if the comply the min distance requirement
    cdef int min_dist = 1000000000
    cdef str barcode1
    cdef str barcode2
    cdef int i
    cdef int dist
    cdef int iter = 0
    for i,barcode1 in enumerate(seqslist):
        for barcode2 in seqslist[(i+1):]:
            dist = hamming_distance(barcode1, barcode2, min_dist)
            if dist < min_dist: min_dist = dist
            iter += 1
            if iter >= max_iters:
                return min_dist
    return min_dist