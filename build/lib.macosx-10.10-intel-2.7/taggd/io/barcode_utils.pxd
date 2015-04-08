cdef class Barcode:
    cdef public str sequence
    cdef public list attributes

cpdef dict read_barcode_file(str infile_path)

cpdef int estimate_min_edit_distance(dict true_barcodes, int max_iters)