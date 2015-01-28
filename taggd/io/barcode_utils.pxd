cdef class Barcode:
    cdef public str sequence
    cdef public list attributes

cpdef dict read_barcode_file(str infile_path)

cpdef unsigned int estimate_min_edit_distance(dict true_barcodes, unsigned int max_iters)