from cpython cimport bool

cdef list get_candidates(str read_barcode)

cdef list get_distances(str read_barcode, list candidates)

cdef list get_top_hits(list qual_hits)