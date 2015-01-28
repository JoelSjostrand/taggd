

cdef dict get_candidates(str read_barcode)

cdef list get_distances(str read_barcode, dict candidates)

cdef list get_top_hits(list qual_hits)