cdef int KILL = -1
cdef int UNMATCHED = 0
cdef int MATCHED_PERFECTLY = 1
cdef int MATCHED_UNAMBIGUOUSLY = 2
cdef int MATCHED_AMBIGUOUSLY = 3

cdef str match_type_to_str(int match_type)