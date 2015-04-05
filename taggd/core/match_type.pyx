"""
Defines different kinds of matches.
"""
cdef int KILL = -1
cdef int UNMATCHED = 0
cdef int MATCHED_PERFECTLY = 1
cdef int MATCHED_UNAMBIGUOUSLY = 2
cdef int MATCHED_AMBIGUOUSLY = 3

cdef str match_type_to_str(int match_type):
    """
    Simple converter from match type to string.
    """
    if match_type == 0:
        return "UNMATCHED"
    elif match_type == 1:
        return "MATCHED_PERFECTLY"
    elif match_type == 2:
        return "MATCHED_UNAMBIGUOUSLY"
    elif match_type == 3:
        return "MATCHED_AMBIGUOUSLY"
    return -1