from cpython cimport bool

cdef int hamming_distance(str seq1, str seq2, int limit=?)
cdef int subglobal_distance(str, str)
cdef int levenshtein_distance(str seq1, str seq2, int limit=?)
