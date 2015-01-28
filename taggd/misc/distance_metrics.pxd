from cpython cimport bool

cdef unsigned int hamming_distance(str seq1, str seq2, unsigned int limit=?)

cdef list subglobal_distance(str, str)

cdef unsigned int levenshtein_distance(str seq1, str seq2, unsigned int limit=?)
