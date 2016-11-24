from cpython cimport bool

cdef object get_kmers_dicts(list seqs, int k, bool round_robin=?)
cdef list get_kmers(str seq, int k, bool round_robin=?)
