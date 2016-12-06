from cpython cimport bool

cdef object get_kmers_dicts(list seqs, int k, bool round_robin=?, int slider_increment=?)
cdef list get_kmers(str seq, int k, bool round_robin=?, int slider_increment=?)
