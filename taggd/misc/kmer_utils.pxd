from cpython cimport bool

cdef list get_kmers_dicts(list seqs, unsigned int k, bool round_robin=?, unsigned int slider_increment=?)

cdef list get_kmers(str seq, unsigned int k, bool round_robin=?, unsigned int slider_increment=?)
