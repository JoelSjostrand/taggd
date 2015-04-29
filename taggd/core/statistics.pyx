class Statistics(object):
    """
    Shorthand for a statistics object. Not cdef-cythonized so as to keep threading compatible.
    """


    def __init__(self, int id, int max_edit_distance):
        self.id = id
        self.time = 0
        self.total_reads = 0
        self.total_reads_wr = 0
        self.perfect_matches = 0
        self.imperfect_unambiguous_matches = 0
        self.imperfect_ambiguous_matches = 0   # Non-unique
        self.unmatched = 0
        self.edit_distance_counts = list()
        cdef int i
        for i in xrange(max_edit_distance + 1):
            self.edit_distance_counts.append(0)

    def __str__(self):
        print self.id
