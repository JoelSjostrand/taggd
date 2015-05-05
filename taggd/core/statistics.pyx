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
        return "# Total execution time in secs: " + str(self.time) +\
            "\n# Total reads: " + str(self.total_reads) +\
            "\n# Total reads written: " + str(self.total_reads_wr) +\
            "\n# Matches: " + str(self.perfect_matches + self.imperfect_unambiguous_matches) +\
            "\n#   - Perfect matches: " + str(self.perfect_matches) +\
            "\n#   - Imperfect matches: " + str(self.imperfect_unambiguous_matches) +\
            "\n# Ambiguous matches: " + str(self.total_reads - self.perfect_matches - self.imperfect_unambiguous_matches - self.unmatched) +\
            "\n#   - Non-unique ambiguous matches: " + str(self.imperfect_ambiguous_matches) +\
            "\n# Unmatched: " + str(self.unmatched) +\
            "\n# Matched edit distance counts for 0,1,...+ " + str(self.edit_distance_counts)

    def __iadd__(self, other):
        self.time += other.time
        self.total_reads += other.total_reads
        self.total_reads_wr += other.total_reads_wr
        self.perfect_matches += other.perfect_matches
        self.imperfect_unambiguous_matches += other.imperfect_unambiguous_matches
        self.imperfect_ambiguous_matches += other.imperfect_ambiguous_matches
        self.unmatched += other.unmatched
        cdef int i
        for i in xrange(len(self.edit_distance_counts)):
            self.edit_distance_counts[i] += other.edit_distance_counts[i]
        return self
