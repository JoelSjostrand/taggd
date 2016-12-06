""" 
A wrapper class to store some statistics
"""

class Statistics(object):
    """
    Shorthand for a statistics object. 
    Not cdef-cythonized so as to keep threading compatible.
    """
    def __init__(self, int id, int max_edit_distance):
        """
        Constructor
        """
        self.id = id
        self.time = 0
        self.total_reads = 0
        self.total_reads_wr = 0
        self.perfect_matches = 0
        self.imperfect_unambiguous_matches = 0
        self.imperfect_ambiguous_matches = 0   # Non-unique
        self.unmatched = 0
        self.edit_distance_counts = [0] * (max_edit_distance + 1)
      
    def __str__(self):
        """
        String representation
        """
        # TODO use proper string formatting for this
        cdef int mat = self.perfect_matches + self.imperfect_unambiguous_matches
        cdef int amb = self.total_reads - self.perfect_matches - self.imperfect_unambiguous_matches - self.unmatched
        return "# Total execution time in secs: " + str(self.time) +\
            "\n# Total reads: " + str(self.total_reads) +\
            "\n# Total reads written: " + str(self.total_reads_wr) +\
            "\n# Matches: " + str(mat) + "   [" + str(mat*100.0/self.total_reads) + "%]" +\
            "\n#   - Perfect matches: " + str(self.perfect_matches) + \
            "   [" + str(self.perfect_matches*100.0/self.total_reads) + "%]" +\
            "\n#   - Imperfect matches: " + str(self.imperfect_unambiguous_matches) + \
            "   [" + str(self.imperfect_unambiguous_matches*100.0/self.total_reads) + "%]" +\
            "\n# Ambiguous matches: " + str(amb) + "   [" + str(amb*100.0/self.total_reads) + "%]" +\
            "\n#   - Non-unique ambiguous matches: " + str(self.imperfect_ambiguous_matches) +\
            "\n# Unmatched: " + str(self.unmatched) + "   [" + str(self.unmatched*100.0/self.total_reads) + "%]" +\
            "\n# Matched edit distance counts for 0,1,...: " + str(self.edit_distance_counts)

    def __iadd__(self, other):
        """
        Overrided += operator
        """
        self.time += other.time
        self.total_reads += other.total_reads
        self.total_reads_wr += other.total_reads_wr
        self.perfect_matches += other.perfect_matches
        self.imperfect_unambiguous_matches += other.imperfect_unambiguous_matches
        self.imperfect_ambiguous_matches += other.imperfect_ambiguous_matches
        self.unmatched += other.unmatched
        self.edit_distance_counts = [x + y for x, y in zip(self.edit_distance_counts, 
                                                           other.edit_distance_counts)]
        return self
