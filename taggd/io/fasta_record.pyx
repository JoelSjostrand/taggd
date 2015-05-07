from taggd.io.record import Record

class FASTARecord(Record):
    """Holds a FASTA record. Not cdef-cythonized so as to be threading compatible."""

    def __init__(self, object fqr):
        """Constructor. Represents a (annotation, seq) tuple as a Record."""
        Record.__init__(self)
        self.annotation = fqr[0]
        self.sequence = fqr[1]
        self.taggdtags = ""

    def add_tags(self, list added):
        """Appends tags"""
        self.taggdtags = ""
        cdef str k
        cdef object v
        for k,v in added:
            self.taggdtags += " " + k + ":" + str(v)

    def unwrap(self):
        """
        Returns(annotation, sequence).
        """
        return (self.annotation + self.taggdtags, self.sequence)


    def __str__(self):
        """String representation"""
        cdef str fa_format = '>{header_comments}\n{sequence}\n'
        return fa_format.format(header_comments=self.annotation, sequence=self.sequence)


