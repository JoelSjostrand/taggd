from taggd.io.record import Record

class FASTQRecord(Record):
    """Holds a FASTQ record. Not cdef-cythonized so as to be threading compatible."""

    def __init__(self, object fqr):
        """Constructor. Represents a (annotation, seq, qual) tuple as a Record."""
        Record.__init__(self)
        self.annotation = fqr[0]
        self.sequence = fqr[1]
        self.attributes["quality"] = fqr[2]

    def add_tags(self, list added):
        cdef str k
        cdef object v
        for k,v in added:
            self.annotation += " " + k + ":" + str(v)

    def unwrap(self):
        """
        :return: (annotation, sequence, quality)
        """
        return (self.annotation, self.sequence, self.attributes["quality"])

    def __str__(self):
        cdef str fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
        return fq_format.format(header_comments=self.annotation, 
                                sequence=self.sequence, quality=self.attributes["quality"])