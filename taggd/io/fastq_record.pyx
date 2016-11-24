""" 
This class inherits from Record and it represents
a FASTQ record
"""
from taggd.io.record import Record

class FASTQRecord(Record):
    """
    Holds a FASTQ record. 
    Not cdef-cythonized so as to be threading compatible.
    """

    def __init__(self, object fqr):
        """
        Constructor. 
        Represents a (annotation, seq, qual) tuple as a Record.
        """
        Record.__init__(self)
        self.annotation = fqr[0]
        self.sequence = fqr[1]
        self.attributes["quality"] = fqr[2]
        self.taggdtags = ""

    def add_tags(self, list added):
        """
        Appends tags for extra information
        :param added a list of tag tuples (name,value)
        """
        cdef str k
        cdef object v
        self.taggdtags = ' '.join(["{}:{}".format(k,v) for k,v in added])

    def unwrap(self):
        """
        Returns (annotation, sequence, quality)
        """
        return ("{} {}".format(self.annotation, self.taggdtags), 
                self.sequence, self.attributes["quality"])

    def __str__(self):
        """
        String representation
        """
        cdef str fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
        return fq_format.format(header_comments=self.annotation, 
                                sequence=self.sequence, quality=self.attributes["quality"])