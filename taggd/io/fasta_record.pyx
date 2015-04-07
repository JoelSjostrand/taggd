from taggd.io.record import Record

class FASTARecord(Record):

    def __init__(self, object fqr):
        """Constructor. Represents a (annotation, seq) tuple as a Record."""
        Record.__init__(self)
        self.annotation = fqr[0]
        self.sequence = fqr[1]

    def add_tags(self, list added):
        #TODO this can be slow, use ''.{} format instead
        for k,v in added:
            self.annotation += " " + k + ":" + str(v)

    def unwrap(self):
        return (self.annotation, self.sequence)


    def __str__(self):
        fa_format = '>{header_comments}\n{sequence}\n'
        return fa_format.format(header_comments=self.annotation, sequence=self.sequence)


