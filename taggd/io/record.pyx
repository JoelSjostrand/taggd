"""
Base classes for different records that will be written to files
"""
class Record(object):
    """
    Shorthand for a record. To keep things simple, we use a single object
    regardless of file format, and let util classes deal with conversions.
    Not cdef-cythonized to keep thread compatible.
    """

    def __init__(self):
        """
        Constructor
        """
        self.annotation = None
        self.sequence = None
        self.attributes = dict()

    def unwrap(self):
        """
        Unwraps an object to its proper original representation.
        """
        return None

    def add_tags(self, list added):
        """
        Adds tags provided as a list of (name, attr) pairs.
        """
        return

    def __str__(self):
        """
        String representation
        """
        return "%s\t%s\t%s" % (str(self.annotation),str(self.sequence),str(self.attributes))