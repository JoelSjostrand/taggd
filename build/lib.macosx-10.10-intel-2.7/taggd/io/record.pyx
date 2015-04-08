
class Record(object):
    """
    Shorthand for a record. To keep things simple, we use a single object
    regardless of file format, and let util classes deal with conversions.
    """

    def __init__(self):
        self.annotation = None
        self.sequence = None
        self.attributes = dict()

    def unwrap(self):
        """Unwraps an object to its proper original representation."""
        return None

    def add_tags(self, list added):
        """Adds tags provided as a list of (name, attr) pairs."""
        return

    def __str__(self):
        return str(self.annotation) + "\t" + str(self.sequence) + "\t" + str(self.attributes)