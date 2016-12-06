"""
This represent a barcode match with all
the necessary information that can be used to
write to a file later
"""
cimport taggd.core.match_type as match_type

cdef str get_match_header():
    """
    Returns a header for the match-to-string conversion.
    """
    return "#Annotation\tMatch_result\tBarcode\tEdit_distance"

class Match(object):
    """
    Shorthand for a match. 
    Not cdef-cythonized so as to keep threading compatible.
    """

    def __init__(self,
                 object record,
                 int match_type,
                 str barcode,
                 int edit_distance = 1):
        """
        Constructor.
        """
        self.record = record
        self.match_type = match_type
        self.barcode = barcode
        self.edit_distance = edit_distance


    def __str__(self):
        """
        String converter.
        """
        return ("%s\t%s\t%s\t%i" % (self.record.annotation, 
                                    match_type.match_type_to_str(self.match_type), 
                                    self.barcode, self.edit_distance))