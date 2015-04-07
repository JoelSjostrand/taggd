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
    return "#Annotation\tMatch_result\tBarcode\tEdit_distance\tAmbiguous_top_hits\t" \
        "Qualified_candidates\tRaw_candidates\tLast_position\tApprox_insertions\tApprox_deletions"

class Match(object):
    """
    Shorthand for a match.
    """
    def __init__(self,
                 object record_,
                 int match_type_,
                 str barcode_,
                 int edit_distance_,
                 int ambiguous_top_hits_,
                 int qualified_candidates_,
                 int raw_candidates_,
                 int last_position_,
                 int insertions_read_,
                 int insertions_barcode_):
        """
        Constructor.
        """
        self.record = record_
        self.match_type = match_type_
        self.barcode = barcode_
        self.edit_distance = edit_distance_
        self.ambiguous_top_hits = ambiguous_top_hits_
        self.qualified_candidates = qualified_candidates_
        self.raw_candidates = raw_candidates_
        self.last_position = last_position_
        self.insertions_read = insertions_read_
        self.insertions_barcode = insertions_barcode_


    def __str__(self):
        """
        String converter.
        TODO this is slow, python recommends to build strings
        in the form ''.{} ....
        """
        return (self.record.annotation + '\t' +
                match_type.match_type_to_str(self.match_type) + '\t' +
                self.barcode + '\t' +
                str(self.edit_distance) + '\t' +
                str(self.ambiguous_top_hits) + '\t' +
                str(self.qualified_candidates) + '\t' +
                str(self.raw_candidates) + '\t' +
                str(self.last_position) + '\t' +
                str(self.insertions_read) + '\t' +
                str(self.insertions_barcode))