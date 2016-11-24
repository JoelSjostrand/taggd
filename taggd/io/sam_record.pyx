"""
A simple class that is meant to mimic a Pysam record
and add some extra attributes so it can later be
written to a sam file
"""
import pysam
from taggd.io.record import Record

class SAMRecord(Record):
    """
    Holds a SAM record.
    Follows pysam's representation.
    """

    def __init__(self, object alseq):
        """
        Constructor. 
        Represents a pysam.AlignedSequence as a Record.
        """
        Record.__init__(self)
        self.annotation = str(alseq.query_name)
        self.sequence = str(alseq.query_sequence)
        self.attributes["flag"] = alseq.flag
        self.attributes["reference_id"] = alseq.reference_id
        self.attributes["reference_start"] = alseq.reference_start
        self.attributes["mapping_quality"] = alseq.mapping_quality
        self.attributes["cigar"] = alseq.cigar
        self.attributes["next_reference_id"] = alseq.next_reference_id
        self.attributes["reference_start"] = alseq.reference_start
        self.attributes["template_length"] = alseq.template_length
        self.attributes["query_qualities"] = alseq.query_qualities
        self.attributes["tags"] = alseq.tags
        self.attributes["taggdtags"] = None

    def add_tags(self, list added):
        """
        Appends extra tags for taggd.
        :param added a list of tag tuples (name,value)
        """
        self.attributes["taggdtags"] = added

    def unwrap(self):
        """
        Returns a pysam alignment.
        """
        cdef object a = pysam.AlignedSegment()
        a.query_name = self.annotation
        a.flag = self.attributes["flag"]
        a.reference_id = self.attributes["reference_id"]
        a.reference_start = self.attributes["reference_start"]
        a.mapping_quality = self.attributes["mapping_quality"]
        a.cigar = self.attributes["cigar"]
        a.next_reference_id = self.attributes["reference_id"]
        a.next_reference_start = self.attributes["reference_start"]
        a.template_length = self.attributes["template_length"]
        a.query_sequence = self.sequence
        a.query_qualities = self.attributes["query_qualities"]
        a.tags = self.attributes["tags"]
        if self.attributes["taggdtags"] != None:
                a.tags += self.attributes["taggdtags"]
        return a