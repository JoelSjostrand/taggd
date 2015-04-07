"""
Interface for writing/reading FASTAQ and SAM files
"""

import pysam as ps
import sys
import os
import taggd.io.fastq_utils as fu
from taggd.io.sam_record import *
from taggd.io.fasta_record import *
from taggd.io.fastq_record import *

# Global variables for match types.
cdef int FASTQ, FASTA, SAM, BAM

class ReadsReaderWriter():
    """
    Provides opening of a reads file, and writing to the same format to one or more
    specified files.
    """

    def __init__(self, str reads_infile_name):
        """Constructor"""

        self.file_type = -1
        self.infile_name = reads_infile_name
        self.infile = None
        self.infile_header = None

        # Supported file types
        global FASTQ
        global FASTA
        global SAM
        global BAM
        FASTQ, FASTA, SAM, BAM = range(4)

        cdef str suffix = self.infile_name.split(".")[-1].lower()
        if suffix == "fa" or suffix == "fasta":
            self.file_type = FASTA
        elif suffix == "fq" or suffix == "fastq":
            self.file_type = FASTQ
        elif suffix == "sam":
            self.file_type = SAM
        elif suffix == "bam":
            self.file_type = BAM
        else:
            raise ValueError("Unsupported reads file format!")

        # Read header.
        if self.file_type == SAM:
            self.infile = ps.AlignmentFile(self.infile_name, "r", check_header=True, check_sq=False)
            self.infile_header = self.infile.header
            self.infile.close()
        elif self.file_type == BAM:
            self.infile = ps.AlignmentFile(self.infile_name, "r", check_header=True, check_sq=False)
            self.infile_header = self.infile.header
            self.infile.close()

    def reader_open(self):
        """Opens the reads file using appropriate format."""

        # Open file.
        if self.file_type == FASTA or self.file_type == FASTQ:
            self.infile = fu.readfq(open(self.infile_name, "r"))
        elif self.file_type == SAM:
            self.infile = ps.AlignmentFile(self.infile_name, "r", check_header=True, check_sq=False)
        elif self.file_type == BAM:
            self.infile = ps.AlignmentFile(self.infile_name, "rb", check_header=True, check_sq=False)
        else:
            raise ValueError("Unsupported reads file format!")

        cdef object rec, last
        while True:
            last = None
            for orig in self.infile:
                if self.file_type == FASTA:
                    rec = FASTARecord(orig)
                    last = rec
                    yield rec
                    break
                elif self.file_type == FASTQ:
                    rec = FASTQRecord(orig)
                    last = rec
                    yield rec
                    break
                elif self.file_type == SAM or self.file_type == BAM:
                    rec = SAMRecord(orig)
                    last = rec
                    yield rec
                    break
                else:
                    raise ValueError("Unsupported reads file format!")
            if last == None:
                break

    def reader_close(self):
        """Closes the infile."""
        if self.infile != None:
            self.infile.close()
            self.infile = None

    def __exit__(self, type, value, tb):
        self.close_read()

    def get_format(self):
        """Returns the file format."""
        if self.file_type == FASTA:
            return "fa"
        if self.file_type == FASTQ:
            return "fq"
        if self.file_type == SAM:
            return "sam"
        if self.file_type == BAM:
            return "bam"
        return None

    def get_writer(self, str outfile_name):
        """
        Returns a writer.
        """
        if os.path.exists(outfile_name):
            os.remove(outfile_name)

        if self.file_type == FASTA or self.file_type == FASTQ:
            return open(outfile_name, "w")
        else:
            if self.infile_header == None:
                raise ValueError("Error: missing header in SAM/BAM file")
            if self.file_type == SAM:
                return ps.AlignmentFile(outfile_name, "wh", header=self.infile_header)
            if self.file_type == BAM:
                return ps.AlignmentFile(outfile_name, "wb", header=self.infile_header)
            else:
                return None

    def write_record(self, outfile, record):
        """Writes a record."""
        if self.file_type == FASTA:
            fu.writefa_record(outfile, record.unwrap())
        elif self.file_type == FASTQ:
            fu.writefq_record(outfile, record.unwrap())
        elif self.file_type == SAM:
            outfile.write(record.unwrap())
        elif self.file_type == BAM:
            outfile.write(record.unwrap())
        else:
            return