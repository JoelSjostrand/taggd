"""
Interface for writing/reading FASTQ/FASTA and SAM files
The idea is to write always the abstract class called Record
"""
import pysam as ps
import sys
import os
from cpython cimport bool
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
        """
        Constructor
        """

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

        # Extract the file type
        cdef str suffix = os.path.splitext(self.infile_name)[1].lower()
        if suffix == ".fa" or suffix == ".fasta":
            self.file_type = FASTA
        elif suffix == ".fq" or suffix == ".fastq":
            self.file_type = FASTQ
        elif suffix == ".sam":
            self.file_type = SAM
        elif suffix == ".bam":
            self.file_type = BAM
        else:
            raise ValueError("Unsupported reads file format!")

        # Read header.
        if self.file_type == SAM or self.file_type == BAM:
            self.infile = ps.AlignmentFile(self.infile_name, "r", 
                                           check_header=True, check_sq=False)
            self.infile_header = self.infile.header
            self.infile.close()

    def reader_open(self):
        """
        Opens the reads file using appropriate format.
        """
        # Ensure to close
        self.reader_close()
        
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

        if self.file_type == FASTA:
            for orig in self.infile:
                rec = FASTARecord(orig)
                yield rec
        elif self.file_type == FASTQ:
            for orig in self.infile:
                rec = FASTQRecord(orig)
                yield rec
        else:
            for orig in self.infile:
                rec = SAMRecord(orig)
                yield rec

    def reader_close(self):
        """
        Closes the input file handler.
        """
        if self.infile != None:
            self.infile.close()
            self.infile = None

    def __exit__(self, type, value, tb):
        """
        Always close the input file.
        """
        self.reader_close()

    def get_format(self):
        """
        Returns the file format.
        """
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
        :param outfile_name the name of the file to create
        :return the file handler so records can be written on it.
        """
        # Remove if exists
        if os.path.exists(outfile_name):
            os.remove(outfile_name)
        # Open the file and returns the handler
        if self.file_type == FASTA or self.file_type == FASTQ:
            return open(outfile_name, "w")
        elif self.file_type == SAM:
            return ps.AlignmentFile(outfile_name, "wh", header=self.infile_header)
        elif self.file_type == BAM:
            return ps.AlignmentFile(outfile_name, "wb", header=self.infile_header)
        else:
            raise ValueError("Unknown file format for writer")

    def write_record(self, outfile, record):
        """
        Writes a record in the filename descriptor given.
        Important out_handler must be a descriptor opened (using get_writer())
        Record's type should be the same as the one this instance was created from
        :param out_handler the outpuf file handler
        :param record the Record object to write
        """
        #TODO record could not be the same type of the file handler(check this)
        if self.file_type == FASTA:
            fu.writefa_record(outfile, record.unwrap())
        elif self.file_type == FASTQ:
            fu.writefq_record(outfile, record.unwrap())
        elif self.file_type == SAM:
            outfile.write(record.unwrap())
        elif self.file_type == BAM:
            outfile.write(record.unwrap())
        else:
            raise ValueError("Unknown file format for record")