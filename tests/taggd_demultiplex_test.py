#! /usr/bin/env python
"""
Unit-test for run-tests
"""

import unittest
import tempfile
import os
import taggd.core.demultiplex as deplex
import filecmp
import time

class TestDemultiplexer(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        # Obtain paths and files.
        self.testdir = str(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
        self.inbarcodes = os.path.join(self.testdir, "testbarcodes.tsv")
        self.insam = os.path.join(self.testdir, "testset.sam")
        self.inbam = os.path.join(self.testdir, "testset.bam")
        self.infq = os.path.join(self.testdir, "testset.fq")
        self.infa = os.path.join(self.testdir, "testset.fa")
        assert (os.path.exists(self.inbarcodes))
        assert (os.path.exists(self.insam))
        assert (os.path.exists(self.inbam))
        assert (os.path.exists(self.infq))
        assert (os.path.exists(self.infa))



    #@classmethod
    #def tearDownClass(self):
        # outcnt = os.listdir(self.outdir)
        # for cnt in outcnt:
        # 	os.remove(os.path.join(self.outdir, cnt))
        # os.removedirs(self.outdir)


    def test_normal_sam_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """

        outdir = tempfile.mkdtemp(prefix="taggd_demultiplex_test_out_sam_")
        print "# Demultiplexer test output directory: " + outdir

        args = ["--k", "7", "--max-edit-distance", "7", "--overhang", "2", "--subprocesses", "3", "--seed", "dsfiogwhgfsaeadsgfADSgsagaagd"]
        args += [self.inbarcodes, self.insam, os.path.join(outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running SAM test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal SAM test", "sam", outdir)
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal SAM test failed\n")


    def test_normal_bam_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """
        outdir = tempfile.mkdtemp(prefix="taggd_demultiplex_test_out_bam_")
        args = ["--k", "6", "--max-edit-distance", "5", "--overhang", "0", "--subprocesses", "3", "--seed", "dsfiogwhgfsaeadsgfADSgsagaagd"]
        args += [self.inbarcodes, self.inbam, os.path.join(outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running BAM test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal BAM test", "bam", outdir)
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal BAM test failed\n")


    def test_normal_fq_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """
        outdir = tempfile.mkdtemp(prefix="taggd_demultiplex_test_out_fastq_")
        args = ["--k", "4", "--max-edit-distance", "8", "--overhang", "3", "--subprocesses", "3", "--seed", "dsfiogwhgfsaeadsgfADSgsagaagd"]
        args += [self.inbarcodes, self.infq, os.path.join(outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running Fastq test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal Fastq test", "fq", outdir)
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal Fastq test failed\n")


    def test_normal_fa_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """
        outdir = tempfile.mkdtemp(prefix="taggd_demultiplex_test_out_fasta_")
        args = ["--k", "4", "--max-edit-distance", "8", "--overhang", "3", "--subprocesses", "3", "--seed", "dsfiogwhgfsaeadsgfADSgsagaagd"]
        args += [self.inbarcodes, self.infa, os.path.join(outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running Fasta test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal Fasta test", "fa", outdir)
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal Fasta test failed\n")


    def compare_to_expected_results(self, suffix, filename, file_description, outdir):
         filepath_from_test = os.path.join(outdir, filename)
         expected_result_dir =  os.path.join(self.testdir, "expected_results")
         # The value of the variable "suffix" is used here as a name for a subdirectory
         expected_result_dir2 =  os.path.join(expected_result_dir, suffix)
         expected_result_filepath = os.path.join(expected_result_dir2, filename)
         self.assertTrue(os.path.exists(expected_result_filepath), expected_result_filepath + " exists")
         self.assertTrue(os.path.exists(filepath_from_test), file_description + " exists")
         comp = filecmp.cmp(filepath_from_test, expected_result_filepath, shallow=False)
         self.assertTrue(comp, file_description + " has the expected contents")

    def validate_output_data(self, expName, suffix, outdir):
        print "# Validating test " + expName

        # Verify existence of output files and temp files
        self.assertNotEqual(os.listdir(outdir), [], "Output folder is not empty")

        time.sleep(2)
        self.compare_to_expected_results(suffix, "outfile_matched." + suffix, "Matched file", outdir)
        self.compare_to_expected_results(suffix, "outfile_unmatched." + suffix, "Unmatched file", outdir)
        self.compare_to_expected_results(suffix, "outfile_ambiguous." + suffix, "Ambiguous file", outdir)
        self.compare_to_expected_results(suffix, "outfile_results.tsv", "Results file", outdir)


if __name__ == '__main__':
    unittest.main()
