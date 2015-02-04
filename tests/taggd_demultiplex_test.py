#! /usr/bin/env python
"""
Unit-test for run-tests
"""

import unittest
import tempfile
import os
import taggd.core.demultiplex as deplex


class TestDemultiplexer(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        # Obtain paths and files.
        testdir = str(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
        self.inbarcodes = os.path.join(testdir, "testbarcodes.tsv")
        self.insam = os.path.join(testdir, "testset.sam")
        self.inbam = os.path.join(testdir, "testset.bam")
        self.infq = os.path.join(testdir, "testset.fq")
        self.infa = os.path.join(testdir, "testset.fa")
        assert (os.path.exists(self.inbarcodes))
        assert (os.path.exists(self.insam))
        assert (os.path.exists(self.inbam))
        assert (os.path.exists(self.infq))
        assert (os.path.exists(self.infa))

        # Output directories
        self.outdir = tempfile.mkdtemp(prefix="taggd_demultiplex_test_out")
        print "# Demultiplexer test output directory: " + self.outdir


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
        args = ["--k", "7", "--max_edit_distance", "7", "--overhang", "2"]
        args += [self.inbarcodes, self.insam, os.path.join(self.outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running SAM test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal SAM test", "sam")
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal SAM test failed\n")


    def test_normal_bam_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """
        args = ["--k", "6", "--max_edit_distance", "5", "--overhang", "0"]
        args += [self.inbarcodes, self.inbam, os.path.join(self.outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running BAM test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal BAM test", "bam")
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal BAM test failed\n")


    def test_normal_fq_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """
        args = ["--k", "8", "--max_edit_distance", "8", "--overhang", "10"]
        args += [self.inbarcodes, self.infq, os.path.join(self.outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running Fastq test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal Fastq test", "fq")
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal Fastq test failed\n")


    def test_normal_fa_run(self):
        """
        Tests taggd demultiplexer on a variety of small files.
        """
        args = ["--k", "4", "--max_edit_distance", "8", "--overhang", "3"]
        args += [self.inbarcodes, self.infa, os.path.join(self.outdir, "outfile")]

        # Start the demultiplexer
        try:
            print "\n# Running Fasta test with parameters: " + " ".join(args)
            deplex.main(args)
            self.validate_output_data("Normal Fasta test", "fa")
        except Exception as e:
            print e
            self.assertTrue(0, "Running Normal Fasta test failed\n")


    def validate_output_data(self, expName, suffix):
        print "# Validating test " + expName

        # Verify existence of output files and temp files
        self.assertNotEqual(os.listdir(self.outdir), [], "Output folder is not empty")

        match = os.path.join(self.outdir, "outfile.matched." + suffix)
        unmatch = os.path.join(self.outdir, "outfile.unmatched." + suffix)
        ambig = os.path.join(self.outdir, "outfile.ambiguous." + suffix)
        res = os.path.join(self.outdir, "outfile.results.tsv")

        self.assertTrue(os.path.exists(match), "Matched file exists")
        self.assertTrue(os.path.exists(unmatch), "Unmatched file exists")
        self.assertTrue(os.path.exists(ambig), "Ambiguous file exists")
        self.assertTrue(os.path.exists(res), "Results file exists")

        #self.assertTrue(os.path.getsize(match) > 16, "Match file is not empty")
        #self.assertTrue(os.path.getsize(unmatch) > 16, "Unmatched file is not empty")
        #self.assertTrue(os.path.getsize(ambig) > 16, "Ambiguous file is not empty")
        #self.assertTrue(os.path.getsize(res) > 16, "Results file is not empty")


if __name__ == '__main__':
    unittest.main()
