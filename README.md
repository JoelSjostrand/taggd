TagGD barcode demultiplexing utilities for Spatial Transcriptomics data.

This is the Python version, which is a generalized and more
up-to-date version of the C++ demultiplexer named "findIndexes"
which you can find here https://github.com/pelinakan/UBD.

See http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0057521
for the peer-reviewed reference to the program.

The main idea is to extract a cDNA barcodes
from the input file (FASTQ, FASTA or BAM) and then
tries to find a match in a file with a list
of reference barcodes using a kmer-based approach.
All the matches will be outputted and the barcode
and spatial information will be added.

##Manual

See Wiki section

##Requirements:

TagGD requires PySam and Numpy.

##Installation:

(Assuming that you have a virtual environment
installed such as Anaconda 2.7)

cd <taggd demultiplexer root>
python setup.py build
python setup.py install

##Run:

taggd_demultiplex.py -h


##Contact: 

{joel.sjostrand, jose.fernandez.navarro}@scilifelab.se.


[![Build Status](https://travis-ci.org/JoelSjostrand/taggd.svg?branch=master)](https://travis-ci.org/JoelSjostrand/taggd)
