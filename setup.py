from setuptools import setup, find_packages
from setuptools.extension import Extension
from glob import glob
import numpy
from Cython.Distutils import build_ext

ext_modules = [
    Extension("taggd.core.demultiplex_core_functions",   ["taggd/core/demultiplex_core_functions.pyx"]),
    Extension("taggd.core.demultiplex_sub_functions",    ["taggd/core/demultiplex_sub_functions.pyx"]),
    Extension("taggd.core.demultiplex_search_functions", ["taggd/core/demultiplex_search_functions.pyx"]),
    Extension("taggd.core.match",                        ["taggd/core/match.pyx"]),
    Extension("taggd.core.match_type",                   ["taggd/core/match_type.pyx"]),
    Extension("taggd.core.statistics",                   ["taggd/core/statistics.pyx"]),
    Extension("taggd.misc.distance_metrics",             ["taggd/misc/distance_metrics.pyx"]),
    Extension("taggd.misc.kmer_utils",                   ["taggd/misc/kmer_utils.pyx"]),
    Extension("taggd.io.fastq_utils",                    ["taggd/io/fastq_utils.pyx"]),
    Extension("taggd.io.barcode_utils",                  ["taggd/io/barcode_utils.pyx"]),
    Extension("taggd.io.record",                         ["taggd/io/record.pyx"]),
    Extension("taggd.io.sam_record",                     ["taggd/io/sam_record.pyx"]),
    Extension("taggd.io.fasta_record",                   ["taggd/io/fasta_record.pyx"]),
    Extension("taggd.io.fastq_record",                   ["taggd/io/fastq_record.pyx"]),
    Extension("taggd.io.reads_reader_writer",            ["taggd/io/reads_reader_writer.pyx"])
]
cmdclass = { 'build_ext': build_ext }

setup(
	name = "taggd",
	version = '0.3.5',
	author = 'Joel Sjostrand, Jose Fernandez',
	author_email = 'joel.sjostrand@scilifelab.se, jose.fernandez.navarro@scilifelab.se',
	license = 'Open BSD',
    description = 'Bioinformatics genetic barcode demultiplexing',
    url = 'https://github.com/SpatialTranscriptomicsResearch/taggd',
    download_url = 'https://github.com/SpatialTranscriptomicsResearch/taggd/0.3.2',
	scripts = glob("scripts/*.py"),
    packages = ['taggd', 'taggd.core', 'taggd.io', 'taggd.misc'],
    package_data = {'': ['*.pyx', '*.pxd', '*.h', '*.c'], },
    setup_requires = ["cython"],
	install_requires = [
	    'setuptools',
	    'pysam',
	   	'numpy',
	],
	test_suite = 'tests',
	cmdclass = cmdclass,
    ext_modules = ext_modules,
	include_dirs = [numpy.get_include(), '.'],
    keywords = ['bioinformatics', 'demultiplexing']
    )