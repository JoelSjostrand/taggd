from setuptools import setup, find_packages
from setuptools.extension import Extension
from glob import glob
import numpy
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = { }
ext_modules = [ ]


if use_cython:
    ext_modules += [
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
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
    Extension("taggd.core.demultiplex_core_functions",   ["taggd/core/demultiplex_core_functions.c"]),
    Extension("taggd.core.demultiplex_sub_functions",    ["taggd/core/demultiplex_sub_functions.c"]),
    Extension("taggd.core.demultiplex_search_functions", ["taggd/core/demultiplex_search_functions.c"]),
    Extension("taggd.core.match",                        ["taggd/core/match.c"]),
    Extension("taggd.core.match_type",                   ["taggd/core/match_type.c"]),
    Extension("taggd.core.statistics",                   ["taggd/core/statistics.c"]),
    Extension("taggd.misc.distance_metrics",             ["taggd/misc/distance_metrics.c"]),
    Extension("taggd.misc.kmer_utils",                   ["taggd/misc/kmer_utils.c"]),
    Extension("taggd.io.fastq_utils",                    ["taggd/io/fastq_utils.c"]),
    Extension("taggd.io.barcode_utils",                  ["taggd/io/barcode_utils.c"]),
    Extension("taggd.io.record",                         ["taggd/io/record.c"]),
    Extension("taggd.io.sam_record",                     ["taggd/io/sam_record.c"]),
    Extension("taggd.io.fasta_record",                   ["taggd/io/fasta_record.c"]),
    Extension("taggd.io.fastq_record",                   ["taggd/io/fastq_record.c"]),
    Extension("taggd.io.reads_reader_writer",            ["taggd/io/reads_reader_writer.c"])
    ]


setup(
	name = "taggd",
	version = '0.3.1',
	author = 'Joel Sjostrand',
	author_email = 'joel.sjostrand@scilifelab.se, jose.fernandez.navarro@scilifelab.se',
	license = 'Open BSD',
    description = 'Bioinformatics genetic barcode demultiplexing',
    url = 'https://github.com/SpatialTranscriptomicsResearch/taggd',
    download_url = 'https://github.com/SpatialTranscriptomicsResearch/taggd/0.3.1',
	scripts = glob("scripts/*.py"),
    packages = ['taggd', 'taggd.core', 'taggd.io', 'taggd.misc'],
    package_data = {'': ['*.pyx', '*.pxd', '*.h', '*.c'], },
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