#!/usr/bin/env python
"""
Starts and runs a TagGD demultiplexing job. Type -h for help.
\n
See https://github.com/SpatialTranscriptomicsResearch/taggd/wiki for manual and the most recent source code.
\n
See http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0057521
for the peer-reviewed reference to the program.
"""


import taggd.core.demultiplex as deplex

if __name__ == "__main__":
    deplex.main()


