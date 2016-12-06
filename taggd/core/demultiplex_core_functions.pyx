"""
Main functions for carrying out the demultiplexing with
multithreading.
"""

import os
import fileinput
import time
import tempfile
import uuid
import pysam as ps
import multiprocessing as mp
from collections import deque
import taggd.core.statistics as statistics
import taggd.core.demultiplex_sub_functions as sub
from cpython cimport bool


def demultiplex(str filename_reads,
                str filename_matched,
                str filename_ambig,
                str filename_unmatched,
                str filename_results,
                int subprocesses):
    """
    Demultiplexes the contents of a reads file by dividing the work into
    parallel subprocesses, each writing its own file, then
    merging them together to coherent files.
    """
    cdef object pool = mp.Pool(subprocesses)
    
    # Fire off workers
    cdef object jobs = deque()
    cdef object job = None
    cdef int i
    cdef str fm = None, fa = None, fu = None, fr = None
    cdef list fn_matched = list()
    cdef list fn_ambig = list()
    cdef list fn_unmatched = list()
    cdef list fn_res = list()
    cdef str frmt = filename_reads.split(".")[-1].lower()
    cdef str tmp = tempfile.gettempdir()
    if tmp == "" or tmp == None: tmp = "."
    for i in xrange(subprocesses):
        if filename_matched != None:
            fm = "{}{}_part{}_{}.{}".format(tmp, 
                                            os.path.join(os.sep, os.path.basename(filename_matched)), 
                                            i, uuid.uuid4(), frmt)
            fn_matched.append(fm)
        if filename_ambig != None:
            fa = "{}{}_part{}_{}.{}".format(tmp, 
                                            os.path.join(os.sep, os.path.basename(filename_ambig)), 
                                            i, uuid.uuid4(), frmt)
            fn_ambig.append(fa)
        if filename_unmatched != None:
            fu = "{}{}_part{}_{}.{}".format(tmp, 
                                            os.path.join(os.sep, os.path.basename(filename_unmatched)), 
                                            i, uuid.uuid4(), frmt)
            fn_unmatched.append(fu)
        if filename_results != None:
            fr = "{}{}_part{}_{}.{}".format(tmp, 
                                            os.path.join(os.sep, os.path.basename(filename_results)), 
                                            i, uuid.uuid4(), frmt)
            fn_res.append(fr)
        job = pool.apply_async(sub.demultiplex_lines_wrapper, (filename_reads, fm, fa, fu, fr, i, subprocesses,))
        jobs.append(job)

    # Collect results from the workers through the pool result queue.
    cdef list statss = list()
    while jobs:
        job = jobs.pop()
        if job.ready():
            statss.append(job.get())
        else:
            jobs.append(job)
            time.sleep(0.01)

    # Wait until finished.
    pool.close()
    pool.join()

    # Merge files.
    if filename_matched != None:
        merge_files(filename_matched, fn_matched)
    if filename_ambig != None:
        merge_files(filename_ambig, fn_ambig)
    if filename_unmatched != None:
        merge_files(filename_unmatched, fn_unmatched)
    if filename_results != None:
        merge_files(filename_results, fn_res)

    # Merge stats
    return merge_stats(statss)

cdef merge_files(str filename, list part_names):
    """
    Merges and deletes temporary files.
    """
    cdef str frmt = filename.split(".")[-1].lower()
    cdef object f
    cdef str part
    cdef object p
    cdef str ln
    cdef object rec
    if frmt == "fa" or frmt == "fasta" or frmt == "fq" or frmt == "fastq" or frmt == "tsv":
        with open(filename, 'w') as f:
            for ln in fileinput.input(part_names):
                f.write(ln)
    elif frmt == "sam" or frmt == "bam":
        write_attrib = "wh" if frmt == "sam" else "wb"
        read_attrib = "r" if frmt == "sam" else "rb" 
        with ps.AlignmentFile(filename, write_attrib, 
                              template=ps.AlignmentFile(part_names[0], 
                                                        read_attrib,
                                                        check_header=True, 
                                                        check_sq=False)) as f:
            for part in part_names:
                with ps.AlignmentFile(part, read_attrib, check_header=True, check_sq=False) as p:
                    for rec in p:
                        f.write(rec)
    else:
        raise ValueError("Unknown file format to merge")

    # Delete parts
    for part in part_names:
        if os.path.isfile(part): os.remove(part)

cdef merge_stats(list statss):
    """
    Merges statistics into one.
    """
    cdef object merged_stats = statss[0]
    for stats in statss[1:]:
        merged_stats += stats
    return merged_stats