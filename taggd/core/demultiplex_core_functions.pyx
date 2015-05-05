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


def demultiplex(str filename_reads,\
                str filename_matched,\
                str filename_ambig,\
                str filename_unmatched,\
                str filename_results,\
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
    cdef str frmt = get_format(filename_reads)
    cdef str tmp = tempfile.gettempdir()
    if tmp == "" or tmp == None:
        tmp = "."
    for i in xrange(subprocesses):
        if filename_matched != None:
            fm = tmp + "/" + get_relative_path(filename_matched) + "_part" + str(i) + "_" + str(uuid.uuid4()) + "." + frmt
            fn_matched.append(fm)
        if filename_ambig != None:
            fa = tmp + "/" + get_relative_path(filename_ambig) + "_part" + str(i) + "_" + str(uuid.uuid4()) + "." + frmt
            fn_ambig.append(fa)
        if filename_unmatched != None:
            fu = tmp + "/" + get_relative_path(filename_unmatched) + "_part" + str(i) + "_" + str(uuid.uuid4()) + "." + frmt
            fn_unmatched.append(fu)
        if filename_results != None:
            fr = tmp + "/" + get_relative_path(filename_results) + "_part" + str(i) + "_" + str(uuid.uuid4()) + "." + frmt
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
    """Merges and deletes temporary files."""

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
    elif frmt == "sam":
        with ps.AlignmentFile(filename, "wh", template=ps.AlignmentFile(part_names[0], "r", check_header=True, check_sq=False)) as f:
            for part in part_names:
                with ps.AlignmentFile(part, "r", check_header=True, check_sq=False) as p:
                    for rec in p:
                        f.write(rec)
    elif frmt == "bam":
        with ps.AlignmentFile(filename, "wb", template=ps.AlignmentFile(part_names[0], "rb", check_header=True, check_sq=False)) as f:
            for part in part_names:
                with ps.AlignmentFile(part, "rb", check_header=True, check_sq=False) as p:
                    for rec in p:
                        f.write(rec)
    else:
        raise ValueError("Unknown file format to merge")

    # Delete parts
    for part in part_names:
        os.remove(part)


cdef merge_stats(list statss):
    """Merges statistics into one."""

    cdef object stats = statss[0]
    cdef int i
    for i in xrange(1,len(statss)):
        stats += statss[i]
    return stats


cdef get_format(str filename):
    """Returns the format of a file."""
    return filename.split(".")[-1]


cdef get_relative_path(str filename):
    """Returns the relative filename of a file"""
    if "/" in filename:
        return filename.split("/")[-1]
    if "\\" in filename:
        return filename.split("\\")[-1]
    return filename