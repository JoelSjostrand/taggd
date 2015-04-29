"""
Main functions for carrying out the demultiplexing with
multithreading.
"""

import multiprocessing as mp
import taggd.core.statistics as statistics
import taggd.core.demultiplex_sub_functions as sub
from cpython cimport bool


def demultiplex(str filename_reads,
                str filename_matched,
                str filename_ambig,
                str filename_unmatched,
                str filename_results,
                int cores):
    """
    Demultiplexes the contents of a reads file.
    """

    if cores == 0:
        cores = mp.cpu_count() - 1

    cdef object pool = mp.Pool(cores)

    # Fire off workers
    cdef object jobs = mp.Queue()
    cdef object job = None
    cdef int i
    cdef str fn_matched = None
    cdef str fn_ambig = None
    cdef str fn_unmatched = None
    cdef str fn_res = None
    cdef str frmt = filename_reads.split(".")[-1]
    for i in xrange(cores):
        if filename_matched != None:
            fn_matched = filename_matched + ".part" + str(i) + "." + frmt
        if filename_ambig != None:
            fn_ambig = filename_ambig + ".part" + str(i) + "." + frmt
        if filename_unmatched != None:
            fn_unmatched = filename_unmatched + ".part" + str(i) + "." + frmt
        if filename_results != None:
            fn_res = filename_results + ".part" + str(i) + "." + frmt
        job = pool.apply_async(sub.demultiplex_lines_wrapper, (fn_matched, fn_ambig, fn_unmatched, fn_res, i, cores,))
        jobs.put(job)

    # Collect results from the workers through the pool result queue.
    cdef list statss = list()
    while not jobs.empty():
        job = jobs.get()
        if job.ready():
            statss.append(job.get())
            print statss[-1]
        else:
            jobs.put(job)

    # Wait until finished.
    pool.close()
    pool.join()

    # Merge files.
    cdef object stats = statistics.Statistics(-1, 0)
    # TODO: Merge stats
    return stats