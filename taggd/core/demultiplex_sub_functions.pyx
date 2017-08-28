"""
Contains functions for demultiplexing records.
It uses a multi-thread approach where several file descriptors
are open and every i line of the input file is processed.
The resulting file are then merged
"""

import random
import time
from cpython cimport bool
import taggd.io.reads_reader_writer as rw
import taggd.core.match as match
cimport taggd.core.match as match
import taggd.core.match_type as match_type
cimport taggd.core.match_type as match_type
import taggd.core.statistics as statistics
import taggd.core.demultiplex_search_functions as srch
cimport taggd.core.demultiplex_search_functions as srch

# User settings.
cdef dict true_barcodes
cdef int start_position
cdef int barcode_length
cdef int pre_overhang
cdef int post_overhang
cdef int max_edit_distance
cdef int homopolymer_filter
cdef list homopolymers
cdef str seed
cdef bool multiple_hits_keep_one
cdef list trim_sequences

def init(dict true_barcodes_,
         int start_position_,
         int pre_overhang_,
         int post_overhang_,
         int max_edit_distance_,
         int homopolymer_filter_,
         str seed_,
         bool multiple_hits_keep_one_,
         list trim_sequences_,
         str barcode_tag_):
    """
    Initializes global variables for matching.
    """

    global true_barcodes
    true_barcodes = true_barcodes_

    # User settings.
    global start_position
    start_position = start_position_
    global pre_overhang
    pre_overhang = pre_overhang_
    global post_overhang
    post_overhang = post_overhang_
    global barcode_length
    barcode_length = len(true_barcodes.keys()[0])
    global max_edit_distance
    max_edit_distance = max_edit_distance_
    global homopolymer_filter
    homopolymer_filter = max(homopolymer_filter_, 0)
    global homopolymers
    homopolymers = list()
    if homopolymer_filter > 0:
        for c in "ACGT":
            homopolymers.append(c * homopolymer_filter)
    global seed
    seed = seed_
    global multiple_hits_keep_one
    multiple_hits_keep_one = multiple_hits_keep_one_
    global trim_sequences
    trim_sequences = trim_sequences_
    global barcode_tag
    barcode_tag = barcode_tag_
    # Adjust the barcode length and the overhang if 
    # we want to trim away helpers from the barcode
    if trim_sequences is not None:
        for start,end in trim_sequences:
            barcode_length += (end - start)
   
def demultiplex_lines_wrapper(str filename_reads,
                              str filename_matched,
                              str filename_ambig,
                              str filename_unmatched,
                              str filename_res,
                              int ln_offset,
                              int ln_mod):
    """
    Non cdef wrapper for cdef:ed subprocess function for demultiplexing parts of a file.
    Demultiplexes every ln_mod line, starting at ln_offset, writing to specified files.
    """
    return demultiplex_lines(filename_reads,
                            filename_matched,
                            filename_ambig,
                            filename_unmatched,
                            filename_res,
                            ln_offset,
                            ln_mod)

cdef object demultiplex_lines(str filename_reads,
                              str filename_matched,
                              str filename_ambig,
                              str filename_unmatched,
                              str filename_res,
                              int ln_offset,
                              int ln_mod):
    """
    Demultiplexes every ln_mod line, starting at ln_offset, writing to specified files.
    """

    random.seed(seed)

    # Time
    cdef int start_time = time.time()

    # Open files
    # TODO check they are not open already
    cdef bool header = (ln_offset == 0)
    cdef object re_wr = rw.ReadsReaderWriter(filename_reads)
    re_wr.reader_open()
    cdef object f_match = re_wr.get_writer(filename_matched) if filename_matched else None 
    cdef object f_ambig = re_wr.get_writer(filename_ambig) if filename_ambig else None
    cdef object f_unmatch = re_wr.get_writer(filename_unmatched) if filename_unmatched else None
    cdef object f_res = open(filename_res, "w") if filename_res else None
    if f_res != None and header: f_res.write(match.get_match_header() + "\n")

    # Start demultiplexing.
    cdef object stats = statistics.Statistics(ln_offset, max_edit_distance)
    cdef int i = 0
    cdef int j
    cdef list mtch = None
    cdef object mt = None
    cdef object rec = None
    cdef str bcseq = None
    cdef object bc = None
    cdef list tags = None

    # Iterate over all input lines.
    for rec in re_wr.reader_open():

        # Process every i-th line.
        if i % ln_mod == ln_offset:
            mtch = demultiplex_record(rec)
            stats.total_reads += 1

            # Iterate over all matches (only more than one if ambiguous)
            for mt in mtch:

                # Write to results file.
                if f_res != None: f_res.write("{}\n".format(mt))

                # No match.
                if mt.match_type == match_type.UNMATCHED:
                    if f_unmatch != None:
                        re_wr.write_record(f_unmatch, mt.record)
                        stats.total_reads_wr += 1
                    stats.unmatched += 1
                    continue

                # Append record with properties. B0:Z:Barcode, B1:Z:Prop1, B2:Z:prop3 ...
                bc = true_barcodes[mt.barcode]
                tags = list()
                tags.append(("B0:Z", mt.barcode))
                for j in xrange(len(bc.attributes)):
                    tags.append(("B{}:Z".format(j+1), bc.attributes[j]))
                mt.record.add_tags(tags)

                # Write to output file.
                if mt.match_type == match_type.MATCHED_PERFECTLY:
                    if f_match != None:
                        re_wr.write_record(f_match, mt.record)
                        stats.total_reads_wr += 1
                    stats.perfect_matches += 1
                    stats.edit_distance_counts[0] += 1
                elif mt.match_type == match_type.MATCHED_UNAMBIGUOUSLY:
                    if f_match != None:
                        re_wr.write_record(f_match, mt.record)
                        stats.total_reads_wr += 1
                    stats.imperfect_unambiguous_matches += 1
                    stats.edit_distance_counts[mt.edit_distance] += 1
                elif mt.match_type == match_type.MATCHED_AMBIGUOUSLY:
                    if f_ambig != None:
                        re_wr.write_record(f_ambig, mt.record)
                        stats.total_reads_wr += 1
                    stats.imperfect_ambiguous_matches += 1
                else:
                    raise ValueError("Invalid match type")

        # Next iteration
        i += 1

    # Close files.
    re_wr.reader_close()
    if f_match != None: f_match.close()
    if f_ambig != None: f_ambig.close()
    if f_unmatch != None: f_unmatch.close()
    if f_res != None: f_res.close()

    # Get finish time
    stats.time = time.time() - start_time
    return stats

cdef str trim_helpers(str seq):
    """Simply helper function to remove
    helper sequences from a barcode"""
    cdef int prev_start = 0
    cdef int prev_end = 0
    for start,end in trim_sequences:
        offset = prev_end - prev_start
        seq = seq[:(start-offset)] + seq[(end-offset):]
        prev_start = start
        prev_end = end
    return seq
            
cdef list demultiplex_record(object rec):
    """
    Demultiplexes a record and returns a list of match objects 
    (only more than one if ambiguous).
    """
    # Define local variables to speed up
    cdef str bcseq = None
    cdef int dist = 0
    cdef str read_barcode
    
    # Try perfect hit first.
    if not barcode_tag:
        read_barcode = rec.sequence[start_position:(start_position+barcode_length)]
    else:
        try:
            read_barcode = {tag:value for tag,value in rec.attributes["tags"]}[barcode_tag]
        except KeyError:
            raise ValueError('Error: cannot demultiplex, the specified SAM/BAM tag ("'+barcode_tag+'") is not present for record '+rec.annotation+'.\n')
    if trim_sequences is not None: read_barcode = trim_helpers(read_barcode)

    if read_barcode in true_barcodes:
        return [match.Match(rec, match_type.MATCHED_PERFECTLY, read_barcode, 0)]

    # Homopolymer filter.
    for filter in homopolymers:
        if filter in read_barcode:
            return [match.Match(rec, match_type.UNMATCHED, "-", -1)]
            
    # Include overhang.
    if pre_overhang != 0 or post_overhang != 0:
        read_barcode = rec.sequence[(start_position - pre_overhang):min(len(rec.sequence), \
                                    (start_position + barcode_length + post_overhang))]
        # NOTE should take care of overhang bases more carefully here
        if trim_sequences is not None: read_barcode = trim_helpers(read_barcode)
          
    # Narrow down hits.
    cdef list candidates = srch.get_candidates(read_barcode)
    cdef qual_hits = srch.get_distances(read_barcode, candidates)
    cdef top_hits = srch.get_top_hits(qual_hits)
    
    if not top_hits:
        # UNMATCHED
        return [match.Match(rec, match_type.UNMATCHED, "-", -1)]

    random.shuffle(top_hits)
    if len(top_hits) == 1 or multiple_hits_keep_one:
        # Unambiguous match
        bcseq, dist = top_hits[0]
        return [match.Match(rec, match_type.MATCHED_UNAMBIGUOUSLY, bcseq, dist)]
    
    # Add the rest as ambiguous match
    return [match.Match(rec, match_type.MATCHED_AMBIGUOUSLY, bcseq, dist) \
            for bcseq,dist in top_hits]
