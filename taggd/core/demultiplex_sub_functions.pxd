

cdef object demultiplex_lines(str filename_reads,
                            str filename_matched,
                            str filename_ambig,
                            str filename_unmatched,
                            str filename_res,
                            int ln_offset,
                            int ln_mod)

cdef list demultiplex_record(object rec)

cdef str trim_helpers(str seq)