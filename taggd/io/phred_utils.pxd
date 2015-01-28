cpdef double probability_to_phred(double p)

cpdef double phred_to_probability(double q)

cpdef str sanger_phred_to_ascii(list phred)

cpdef list sanger_ascii_to_phred(str ascii_seq)

cpdef list sanger_ascii_to_probability(str ascii_seq)

cpdef list illumina_1_3_ascii_to_phred(str ascii_seq)

cpdef list illumina_1_3_ascii_to_probability(str ascii_seq)

cpdef list illumina_1_8_casava_ascii_to_phred(str ascii_seq)

cpdef list illumina_1_8_casava_ascii_to_probability(str ascii_seq)

cpdef list solid_ascii_to_phred(str ascii_seq)

cpdef list solid_ascii_to_probability(str ascii_seq)