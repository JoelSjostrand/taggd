
def coroutine(func):
    """
    Coroutine decorator, starts coroutines upon initialization.
    """
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start



def readfq(fp): # this is a generator function
    """
    Heng Li's fasta/fastq reader function.
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        #name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


@coroutine
def writefq(fp):  # This is a coroutine
    """
    Fastq writing generator sink.
    Send a (header_comments, sequence, quality) triple to the instance to write it to
    the specified file pointer.
    """
    fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
    try:
        while True:
            record = yield
            read = fq_format.format(header_comments=record[0], sequence=record[1], quality=record[2])
            fp.write(read)

    except GeneratorExit:
        return


@coroutine
def writefa(fp):  # This is a coroutine
    """
    Fasta writing generator sink.
    Send a (header_comments, sequence) double to the instance to write it to
    the specified file pointer.
    """
    fa_format = '>{header_comments}\n{sequence}\n'
    try:
        while True:
            record = yield
            read = fa_format.format(header_comments=record[0], sequence=record[1])
            fp.write(read)

    except GeneratorExit:
        return


def writefq_record(fp, record):
    """
    Send a (header_comments, sequence, quality) triple to the instance to write it to
    the specified file pointer.
    """
    fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
    read = fq_format.format(header_comments=record[0], sequence=record[1], quality=record[2])
    fp.write(read)



def writefa_record(fp, record):
    """
    Send a (header_comments, sequence) double to the instance to write it to
    the specified file pointer.
    """
    fa_format = '>{header_comments}\n{sequence}\n'
    read = fa_format.format(header_comments=record[0], sequence=record[1])
    fp.write(read)