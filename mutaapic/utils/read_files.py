from Bio import SeqIO


def read_fasta(path):
    """ 
    reads given sequences
    returns them always in upper string
    """
     
    record = next(SeqIO.parse(path, "fasta"))
    sequence = str(record.seq).upper()

    return sequence


def read_txt(path):
    """
    if changes_file provided
    reads given file with written mutations
    returns them in a list
    """
    changes = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line:
                changes.append(line)

    return changes