from Bio import SeqIO

from mutaapic.orf.validate_sequence import validate_nt_sequence

'''
Module for reading FASTA and mutation text files
'''

def read_fasta(path:str) -> str:
    """ 
    reads given sequences
    returns them always in upper string
    """
     
    parser = SeqIO.parse(path, "fasta")
    try:
        record = next(parser)
    except StopIteration:
        return None, "[ERROR] FASTA file is empty"

    sequence = str(record.seq).upper()

    if not validate_nt_sequence(sequence):
        return None, "[ERROR] Invalid characters in sequence"

    return sequence


def read_txt(path:str) -> list:
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
    if len(changes) == 0:
        print('[INFO] Text file with changes is empty.')

    return changes