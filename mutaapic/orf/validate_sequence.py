from Bio.Data import CodonTable
from Bio.Seq import Seq

'''
Module for validating nucleotide sequences and ORF correctness using codon translation rules.
'''

def validate_nt_sequence(sequence :str) -> bool:
    valid_nucleotides = set("ACGTU")
    return all(nucleotide in valid_nucleotides for nucleotide in sequence.upper())


def isit_orf(sequence: str, translation_table: int) -> str:
    """
    checks if given sequence is a valid ORF:
    first three letters are START codon
    last three letters are STOP codon
    STOP codon does not appear in the middle
    length is divisible by three
    """
    seq = Seq(sequence)

    if len(seq) % 3 != 0:
        raise ValueError("[ERROR] Sequence length is not divisible by 3")
    
    table = CodonTable.unambiguous_dna_by_id[translation_table]

    start_codons = table.start_codons
    stop_codons = table.stop_codons

    if sequence[:3] not in start_codons:
        raise ValueError("[ERROR] Sequence does not start with valid start codon")
    
    if sequence[-3:] not in stop_codons:
        raise ValueError("[ERROR] Sequence does not end with valid stop codon")

    try:
        protein = seq.translate(table=translation_table, cds=True)
    except Exception as e:
        raise ValueError(f"[ERROR] Sequence is not a valid coding sequence: {e}")
    
    return str(protein)


def validate_no_internal_stop(protein_sequence: str) -> None:
    """
    raises an error if codon STOP apperas inside the protein
    """
    if "*" in protein_sequence[:-1]:
        raise ValueError("[ERROR] Internal stop codon detected - protein is non-functional")