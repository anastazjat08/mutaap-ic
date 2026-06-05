'''
This module provides functions to validate amino acid and nucleotide sequences.
'''

def validate_aa_sequence(sequence :str) -> bool:
    if not sequence or not isinstance(sequence, str):
        raise ValueError("Input sequence is empty or not a string.")
    
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    return all(residue in valid_amino_acids for residue in sequence.upper())

def validate_nt_sequence(sequence :str) -> bool:
    if not sequence or not isinstance(sequence, str):
        raise ValueError("Input sequence is empty or not a string.")
    
    valid_nucleotides = set("ACGTU")
    return all(nucleotide in valid_nucleotides for nucleotide in sequence.upper())

