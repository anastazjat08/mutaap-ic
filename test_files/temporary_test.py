#!/usr/bin/env python3

import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Align import PairwiseAligner

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("original_file", 
                        help="FASTA file with original ORF - DNA sequence")
    parser.add_argument("modified_file", 
                        help="FASTA file with modified ORF - DNA sequence")
    
    parser.add_argument("changes_file", nargs="?", default=None,
                        help="user-provided modifications in given ORF, format should be:" \
                        "one change in one line, eg.:\n" \
                        "8\n" \
                        "28-32")
    parser.add_argument("-a", "--auto-alignment", action="store_true",
                        help="do automatic alignment")

    parser.add_argument("-o", "--output", default="report", 
                        help="output directory")
    parser.add_argument("-t", "--table", type=int, default=11, 
                        help="NCBI translation code table, default is bacterial")

    return parser.parse_args()


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
 
             
def isit_orf(sequence, translation_table=11):
    """
    checks if given sequence is a valid ORF:
    first three letters are START codon
    last three letters are STOP codon
    STOP codon does not appear in the middle
    length is divisible by three

    """
    seq = Seq(sequence)

    if len(seq) % 3 != 0:
        raise ValueError("Sequence length is not divisible by 3")
    
    table = CodonTable.unambiguous_dna_by_id[translation_table]

    start_codons = table.start_codons
    stop_codons = table.stop_codons

    if sequence[:3] not in start_codons:
        raise ValueError("Sequence does not start with valid start codon")
    
    if sequence[-3:] not in stop_codons:
        raise ValueError("Sequence does not end with valid stop codon")

    try:
        protein = seq.translate(table=translation_table, cds=True)
    except Exception as e:
        raise ValueError(f"Sequence is not a valid coding sequence: {e}")
    
    return str(protein)

# def automatic_alignment(original_seq, modified_seq):
    
#     Returns mutation positions in user-compatible format:
#     e.g. ["1", "4-6"]

#     Handles:
#     - substitutions
#     - deletions
#     - insertions
#     and merges contiguous positions correctly.
    """

    aligner = PairwiseAligner()
    alignment = aligner.align(original_seq, modified_seq)[0]

    aligned_orig = alignment[0]
    aligned_mod = alignment[1]

    changes = []
    seq_pos = 0  # position in ORIGINAL sequence

    for o, m in zip(aligned_orig, aligned_mod):

        # advance only on real original base
        if o != "-":
            seq_pos += 1

        # SUBSTITUTION
        if o != "-" and m != "-" and o != m:
            changes.append(seq_pos)

        # DELETION (base exists in original but missing in modified)
        elif o != "-" and m == "-":
            changes.append(seq_pos)

        # INSERTION (gap in original)
        elif o == "-" and m != "-":
            # insertion is attributed to previous position context
            # so we mark current seq_pos (or 1 if at start)
            changes.append(seq_pos if seq_pos > 0 else 1)

    # ---- merge into ranges ----
    if not changes:
        return []

    changes = sorted(set(changes))

    merged = []
    start = prev = changes[0]

    for p in changes[1:]:
        if p == prev + 1:
            prev = p
        else:
            merged.append(str(start) if start == prev else f"{start}-{prev}")
            start = prev = p

    merged.append(str(start) if start == prev else f"{start}-{prev}")

    return merged
"""

    """
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)

    original_seq = read_fasta(args.original_file)
    modified_seq = read_fasta(args.modified_file)

    original_protein = isit_orf(original_seq, args.table)

    modified_protein, length_status, frameshift = (modified_analysis(
                                            modified_seq, original_seq, 
                                            modified_seq, args.table
                                            ))

    mutations = compare_proteins(original_protein, modified_protein)

    internal_stop = "*" in modified_protein[:-1]

    too_short = len(modified_protein) < 20

    functional = (not internal_stop
                    and not frameshift
                    and len(modified_protein) >= 20)


    if args.auto_alignment:
        alignment = automatic_alignment(original_seq, modified_seq)
        nt_changes = nucleotide_changes(original_seq, modified_seq)


    elif args.changes_file:
        nt_changes = read_txt(args.changes_file)
    
    else:
        raise ValueError("Provide -a for automatic alignment or your own file")

    report_file(
        output_dir=args.output,
        mutations=mutations,
        original_protein=str(original_protein),
        modified_protein=str(modified_protein),
        length_status=length_status,
        frameshift=frameshift,
        internal_stop=internal_stop,
        too_short=too_short,
        nt_changes=nt_changes,
        functional=functional
    )

    print("Original protein length:", len(original_protein))
    print("Modified protein length:", len(modified_protein))
    print("Length_status:", length_status)
    print("Frameshift:", frameshift)
    print("Internal stop:", internal_stop) 
    
    """
# not used
# def nucleotide_changes(original_seq, modified_seq):
#     """
#     finds nucleotide subsitutions 
#     """
#     changes = []

#     for pos, (nt1, nt2) in enumerate(
#         zip(original_seq, modified_seq), start=1):
#         if nt1 != nt2:
#             changes.append({
#                 "position": pos,
#                 "original": nt1,
#                 "modified": nt2
#             })

#     return changes
# not used
# def compare_proteins(original, modified):
#     """ 
#     checks if aminoacids stay at their
#     original position 
#     """

#     mutations = []

#     for pos, (aa1, aa2) in enumerate(
#         zip(original, modified), start=1):
        
#         if aa1 != aa2:
#             mutation_type = ("nonsynonymous")

#             if aa2 == "*":
#                 mutation_type = ("stop_gained")
                
#             mutations.append(
#                 {
#                     "position": pos,
#                     "original": aa1,
#                     "modified": aa2,
#                     "type": mutation_type
#                 }
#             )

#     return mutations

def main():
    args = parse_args()

    original_seq = read_fasta(args.original_file)
    modified_seq = read_fasta(args.modified_file)

    original_protein = isit_orf(original_seq, args.table)

    if args.auto_alignment:
        changes = mismatches_to_tsv(original_seq, modified_seq)

        print(changes)


    elif args.changes_file:
        nt_changes = read_txt(args.changes_file)
        print(nt_changes)
    
    else:
        raise ValueError("Provide -a for automatic alignment or your own file")

if __name__ == "__main__":
    main()