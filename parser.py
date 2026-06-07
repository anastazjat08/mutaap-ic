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

# 1.1
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
###

def automatic_alignment(original_seq, modified_seq):
    """
    if chosen --auto-alignment
    alignemnt of DNA sequences
    returns changes in the same format as read_txt

    needs normalization
    """
    aligner = PairwiseAligner()

    aligner.mode = "global"

    alignment = aligner.align(original_seq, modified_seq)[0]

    aligned_original = alignment[0]
    aligned_modified = alignment[1]

    changes = []
    seq_pos = 0

    for o, m in zip(aligned_original, aligned_modified):

        if o != "-":
            seq_pos += 1

        # substitution
        if o != "-" and m != "-" and o != m:
            changes.append(str(seq_pos))

        # deletion
        elif m == "-" and o != "-":
            changes.append(str(seq_pos))

        # insertion
        if o == "-" and m != "-":
            changes.append(str(seq_pos))

    return changes

def normalize_changes(changes):
    """
    merges overlapiing sequences from automatic_alignment into ranges
    """
    points = []

    for c in changes:
        if "-" in c:
            start, end = map(int, c.split("-"))
            points.extend(range(start, end + 1))
        else:
            points.append(int(c))

    points = sorted(set(points))

    if not points:
        return []
    
    merged = []
    start = prev = points[0]
    
    for p in points[1:]:
        if p == prev + 1:
            prev = p
        else:
            merged.append(str(start) if start == prev else f"{start}-{prev}")
            start = prev = p
        
    merged.append(str(start) if start == prev else f"{start}-{prev}")

    return merged

# 1.2
def validate_no_internal_stop(protein_sequence):
    """
    raises an error if codon STOP apperas inside the protein
    """
    if "*" in protein_sequence[:-1]:
        raise ValueError("Internal stop codon detected - protein is non-functional")

def compare_proteins(original, modified):
    """
    checks aminoacids differences
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignment = aligner.align(original, modified)[0]

    aligned_original = alignment[0]
    aligned_modified = alignment[1]

    mutations = []

    aa_position = 0

    for aa1, aa2 in zip(aligned_original, aligned_modified):

        if aa1 != "-":
            aa_position += 1

        if aa1 == aa2:
            continue

        if aa2 == "-":
            mutation_type = "deletion"

        else:
            mutation_type = "nonsynonymous"

        mutations.append(
            {"position": aa_position,
             "original": aa1,
             "modified": aa2,
             "type": mutation_type}
        )
    return mutations


def modified_analysis(
        sequence, original_seq, modified_seq, translation_table=11):
    """
    detects:
    size changes
    frameshifts' consequences
    """
    seq = Seq(sequence)

    diff = (len(modified_seq) - len(original_seq))

    frameshift = (diff % 3 != 0)

    return str(seq.translate(table=translation_table)), frameshift


def generate_input_summary(
        output_dir,
        mutations,
        original_protein,
        modified_protein,
        frameshift,
        nt_changes
        ):
    """
    if directory doesn't exist, creates one: report 
    name of file with report: changes_report.tsv
    """
    os.makedirs(output_dir, exist_ok=True)

    report_path = os.path.join(
        output_dir, "changes_report.tsv"
    )

    with open(report_path, "w") as report:
        
        report.write("SUMMARY\n")

        report.write(f"original_protein_length:\t"
                     f"{len(original_protein)}\n")
        
        report.write(f"modified_protein_length:\t"
                     f"{len(modified_protein)}\n")

        report.write(f"frameshift:\t"
                     f"{frameshift}\n\n")
        

        # aminoacid mutations
        report.write("AMINO_ACID_CHANGES\n")

        report.write("aa_position\t"
                     "original_aa\t"
                     "modified_aa\t"
                     "mutation_type\n"
                     )

        for mutation in mutations:
            report.write(
                f"{mutation['position']}\t"
                f"{mutation['original']}\t"
                f"{mutation['modified']}\t"
                f"{mutation['type']}\n\n"
            )

        # nucleotide changes
        report.write("\nNUCLEOTIDE_CHANGES\n")
        
        for change in nt_changes:
            report.write(
                f"{change}\n")

    return report_path
            
# to main in cli
def main():
    args = parse_args()

    # original_seq = read_fasta(args.original_file)
    # modified_seq = read_fasta(args.modified_file)

    # original_protein = isit_orf(original_seq, args.table)

    # modified_protein, frameshift = (
    #     modified_analysis(
    #         modified_seq, original_seq, modified_seq, args.table))
    
    # validate_no_internal_stop(modified_protein)


    # if args.auto_alignment:
    #     nt_changes = normalize_changes(
    #         automatic_alignment(
    #             original_seq, modified_seq))
        

    # elif args.changes_file:
    #     nt_changes = read_txt(args.changes_file)
    
    # else:
    #     raise ValueError("Provide -a for automatic alignment or your own file")


    # mutations = compare_proteins(original_protein, modified_protein)

    # generate_input_summary(
    #     output_dir=args.output,
    #     mutations=mutations,
    #     nt_changes=nt_changes,
    #     original_protein=original_protein,
    #     modified_protein=modified_protein,
    #     frameshift=frameshift
    # )

    print("First part of analysis is finished")


if __name__ == "__main__":
    main()