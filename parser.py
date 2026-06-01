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
                        help="FASTA file with ORF in DNA sequence")
    parser.add_argument("modified_file", 
                        help="modified ORF")
    parser.add_argument("-a", "--allignment",
                        help=("do automatic allignment, instead of flag you can put your file:" \
                        "one change in one line, eg.:\n" \
                        "8\n" \
                        "28-32"))
    parser.add_argument("-o", "--output", default="report", 
                        help="output directory")
    parser.add_argument("-t", "--table", type=int, default=11, 
                        help="NCBI translation code table, default is bacterial")

    return parser.parse_args()


def read_fasta(path):
    """ 
    reads first sequence
    returns DNA always in upper string
    """
     
    record = next(SeqIO.parse(path, "fasta"))
    sequence = str(record.seq).upper()

    return sequence


def isit_orf(sequence, translation_table=11):
    """
    checks if given sequence is a valid ORF:
    first three letters are START codon
    last three letters are STOP codon
    STOP codon appears in the middle
    length is divisible by three

    if so, error
    """
    seq = Seq(sequence)


    if len(seq) % 3 != 0:

        raise ValueError("Sequence length is not divisible by 3")
    

    try:
        protein = seq.translate(
            table=translation_table,
            cds=True
        )

    except Exception as e:
        raise ValueError(f"Sequence is not a valid coding sequence: {e}")
    
    return protein

def parse_changes_file(path):
    """
    reads mutation positions from user's file (-a)

    accepts format:
    8
    28-32
    150

    returns a list of positions and/or ranges
    """

    changes = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if "-" in line:
                start, end = map(int, line.split("-"))
                changes.append(
                    {
                        "type": "range",
                        "start": start,
                        "end": end
                    }
                )

            else:
                changes.append(
                    {
                        "type": "single",
                        "position": int(line)
                    }
                )

    return changes


def compare_proteins(original, modified):
    """ 
    checks if aminoacids stay at their
    original position 
    """

    mutations = []

    for pos, (aa1, aa2) in enumerate(
        zip(original, modified),
        start=1
    ):
        if aa1 != aa2:
            mutations.append(
                {
                    "position": pos,
                    "original": aa1,
                    "modified": aa2,
                    "type": "aa_change"
                }
            )

    return mutations


def report_file(
        output_dir,
        mutations,
        original_protein,
        modified_protein,
        functional=True
        ):
    """
    if directory doesn't exist, creates one: report
    
    name is always: changes_report.csv
    
    information to write:
    where AMINOACIDS change
    is protein same length
    functional or has STOP in the middle
    
    """
    os.makedirs(output_dir, exist_ok=True)

    report_path = os.path.join(
        output_dir, "changes_report.tsv"
    )

    same_length = (len(original_protein) == len(modified_protein))

    with open(report_path, "w") as report:

        # header
        report.write("aa_position\t"
                     "original_aa\t"
                     "modified_aa\t"
                     "mutation_type\n"
                     )
        
        # summary
        report.write("\n")
        report.write("SUMMARY\n")

        report.write(f"Original protein length:\t"
                     f"{len(original_protein)}\n")
        report.write(f"modified_protein_length:\t"
                     f"{len(modified_protein)}\n")
        report.write(f"same_length\t"
                     f"{same_length}\n")
        report.write(f"functional\t"
                     f"{functional}\n")

        # mutations
        for mutation in mutations:
            report.write(
                f"{mutation['position']}\t"
                f"{mutation['original']}\t"
                f"{mutation['modified']}\t"
                f"{mutation['type']}\n"
            )

    return report_path
            

def main():
    args = parse_args()
    
    os.makedirs(args.output, exist_ok=True) # if output file doesn't exist it creates it

    original_seq = read_fasta(args.original_file)
    modified_seq = read_fasta(args.modified_file)

    original_protein = isit_orf(original_seq, args.table)
    modified_protein = isit_orf(modified_seq, args.table)

    mutations = compare_proteins(
        original_protein,
        modified_protein
    )

    report_file(
        output_dir=args.output,
        mutations=mutations,
        original_protein=str(original_protein),
        modified_protein=str(modified_protein),
        functional=True
    )

    print("Original protein length:", len(original_protein))
    print("Modified protein length:", len(modified_protein))


if __name__ == "__main__":
    main()