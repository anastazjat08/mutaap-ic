#!/usr/bin/env python3

import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
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
    reads first sequence
    returns DNA always in upper string
    """
     
    record = next(SeqIO.parse(path, "fasta"))
    sequence = str(record.seq).upper()

    return sequence


def find_orfs(sequence, translation_table=11, min_length=90):
    """
    predicts ORF's
    """
    orfs = []
    
    for frame in range(3):
        for start in range(frame, len(sequence)-2, 3):
            codon = sequence[start:start + 3]

            try:
                if Seq(codon).translate(table=translation_table) != "M":
                    continue
            except Exception:
                continue

            for stop in range(start + 3, len(sequence) - 2, 3):
                stop_codon = sequence[stop:stop + 3]

                if stop_codon in {"TAA", "TAG", "TGA"}:
                    length = stop + 3 - start

                    if length >= min_length:
                        orfs.append(
                            {
                                "start": start + 1,
                                "end": stop + 3,
                                "length": length
                            }
                        )
                    break
    return orfs
                

def isit_orf(sequence, translation_table=11):
    """
    checks if given sequence is a valid ORF:
    first three letters are START codon
    last three letters are STOP codon
    STOP codon appears in the middle
    length is divisible by three

    """
    seq = Seq(sequence)

    if len(seq) % 3 != 0:
        raise ValueError("Sequence length is not divisible by 3")

    try:
        protein = seq.translate(table=translation_table, cds=True)
    except Exception as e:
        raise ValueError(f"Sequence is not a valid coding sequence: {e}")
    
    return str(protein)


def validate_modified_orf(sequence, translation_table=11):
    pass


def automatic_alignment(original_seq, modified_seq):
    """
    alignemnt of DNA sequences if user did not give their's
    """
    aligner = PairwiseAligner()

    alignment = aligner.align(original_seq, modified_seq)[0]

    return alignment


def nucleotide_changes(original_seq, modified_seq):
    """
    finds nucleotide subsitutions 
    """
    changes = []

    for pos, (nt1, nt2) in enumerate(
        zip(original_seq, modified_seq), start=1):
        if nt1 != nt2:
            changes.append({
                "position": pos,
                "original": nt1,
                "modified": nt2
            })

    return changes


def compare_proteins(original, modified):
    """ 
    checks if aminoacids stay at their
    original position 
    """

    mutations = []

    for pos, (aa1, aa2) in enumerate(
        zip(original, modified), start=1):
        
        if aa1 != aa2:
            mutation_type = ("nonsynonymous")

            if aa2 == "*":
                mutation_type = ("stop_gained")
                
            mutations.append(
                {
                    "position": pos,
                    "original": aa1,
                    "modified": aa2,
                    "type": mutation_type
                }
            )

    return mutations


def modified_analysis(  sequence,
                        original_seq, 
                        modified_seq,
                        translation_table=11):
    """
    detects:
    stops in the middle of sequence
    gained stops
    frameshifts' consequences
    size changes
    """
    seq = Seq(sequence)

    diff = (len(modified_seq) - len(original_seq))

    if diff == 0:
        status = "same"
    elif diff > 0:
        status = "longer"
    else:
        status = "shorter"

    frameshift = ( diff % 3 != 0)

    return str(seq.translate(table=translation_table)), status, frameshift


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


# def indels_detection(alignment):
#     """
#     extracts insertions and deletions from PairwiseAligner output
#     """
#     aligned_orig, aligned_mod = alignment[0]

#     indels = []

#     for i, (o, m) in enumerate(zip(aligned_orig, aligned_mod), start=1):
        
#         if o == "-" and m != "-":
#             indels.append({
#                 "position": i,
#                 "type": "insertion",
#                 "base": m
#             })

#         elif m == "-" and o != "-":
#             indels.append({
#                 "position": i,
#                 "type": "deletion",
#                 "base": o
#             })

#     return indels


# def apply_frameshift(sequence, indels):
    #     """
    #     applies indels to sequence and returns shifted sequence
    #     """

    #     seq = list(sequence)

    #     shift = 0

    #     for indel in sorted(indels, key=lambda x: x["position"]):
    #         pos = indel["position"] + shift

    #         if indel["type"] == "deletion":
    #             if pos < len(seq):
    #                 seq.pop(pos)
    #                 shift -= 1

    #         elif indel["type"] == "insertion":
    #             seq.insert(pos, indel["base"])
    #             shift += 1

    #     return "".join(seq)



# def translate_with_frameshift(sequence, table=11):
#     """
#     translates sequence after indel modifications
#     """
#     seq = Seq(sequence)
#     return str(seq.translate(table=table))


def report_file(
        output_dir,
        mutations,
        original_protein,
        modified_protein,
        length_status,
        frameshift,
        internal_stop,
        too_short,
        nt_changes,
        indels,
        functional=True
        ):
    """
    if directory doesn't exist, creates one: report
    
    name is always: changes_report.csv
    
    information to write:
    where AMINOACIDS change
    is protein same length
    functional or has STOP in the middle
    etc.
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
        report.write(f"length_status\t"
                     f"{length_status}\n")
        report.write(f"frameshift\t"
                     f"{frameshift}\n")
        report.write(f"internal_stop\t"
                     f"{internal_stop}\n")
        report.write(f"functional\t"
                     f"{functional}\n")
        report.write(f"protein_too_short\t"
                     f"{too_short}\n")        

        # frameshift simulation
        # report.write("\nFRAMESHIFT_SIMULATION\n")
        # report.write(f"shifted_protein_length\t{len(shifted_protein)}\n")
        # report.write(f"frameshift_effect\t{frameshift_protein_diff}\n")

        # mutations
        for mutation in mutations:
            report.write(
                f"{mutation['position']}\t"
                f"{mutation['original']}\t"
                f"{mutation['modified']}\t"
                f"{mutation['type']}\n"
            )

        # nucleotide changes
        report.write("\nNUCLEOTIDE_CHANGES\n")
        report.write("position\toriginal_nt\tmodified_nt\n")
        
        for change in nt_changes:
            report.write(
                f"{change['position']}\t"
                f"{change['original']}\t"
                f"{change['modified']}\n")
            
        # indels
        report.write("\nINDELS\n")
        report.write("position\ttype\tbase\n")
        
        for indel in indels:
            report.write(
                f"{indel['position']}\t"
                f"{indel['type']}\t"
                f"{indel['base']}\n")

    return report_path
            

def main():
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

        # indels = indels_detection(alignment)
        # shifted_seq = apply_frameshift(original_seq, indels)
        # shifted_protein = translate_with_frameshift(shifted_seq, args.table)

        indels = []

    elif args.changes_file:
        nt_changes = parse_changes_file(args.changes_file)
        indels = []
    
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
        indels=indels, 
        functional=functional
    )

    print("Original protein length:", len(original_protein))
    print("Modified protein length:", len(modified_protein))
    print("Length_status:", length_status)
    print("Frameshift:", frameshift)
    print("Internal stop:", internal_stop) 


if __name__ == "__main__":
    main()