from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


'''
This module is for analyzing amino acid sequences by aligning proteins and detecting modified sequence consequences.
'''

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


def modified_analysis(sequence, original_seq, modified_seq, translation_table):
    """
    detects:
    size changes
    frameshifts' consequences
    """
    seq = Seq(sequence)

    diff = (len(modified_seq) - len(original_seq))

    frameshift = (diff % 3 != 0)

    return str(seq.translate(table=translation_table)), frameshift
