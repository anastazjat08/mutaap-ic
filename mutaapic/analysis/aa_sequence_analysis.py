from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from typing import Any

'''
This module is for analyzing amino acid sequences by aligning proteins and detecting modified sequence consequences.
'''

def compare_proteins(original: str, modified:str) -> list[dict[str, Any]]:
    """
    checks aminoacids differences
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignment = aligner.align(original, modified)[0]

    aligned_original = alignment[0]
    aligned_modified = alignment[1]

    mutations = []
    new_pos = 0

    for aa1, aa2 in zip(aligned_original, aligned_modified):

        if aa2 != "-":
            new_pos += 1

        if aa1 == aa2:
            continue

        if aa2 == "-":
            mutation_type = "deletion"

        else:
            mutation_type = "nonsynonymous"

        mutations.append(
            {"position": new_pos,
             "original": aa1,
             "modified": aa2,
             "type": mutation_type}
        )
    return mutations

def modified_analysis(
        original_seq: str, 
        modified_seq: str,
        translation_table: int
        ) -> tuple[str, bool]:
    """
    detects:
    size changes
    frameshifts' consequences
    """

    diff = (len(modified_seq) - len(original_seq))

    frameshift = (diff % 3 != 0)

    return str(Seq(modified_seq).translate(table=translation_table)), frameshift

def merge_adjacent_aa_changes(mutations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    merges neighbour amino acid changes that map to the same or consecutive aa positions

    """

    if not mutations:
        return []

    merged = []
    block = [mutations[0]]

    for m in mutations[1:]:
        prev = block[-1]

        if m["position"] == prev["position"] or m["position"] == prev["position"] + 1:
            block.append(m)
        else:
            merged.append(block)
            block = [m]

    merged.append(block)

    # change to singular records
    final = []
    for block in merged:
        if len(block) == 1:
            final.append(block[0])
        else:
            start = block[0]["position"]
            end = block[-1]["position"]

            original = "".join(m["original"] for m in block)
            modified = "".join(m["modified"] for m in block)
            mutation_type = block[0]["type"]

            final.append({
                "position": f"{start}-{end}",
                "original": original,
                "modified": modified,
                "type": mutation_type
            })

    return final
