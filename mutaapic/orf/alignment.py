from Bio.Align import PairwiseAligner

'''
This module is for aligning DNA sequences automatically and finding change ranges.
'''

def automatic_alignment(original_seq, modified_seq) -> list:
    """
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

def normalize_changes(changes: list) -> list:
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
