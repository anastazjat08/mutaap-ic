from Bio.Align import PairwiseAligner

'''
This module is for aligning DNA sequences automatically and finding change ranges.
'''

def automatic_alignment(original_seq: str, modified_seq: str) -> list[str]:
    """
    if chosen --auto-alignment
    alignemnt of DNA sequences
    returns changes in the string tabular format

    needs normalization
    """
    aligner = PairwiseAligner()

    aligner.mode = "global"

    alignment = aligner.align(original_seq, modified_seq)[0]

    aligned_original = alignment[0]
    aligned_modified = alignment[1]


    return extract_nt_mismatches(aligned_original, aligned_modified)

def extract_nt_mismatches(aligned_original: str, aligned_modified: str) -> list[str]:
    report = []
    new_pos = 0

    for o, m in zip(aligned_original, aligned_modified):

        if m != "-":
            new_pos += 1

        if o == m:
            continue

        if o == "-" and m != "-":
            change_type = "insertion"
        elif o != "-" and m == "-":
            change_type = "deletion"
        else:
            change_type = "mismatch"

        report.append(f"{new_pos}\t{o}\t{m}\t{change_type}")

    return report

def merge_adjacent_changes(nt_report: list[str]) -> list[str]:
    """
    merges neighbour changes that map to the same ORF position progression
    
    """

    merged = []
    block = []

    def flush_block() -> None:
        if not block:
            return

        if len(block) == 1:
            merged.append(block[0])
            return

        # extract fields
        positions = [int(x.split("\t")[0]) for x in block]
        originals = [x.split("\t")[1] for x in block]
        modifieds = [x.split("\t")[2] for x in block]
        change_type = block[0].split("\t")[3]

        start = positions[0]
        end = positions[-1]

        merged.append(
            f"{start}-{end}\t{''.join(originals)}\t{''.join(modifieds)}\t{change_type}"
        )

    prev_pos = None
    prev_type = None

    for line in nt_report:
        pos, orig, mod, ctype = line.split("\t")

        pos = int(pos)

        if prev_pos is None:
            block = [line]
        else:
            # if this change continues the block
            if pos == prev_pos or pos == prev_pos + 1 and ctype == prev_type:
                block.append(line)
            else:
                flush_block()
                block = [line]

        prev_pos = pos
        prev_type = ctype

    flush_block()
    return merged
