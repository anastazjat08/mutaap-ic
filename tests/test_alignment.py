from mutaapic.orf.alignment import automatic_alignment, extract_nt_mismatches, merge_adjacent_changes

def test_automatic_alignment_detects_mismatch():
    result = automatic_alignment(
        "ATGC",
        "ATTC"
    )

    assert len(result) == 1
    assert result[0].endswith("mismatch")

def test_automatic_alignment_substitution():

    original = "ATGAAA"
    modified = "GTGAAA"

    result = automatic_alignment(original, modified)

    assert result == [
        "1\tA\tG\tmismatch"
    ]

def test_automatic_alignment_deletion():

    original = "ATGCCCAAA"
    modified = "ATGAAA"

    result = automatic_alignment(original, modified)

    assert any("deletion" in row for row in result)

def test_extract_nt_mismatches():

    result = extract_nt_mismatches(
        "ATGC",
        "ATTC"
    )

    assert result == [
        "3\tG\tT\tmismatch"
    ]

def test_merge_adjacent_changes():

    data = [
        "10\tA\tG\tmismatch",
        "11\tT\tC\tmismatch",
        "12\tG\tA\tmismatch"
    ]

    result = merge_adjacent_changes(data)

    assert result == [
        "10-12\tATG\tGCA\tmismatch"
    ]