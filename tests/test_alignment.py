from parser import automatic_alignment, normalize_changes

def test_automatic_alignment_substitution():

    original = "ATGAAA"
    modified = "GTGAAA"

    changes = automatic_alignment(original, modified)

    assert normalize_changes(changes) == ["1"]

def test_automatic_alignment_deletion():
    original = "ATGCCCAAA"
    modified = "ATGAAA"

    changes = automatic_alignment(original, modified)

    result = normalize_changes(changes)

    assert result == ["4-6"]