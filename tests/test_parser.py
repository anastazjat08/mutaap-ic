import pytest

from parser import(
    isit_orf,
    normalize_changes,
    validate_no_internal_stop,
    compare_proteins,
    modified_analysis
)

def test_isit_orf_valid():
    """
    valid CDS:
    ATG AAA TAA
    M K
    """
    seq = "ATGAAATAA"
    protein = isit_orf(seq)
    assert protein == "MK"

def test_isit_orf_no_start():
    seq = "AAAAAATAA"

    with pytest.raises(ValueError):
        isit_orf(seq)

def test_isit_orf_no_stop():
    seq = "ATGAAAAAA"

    with pytest.raises(ValueError):
        isit_orf(seq)

def test_isit_orf_internal_stop():
    seq = "ATGAAATAATAA"
    with pytest.raises(ValueError):
        isit_orf(seq)

def test_normalize_changes():
    changes = ["1", "4", "5", "6"]
    assert normalize_changes(changes) == ["1", "4-6"]

def test_validate_no_internal_stop_ok():
    validate_no_internal_stop("MKLQ")

def test_validate_no_internal_stop_error():
    with pytest.raises(ValueError):
        validate_no_internal_stop("MK*LQ")

def test_compare_proteins_one_change():
    mutations = compare_proteins("MKT", "MRT")
    assert len(mutations) == 1
    assert mutations[0]["position"] == 2   
    assert mutations[0]["original"] == "K"
    assert mutations[0]["modified"] == "R"

def test_modified_analysis_no_frameshift():
    protein, frameshift = modified_analysis(
        "ATGAAATAA",
        "ATGAAATAA",
        "ATGAAATAA"
    )

    assert frameshift is False

def test_modified_analysis_frameshift():
    protein, frameshift = modified_analysis(
        "ATGAATAA",
        "ATGAAATAA",
        "ATGAATAA"
    )

    assert frameshift is True