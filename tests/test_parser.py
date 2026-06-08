import pytest

from mutaapic.orf.validate_sequence import isit_orf, validate_no_internal_stop
from mutaapic.analysis.aa_sequence_analysis import compare_proteins, modified_analysis, merge_adjacent_aa_changes

def test_isit_orf_valid():
    """
    valid CDS:
    ATG AAA TAA
    M K
    """
    seq = "ATGAAATAA"
    protein = isit_orf(seq, 1)
    assert protein == "MK"

def test_isit_orf_no_start():
    seq = "AAAAAATAA"

    with pytest.raises(ValueError):
        isit_orf(seq, 1)

def test_isit_orf_no_stop():
    seq = "ATGAAAAAA"

    with pytest.raises(ValueError):
        isit_orf(seq, 1)

def test_isit_orf_internal_stop():
    seq = "ATGAAATAATAA"
    with pytest.raises(ValueError):
        isit_orf(seq, 1)

def test_validate_no_internal_stop_ok():
    validate_no_internal_stop("MKLQ")

def test_validate_no_internal_stop_error():
    with pytest.raises(ValueError):
        validate_no_internal_stop("MK*LQ")

def test_compare_proteins():

    result = compare_proteins(
        "MAAAA",
        "MAVAA"
    )

    assert result == [
        {
            "position": 3,
            "original": "A",
            "modified": "V",
            "type": "nonsynonymous"
        }
    ]

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

def test_merge_adjacent_aa_changes():

    mutations = [
        {
            "position": 5,
            "original": "A",
            "modified": "V",
            "type": "nonsynonymous"
        },
        {
            "position": 6,
            "original": "T",
            "modified": "G",
            "type": "nonsynonymous"
        }
    ]

    result = merge_adjacent_aa_changes(mutations)

    assert result == [
        {
            "position": "5-6",
            "original": "AT",
            "modified": "VG",
            "type": "nonsynonymous"
        }
    ]