import mutaapic.validation.validate_sequence as vs
import pytest

def test_valid_sequence():
    assert vs.validate_aa_sequence("ACDEFGHIKLMNPQRSTVWY")

def test_invalid_sequence():
    assert vs.validate_aa_sequence("AXYZ123") is False