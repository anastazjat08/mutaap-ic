import pytest

from mutaapic.utils.read_files import read_fasta, read_txt

def test_read_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "atgaaataa\n"
    )

    seq = read_fasta(fasta)

    assert seq == "ATGAAATAA"
    
def test_read_txt(tmp_path):
    
    txt = tmp_path / "changes.txt"
    
    txt.write_text(
        "1\n"
        "4-6\n"
    )

    changes = read_txt(txt)

    assert changes == ["1", "4-6"]