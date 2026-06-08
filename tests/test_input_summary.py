from pathlib import Path

from mutaapic.analysis.input_summary import generate_input_summary


def test_generate_input_summary(tmp_path):

    output_dir = tmp_path / "report"

    path = generate_input_summary(
        output_dir=str(output_dir),
        mutations=[],
        nt_changes=[],
        original_protein="MAAA",
        modified_protein="MVAA",
        frameshift=False,
    )

    assert Path(path).exists()

    content = Path(path).read_text()

    assert "SUMMARY" in content
    assert "AMINO_ACID_CHANGES" in content
    assert "NUCLEOTIDE_CHANGES" in content