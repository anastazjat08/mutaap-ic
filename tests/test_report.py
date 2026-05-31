from pathlib import Path

import mutaapic.reporting.report as report


def test_generate_report_includes_structure_viewer(tmp_path):
    out_dir = tmp_path / "results"
    out_dir.mkdir()

    orig_pdb = tmp_path / "orig_esmfold_v1.pdb"
    mut_pdb = tmp_path / "mut_esmfold_v1.pdb"
    orig_pdb.write_text(
        "ATOM      1  CA  ALA A   1      10.000  11.000  12.000  1.00 20.00           C\nEND\n",
        encoding="utf-8",
    )
    mut_pdb.write_text(
        "ATOM      1  CA  ALA A   1      20.000  21.000  22.000  1.00 20.00           C\nEND\n",
        encoding="utf-8",
    )

    report_path = report.generate_report(
        str(out_dir),
        str(orig_pdb),
        str(mut_pdb),
        {
            "tm_score": 0.91234,
            "rmsd": 1.234,
            "alignment": {"seq1": "AAA", "similarity": ":.:", "seq2": "BBB"},
            "superposition": [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
        },
    )

    report_html = Path(report_path).read_text(encoding="utf-8")

    assert "Preliminary interpretation" in report_html
    assert "probably tolerated" in report_html
    # New combined viewer and legend
    assert "Structure View" in report_html
    assert "structure_viewer" in report_html
    assert "3Dmol-min.js" in report_html
    assert "Original (blue)" in report_html
    assert "Mutant (red)" in report_html
    assert "0.9123" in report_html
    assert "1.234" in report_html