import subprocess
import sys
from unittest.mock import patch
from mutaapic.cli import main

def test_cli_help():
    """CLI should run and show help."""
    result = subprocess.run(
        [sys.executable, "-m", "mutaapic.cli", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "MutAAP-IC" in result.stdout


@patch("mutaapic.structure.predict_structure.predictESM")
@patch("mutaapic.analysis.compare_structures.compare_structures", return_value="OK")

def test_cli_runs_minimal(mock_compare, mock_predict, tmp_path):
    out_dir = tmp_path / "mutaap_results"
    out_dir.mkdir()

    fake_pdb = out_dir / "orig_esmfold_v1.pdb"
    fake_pdb.write_text("MODEL")

    # set mock predictESM
    mock_predict.return_value = str(fake_pdb)

    orig_fa = out_dir / "orig.fasta"
    mod_fa = out_dir / "mod.fasta"
    orig_fa.write_text(">a\nATGAAATAA")
    mod_fa.write_text(">b\nATGTTTTAA")

    sys.argv = [
        "mutaapic.cli",
        str(orig_fa),
        str(mod_fa),
        "--exclude_pdb",
        "--exclude_af",
        "--out_dir", str(out_dir)
    ]

    main()