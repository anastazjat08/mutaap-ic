import subprocess
import sys
from unittest.mock import patch

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
    """CLI should run with minimal arguments."""

    # 1. Create an output directory
    out_dir = tmp_path / "mutaap_results"
    out_dir.mkdir()

    # 2. Create a fake PDB file
    fake_pdb = out_dir / "orig_esmfold_v1.pdb"
    fake_pdb.write_text("MODEL")

    # 3. Mock predictESM to return the path to the fake PDB file
    mock_predict.return_value = str(fake_pdb)

    result = subprocess.run(
        [
            sys.executable,
            "-m", "mutaapic.cli",
            "--exclude_pdb",
            "--exclude_af",
            "--out_dir", str(out_dir)
        ],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0