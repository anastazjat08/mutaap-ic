import pytest
from unittest.mock import patch, MagicMock
import mutaapic.analysis.compare_structures as cs
from subprocess import CalledProcessError


@patch("shutil.which", return_value="/usr/bin/TMalign")
@patch("subprocess.run")
def test_compare_structures_success(mock_run, mock_which, tmp_path):
    """Test successful TM-align run and correct parsing."""
    
    # Fake TM-align output
    fake_output = """
    TM-score= 0.8765
    RMSD= 1.23
    """

    mock_result = MagicMock()
    mock_result.stdout = fake_output
    mock_run.return_value = mock_result

    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text("MODEL")
    pdb2.write_text("MODEL")

    result = cs.compare_structures(str(pdb1), str(pdb2))

    assert result["tm_score"] == 0.8765
    assert result["rmsd"] == 1.23
    assert "TM-score" in result["raw_output"]


@patch("shutil.which", return_value=None)
def test_compare_structures_no_tmalign(mock_which):
    """Should raise error when TMalign is missing."""
    with pytest.raises(EnvironmentError):
        cs.compare_structures("a.pdb", "b.pdb")


@patch("shutil.which", return_value="/usr/bin/TMalign")
@patch("subprocess.run", side_effect=CalledProcessError(1, "TMalign"))
def test_compare_structures_failure(mock_run, mock_which):
    """Should raise RuntimeError when TM-align crashes."""
    with pytest.raises(RuntimeError):
        cs.compare_structures("a.pdb", "b.pdb")