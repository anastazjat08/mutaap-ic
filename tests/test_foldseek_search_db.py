import pytest
from unittest.mock import patch, MagicMock
import pandas as pd
import mutaapic.analysis.compare_structures as cs


@patch("shutil.which", return_value="/usr/bin/foldseek")
@patch("subprocess.run")
def test_foldseek_search_db_success(mock_run, mock_which, tmp_path):
    """Test Foldseek search and correct parsing of results."""

    # Filename for the mock result TSV
    result_tsv = tmp_path / "result.tsv"

    result_tsv.write_text(
        "query1\ttargetA\t0.9\t1.0\t0.8\t1e-5\n"
        "query1\ttargetB\t0.7\t2.0\t0.7\t1e-3\n"
    )

    # Patch TemporaryDirectory returns tmp_path
    with patch("tempfile.TemporaryDirectory") as mock_tmp:
        mock_tmp.return_value.__enter__.return_value = tmp_path

        df, ids = cs.foldseek_search_db("query.pdb", "db", k=1)

        assert len(df) == 1
        assert df.iloc[0]["target"] == "targetA"
        assert ids == ["targetA"]


@patch("shutil.which", return_value=None)
def test_foldseek_missing(mock_which):
    """Should raise error when Foldseek is missing."""
    with pytest.raises(EnvironmentError):
        cs.foldseek_search_db("q.pdb", "db")