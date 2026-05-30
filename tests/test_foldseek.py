from unittest.mock import patch
import mutaapic.utils.foldseek as fs

@patch("subprocess.run")
def test_foldseek_db_exists(mock_run, tmp_path):
    """Check DB detection logic."""
    db = tmp_path / "db"
    db.mkdir()
    for ext in ["dbtype", "index", "lookup", "source"]:
        (db / f"db.{ext}").touch()

    assert fs.foldseek_db_exists(str(db))


@patch("subprocess.run")
def test_create_own_db(mock_run, tmp_path):
    """Check that Foldseek createdb is called."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "a.pdb").touch()

    out = tmp_path / "db" / "mydb"
    fs.create_own_db(str(input_dir), str(out))

    mock_run.assert_called_once()