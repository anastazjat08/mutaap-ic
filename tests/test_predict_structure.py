import os
from unittest.mock import patch, MagicMock
import mutaapic.structure.predict_structure as ps

def test_predict_structure_success(tmp_path):
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.text = "PDBDATA"

    with patch("requests.post", return_value=mock_response):
        out = ps.predictESM("orig", "AAAA", tmp_path)

        assert os.path.exists(out)

        with open(out) as f:
            assert f.read() == "PDBDATA"


def test_predict_structure_failure(tmp_path):
    """Should return None when API fails."""
    mock_response = MagicMock()
    mock_response.status_code = 504
    mock_response.text = "timeout"

    with patch("requests.post", return_value=mock_response):
        out = ps.predictESM("orig", "AAAA", tmp_path)
        assert out is None