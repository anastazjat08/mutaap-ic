import mutaapic.utils.filesystem as fs

def test_get_db_name_from_path():
    assert fs.get_db_name_from_path("/path/to/db") == "db"
    assert fs.get_db_name_from_path("folder/abc") == "abc"