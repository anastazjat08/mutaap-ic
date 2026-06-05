import os

'''Utility functions for filesystem operations such as path handling,
folder inspection, and file extension checks.'''


def get_db_name_from_path(folder_path: str) -> str:
    """
    Extracts the last folder name from a given path.
    """
    folder_path = os.path.normpath(folder_path)
    return os.path.basename(folder_path)


