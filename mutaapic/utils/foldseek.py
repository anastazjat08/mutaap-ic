import os
import subprocess
import shutil
import tempfile

'''Utility functions for handling Foldseek databases, including checking for existing databases,
and creating new databases.'''

def foldseek_db_exists(db_folder: str) -> bool:
    """
    Checks whether a Foldseek database exists in the given folder by verifying
    that required file extensions are present among the files.

    Parameters
    ----------
    db_folder : str
        Path to the folder containing Foldseek database files.

    Returns
    -------
    bool
        True if all required Foldseek database file types are present.
    """

    required_ext = [".dbtype", ".index", ".lookup", ".source"]

    if not os.path.isdir(db_folder):
        return False

    files = os.listdir(db_folder)

    # Check if all required file types are present
    for ext in required_ext:
        if not any(f.endswith(ext) for f in files):
            return False

    return True

def create_own_db(input_folder: str, db_path: str) -> str:
    '''
    Creates a Foldseek database from PDB files in the specified input folder.
    Parameters
    ----------
    input_folder : str
        The folder containing PDB files to be included in the database.
    db_path : str
        The path where the Foldseek database will be created.
    '''
    if not shutil.which("foldseek"):
        raise EnvironmentError("[ERROR] Foldseek is not installed or not in the system PATH.")
    
    # Create the database
    if not os.path.isdir(input_folder):
        raise ValueError(f"[ERROR] Folder does not exist: {input_folder}")
    
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    cmd = [
        "foldseek", "createdb",
        input_folder,
        db_path
    ]

    subprocess.run(cmd, check=True)
    print(f"[OK] Foldseek DB '{db_path}' created.")

    return db_path

def create_supported_foldseek_db(db_name: str, out_dir: str, prefix: str) -> str:
    '''Creates a Foldseek database from a supported database name and returns the path to the created database.
    Parameters
    ----------
    db_name : str
        The name of the supported Foldseek database to create (Alphafold/Swiss-Prot, PDB).
    out_dir : str
        The directory where the Foldseek database will be created. If the directory does not exist, it will be created.
    prefix : str
        The prefix for the created database file.
    Returns
    -------
    str
        The file path to the created Foldseek database.
    '''
    if not shutil.which("foldseek"):
        raise EnvironmentError("Foldseek not installed")

    os.makedirs(out_dir, exist_ok=True)

    db_path = os.path.join(out_dir, prefix)
    tmp_path = os.path.join(out_dir, "tmp")

    os.makedirs(tmp_path, exist_ok=True)

    subprocess.run([
        "foldseek",
        "databases",
        db_name,
        db_path,
        tmp_path
    ], check=True)

    return db_path