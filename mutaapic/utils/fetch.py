import os
import subprocess

'''Utility functions for fetching data from external sources, such as downloading PDB files from the RCSB PDB database.'''

def download_pdb_structures(pdb_ids: list, out_folder: str):
    '''Downloads PDB files for a list of PDB IDs and saves them to a specified folder.
    Parameters
    ----------
    pdb_ids : list
        A list of PDB IDs to download.
    out_folder : str
        The folder path where the downloaded PDB files will be saved.
    Returns
    -------
    dict
        A dictionary mapping each PDB ID to the file path of the downloaded PDB file. If a download fails, the value will be None.
    '''
    os.makedirs(out_folder, exist_ok=True)
    results = {}

    for pdb_id in pdb_ids:
        pdb_code = pdb_id.split("-")[0].upper()
        url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
        out_path = os.path.join(out_folder, f"{pdb_code}.pdb")

        print(f"[INFO] wget {url}")

        try:
            subprocess.run(
                ["wget","-q", "-O", out_path, url],
                check=True
            )
            print(f"[OK] {pdb_code} -> {out_path}")
            results[pdb_id] = out_path

        except subprocess.CalledProcessError:
            print(f"[FAIL] Could not download {pdb_code}")
            results[pdb_id] = None

    return results

def download_alfaphold_structures(af_ids: list, out_folder: str):
    '''Downloads AlphaFold PDB files for a list of AlphaFold IDs and saves them to a specified folder.
    Parameters
    ----------
    af_ids : list
        A list of AlphaFold IDs to download.
    out_folder : str
        The folder path where the downloaded PDB files will be saved.
    Returns
    -------
    dict
        A dictionary mapping each AlphaFold ID to the file path of the downloaded PDB file. If a download fails, the value will be None.
    '''
    os.makedirs(out_folder, exist_ok=True)
    results = {}

    for af_id in af_ids:
        url = f"https://alphafold.ebi.ac.uk/files/{af_id}.pdb"
        out_path = os.path.join(out_folder, f"{af_id}.pdb")

        print(f"[INFO] wget {url}")

        try:
            subprocess.run(
                ["wget", "-O", out_path, url],
                check=True
            )
            print(f"[OK] {af_id} -> {out_path}")
            results[af_id] = out_path

        except subprocess.CalledProcessError:
            print(f"[FAIL] Could not download {af_id}")
            results[af_id] = None

    return results