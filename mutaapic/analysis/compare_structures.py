import subprocess
import shutil
import re
import os
import tempfile
import pandas as pd

'''
This module provides functions for comparing protein structures.
'''

def compare_structures(pdb1: str, pdb2: str) -> dict:
    '''Compares two PDB files using TM-align and returns the TM-score and RMSD.
    Parameters
    ----------
    pdb1 : str
        The file path to the first PDB file.
    pdb2 : str
        The file path to the second PDB file.
    Returns
    -------
    dict
        A dictionary containing the TM-score, RMSD, and raw output from TM-align.
    '''
    # Check if TM-align is available
    if not shutil.which("TMalign"):
        raise EnvironmentError("[ERROR] TM-align is not installed or not in the system PATH.")
    
    # Run TM-align
    try:
        result = subprocess.run(
            ["TMalign", pdb1, pdb2],
            capture_output=True,
            text=True,
            check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"[ERROR] TM-align failed: {e.stderr}")
        
    output = result.stdout

    # Parse TM-score
    tm_score = None
    rmsd = None

    # TM-score for structure 1 -> 2
    m = re.search(r"TM-score=\s*([0-9.]+)", output)
    if m:
        tm_score = float(m.group(1))

    # RMSD
    m = re.search(r"RMSD=\s*([0-9.]+)", output)
    if m:
        rmsd = float(m.group(1))

    return {
        "tm_score": tm_score,
        "rmsd": rmsd,
        "raw_output": output
    }


def foldseek_search_db(query_pdb: str, db_path: str, k: int = 10):
    '''Searches a query PDB file against a Foldseek database and returns the top-k results as a DataFrame and a list of structure IDs.
    Parameters
    ----------
    query_pdb : str
        The file path to the query PDB file.
    db_path : str
        The file path to the Foldseek database to search against.
    k : int
        The number of top results to return (default is 10).
    Returns
    -------
    tuple
        A tuple containing a DataFrame of the top-k search results and a list of structure IDs corresponding to the top-k results.
    '''


    if not shutil.which("foldseek"):
        raise EnvironmentError("Foldseek not found in PATH")

    query_pdb = os.path.abspath(query_pdb)
    db_path = os.path.abspath(db_path)

    with tempfile.TemporaryDirectory() as tmpdir:

        result_tsv = os.path.join(tmpdir, "result.tsv")
        tmp_fs = os.path.join(tmpdir, "tmp")

        cmd = [
            "foldseek",
            "easy-search",
            query_pdb,
            db_path,
            result_tsv,
            tmp_fs,
            "--format-output",
            "query,target,alntmscore,rmsd,lddt,evalue"
        ]

        subprocess.run(cmd, check=True, cwd=tmpdir)

        df = pd.read_csv(
            result_tsv,
            sep="\t",
            names=[
                "query",
                "target",
                "tm",
                "rmsd",
                "lddt",
                "evalue"
            ]
        )

        n_hits = len(df)
        if k > n_hits:
            k = n_hits

        df_top = (
            df.sort_values("tm", ascending=False)
              .head(k)
              .reset_index(drop=True)
        )

        structure_ids = df_top["target"].tolist()

        return df_top, structure_ids





# TESTS
# comp_2_structures = compare_structures("/home/nastka/ADP/mutaap_test/pdb_structures/1JXU.pdb", "/home/nastka/ADP/mutaap_test/pdb_structures/1JXX.pdb")
# print(comp_2_structures)

# create_own_db("/home/nastka/ADP/mutaap_test/pdb_structures", "/home/nastka/ADP/mutaap_test/test_db/pdb_db")
# sup_db = create_supported_foldseek_db("PDB", "/home/nastka/ADP/mutaap_test/test_db", "pdb_db")
# print(sup_db)
# df_results, structure_ids = foldseek_search_db("/home/nastka/ADP/mutaap_test/1CRN.pdb", "/home/nastka/ADP/mutaap_test/test_db/pdb_db")
# print(df_results)
# print(structure_ids)

# sup_db = create_supported_foldseek_db(
#     "Alphafold/Proteome",
#     "/home/nastka/ADP/mutaap_test/alfa_test_db",
#     "alfa_db"
# )
# print(sup_db)
# df_results, structure_ids = foldseek_search_db("/home/nastka/mutaap/mutaap_test/1CRN.pdb", sup_db)
# print(df_results)
# print(structure_ids)

# pdb_list = ['1crn-assembly1_A', '4fc1-assembly1_A', '1ejg-assembly1_A', '1ab1-assembly1_A', '3nir-assembly1_A', '2fd7-assembly1_A', '2fd9-assembly1_A', '2v9b-assembly1_B', '1okh-assembly1_A', '1wuw-assembly1_A']
# af_list = ['AF-Q5Z4W6-F1-model_v6', 'AF-A0A0P0WX99-F1-model_v6', 'AF-Q8VZK8-F1-model_v6', 'AF-Q9C8D6-F1-model_v6', 'AF-Q42597-F1-model_v6', 'AF-Q6ZL66-F1-model_v6', 'AF-O25097-F1-model_v6', 'AF-Q42596-F1-model_v6', 'AF-A0A1C1CIK4-F1-model_v6', 'AF-A0A0J9YTU8-F1-model_v6']

# # pdb_downloads = download_pdb_structures(pdb_list, "/home/nastka/mutaap/mutaap_test/pdb_structures")
# # print(pdb_downloads)
# af_downloads = download_alfaphold_structures(af_list, "/home/nastka/mutaap/mutaap_test/af_structures")
# print(af_downloads)