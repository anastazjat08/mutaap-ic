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
        The file path to the first PDB file (reference).
    pdb2 : str
        The file path to the second PDB file (mutated).
    Returns
    -------
    dict
        A dictionary containing the TM-score, RMSD, and raw output from TM-align.
    '''
    # Check if TM-align is available
    if not shutil.which("TMalign"):
        raise EnvironmentError("[ERROR] TM-align is not installed or not in the system PATH.")
    
    def _parse_tmalign_matrix(matrix_path: str):
        if not os.path.exists(matrix_path):
            return None

        rows = []
        with open(matrix_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                numbers = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
                if len(numbers) >= 4:
                    rows.append([float(value) for value in numbers[:4]])
                if len(rows) == 3:
                    break

        return rows if len(rows) == 3 else None

    with tempfile.TemporaryDirectory() as tmpdir:
        matrix_path = os.path.join(tmpdir, "tmalign_matrix.txt")

        # Run TM-align and ask it to write the rotation matrix.
        try:
            result = subprocess.run(
                ["TMalign", pdb1, pdb2, "-m", matrix_path],
                capture_output=True,
                text=True,
                check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"[ERROR] TM-align failed: {e.stderr}")

        output = result.stdout
        matrix = _parse_tmalign_matrix(matrix_path)

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

    # Parse alignment block
    seq1 = None
    similarity = None
    seq2 = None

    lines = [line.strip() for line in output.splitlines()]

    for i, line in enumerate(lines):
        if '":" denotes residue pairs' in line:
            if i + 3 < len(lines):
                seq1 = lines[i + 1]
                similarity = lines[i + 2]
                seq2 = lines[i + 3]
            break

    return {
        "tm_score": tm_score,
        "rmsd": rmsd,
        "alignment": {
            "seq1": seq1,
            "similarity": similarity,
            "seq2": seq2,
        },
        "superposition": matrix,
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