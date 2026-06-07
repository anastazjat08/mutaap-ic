import requests
import os
import re
import time
import shutil
import subprocess
import tempfile
import pandas as pd
import argparse as ap

import mutaapic.structure.predict_structure as predict_structure
import mutaapic.analysis.compare_structures as compare_structures
import mutaapic.validation.validate_sequence as validate_sequence
import mutaapic.utils.filesystem as filesystem
import mutaapic.utils.foldseek as foldseek
import mutaapic.utils.fetch as fetch
import mutaapic.reporting.report as report
from mutaapic.function.predict_function import predict_function, compare_function



def main():
    parser = ap.ArgumentParser(description="MutAAP-IC: A tool for predicting and analyzing the effects of mutations on protein structure.")
    # parser.add_argument("fasta_orig", type=str, help="The FASTA file containing the original nucleotide sequence to analyze.")
    # parser.add_argument("fasta_mut", type=str, help="The FASTA file containing the mutated nucleotide sequence to analyze.")
    parser.add_argument("--db_path", type=str, default="./mutaap_db", help="Path where the default Foldseek databases are stored and downloaded if missing.")
    parser.add_argument("--custom_db", type=str, help="Path to a custom database of structures for comparison (optional).")
    parser.add_argument("--exclude_pdb", action="store_true", help="Whether to exclude comparison with structures from PDB DB supported by Foldseek (default: False).")
    parser.add_argument("--exclude_af", action="store_true", help="Whether to exclude comparison with structures from Alphafold DB supported by Foldseek (default: False).")
    parser.add_argument("--top_k", type=int, default=10, help="Number of top structures found by Foldseek in specific DB to take into statistical analysis (default: 10).")
    parser.add_argument("--out_dir", type=str, default="./mutaap_results", help="Directory to save results (default: mutaap_results).")
    parser.add_argument("--email", default="example@example.com", help="Email for InterProScan API (required by EBI)")

    args = parser.parse_args()

    # initialize result placeholders for reporting
    custom_results_df = None
    custom_db = None
    custom_db_structures = None
    custom_structure_ids = None

    pdb_results_df = None
    pdb_structure_ids = None
    pdb_db_structures = None

    af_results_df = None
    af_structure_ids = None
    af_db_structures = None

    def _resolve_db_root(db_root: str) -> str:
        """Return the directory that actually contains the supported Foldseek databases."""
        candidates = [db_root, os.path.join(db_root, "db")]
        for candidate in candidates:
            pdb_db_path = os.path.join(candidate, "pdb_db")
            af_db_path = os.path.join(candidate, "af_db")

            pdb_exists = foldseek.foldseek_db_exists(pdb_db_path)
            af_exists = foldseek.foldseek_db_exists(af_db_path)

            if pdb_exists or af_exists:
                return candidate

        return db_root

    db_root = _resolve_db_root(args.db_path)

    # ================= SEQUENCE EXTRACTION AND VALIDATION =================
    # HERE PAULA'S PART WITH SEQUENCE EXTRACTION FROM FASTA FILES AND ORF ANALYSIS, PARSERS
    # I wrote simple validation functions in validate_sequence.py - use here, the remove from predict_structure.py


    # I assume we have sequences extracted as strings at this point, e.g.:
    orig_sequence = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVHAAKFESNFNTYKTNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNYTRSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
    mut_sequence = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"


    # ================= STRUCTURE PREDICTION =================
    orig_pdb = predict_structure.predictESM('orig', orig_sequence, args.out_dir)
    time.sleep(5)
    mut_pdb = predict_structure.predictESM('mut', mut_sequence, args.out_dir)

    # ================= STRUCTURE COMPARISON =================
    comparison_results = compare_structures.compare_structures(orig_pdb, mut_pdb)

    # ================= CUSTOM DATABASE COMPARISON =================
    custom_db_source = args.custom_db
    if custom_db_source:
        # Check if Foldseek DB already exists for the given folder, if not create it
        custom_db_foldername = filesystem.get_db_name_from_path(custom_db_source)
        custom_db_path = f"{db_root}/custom_db/{custom_db_foldername}"

        if not foldseek.foldseek_db_exists(custom_db_path):
            print(f"[INFO] Creating Foldseek DB from custom structures in '{custom_db_source}'...")
            custom_db = foldseek.create_own_db(custom_db_source, f"{db_root}/custom_db/{custom_db_foldername}/{custom_db_foldername}") # the last after "/" is the name used for filenames of db
        else:
            print(f"[INFO] Foldseek DB already exists for '{custom_db_source}'. Using existing DB.")
            custom_db = f"{db_root}/custom_db/{custom_db_foldername}/{custom_db_foldername}"
        
        # Perform comparison against the custom database
        custom_results_df, custom_structure_ids = compare_structures.foldseek_search_db(mut_pdb, custom_db, k=args.top_k)



    # ================= SUPPORTED DATABASE COMPARISON =================
    if not args.exclude_pdb:
        # Check if Foldseek DB already exists for the given folder, if not create it
        pdb_db_path = f"{db_root}/pdb_db"

        if not foldseek.foldseek_db_exists(pdb_db_path):
            print(f"[INFO] Creating Foldseek DB from PDB structures...")
            pdb_db = foldseek.create_supported_foldseek_db("PDB", db_root, "pdb_db")
        else:
            print(f"[INFO] Foldseek DB already exists for PDB. Using existing DB.")
            pdb_db = f"{db_root}/pdb_db"

        # Perform comparison against the PDB database
        pdb_results_df, pdb_structure_ids = compare_structures.foldseek_search_db(mut_pdb, pdb_db, k=args.top_k)
        pdb_db_structures = fetch.download_pdb_structures(pdb_structure_ids, f"{args.out_dir}/pdb_db_structures")

    if not args.exclude_af:
        # Check if Foldseek DB already exists for the given folder, if not create it
        af_db_path = f"{db_root}/af_db"

        if not foldseek.foldseek_db_exists(af_db_path):
            print(f"[INFO] Creating Foldseek DB from Alphafold/Swiss-Prot structures...")
            af_db = foldseek.create_supported_foldseek_db("Alphafold/Swiss-Prot", db_root, "af_db")
        else:
            print(f"[INFO] Foldseek DB already exists for Alphafold/Swiss-Prot. Using existing DB.")
            af_db = f"{db_root}/af_db"
        
        
        af_results_df, af_structure_ids = compare_structures.foldseek_search_db(mut_pdb, af_db, k=args.top_k)
        af_db_structures = fetch.download_af_structures(af_structure_ids, f"{args.out_dir}/af_db_structures")


    # ================= FUNCTIONAL PREDICTION =================
    function_results = None
    try:
        orig_go = predict_function(orig_sequence, email=args.email, label="orig")
        mut_go  = predict_function(mut_sequence,  email=args.email, label="mut")
        print("[DEBUG] orig_go:", orig_go)
        print("[DEBUG] mut_go:", mut_go)

        if orig_go and mut_go:
            function_results = compare_function(orig_go, mut_go)

            # Save per-ontology CSVs for downstream use
            for ont, df in function_results.items():
                df.to_csv(f"{args.out_dir}/go_{ont}_comparison.csv", index=False)
        else:
            print("[Function] Skipping comparison - one or both sequences failed annotation.")

    except Exception as e:
        import traceback
        traceback.print_exc() 
        print(f"[Function] Skipping functional annotation: {e}")

    # ================= STATISTICS, REPORT GENERATION =================

    # Generate a report HTML summarizing results
    try:
        report_path = report.generate_report(
            args.out_dir,
            orig_pdb,
            mut_pdb,
            comparison_results,
            function_results=function_results,
            custom_results_df=custom_results_df,
            custom_structure_ids=custom_structure_ids,
            pdb_results_df=pdb_results_df,
            af_results_df=af_results_df,
            custom_downloads=custom_db_structures,
            pdb_downloads=pdb_db_structures,
            af_downloads=af_db_structures,
        )
        print(f"[INFO] Report written to: {report_path}")
    except Exception as e:
        print(f"[WARN] Failed to generate report: {e}")



if __name__ == "__main__":
    main()