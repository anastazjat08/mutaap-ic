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



def main():
    parser = ap.ArgumentParser(description="MutAAP-IC: A tool for predicting and analyzing the effects of mutations on protein structure.")
    # parser.add_argument("fasta_orig", type=str, help="The FASTA file containing the original nucleotide sequence to analyze.")
    # parser.add_argument("fasta_mut", type=str, help="The FASTA file containing the mutated nucleotide sequence to analyze.")
    parser.add_argument("--own_db", type=str, help="Path to a custom database of structures for comparison (optional).")
    parser.add_argument("--exclude_pdb", action="store_true", help="Whether to exclude comparison with structures from PDB DB supported by Foldseek (default: False).")
    parser.add_argument("--exclude_af", action="store_true", help="Whether to exclude comparison with structures from Alphafold DB supported by Foldseek (default: False).")
    parser.add_argument("--top_k", type=int, default=10, help="Number of top structures found by Foldseek in specific DB to take into statistical analysis (default: 10).")
    parser.add_argument("--out_dir", type=str, default="./mutaap_results", help="Directory to save results (default: mutaap_results).")
    
    args = parser.parse_args()

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
    comparison_results = compare_structures.compare_structures(orig_pdb, mut_pdb) # TODO: AGA intercept the result to the report



    # ================= CUSTOM DATABASE COMPARISON =================
    if args.own_db:
        # Check if Foldseek DB already exists for the given folder, if not create it
        custom_db_foldername = filesystem.get_db_name_from_path(args.own_db)
        custom_db_path = f"{args.out_dir}/db/custom_db/{custom_db_foldername}"

        if not foldseek.foldseek_db_exists(custom_db_path):
            print(f"[INFO] Creating Foldseek DB from custom structures in '{args.own_db}'...")
            custom_db = foldseek.create_own_db(args.own_db, f"{args.out_dir}/db/custom_db/{custom_db_foldername}/{custom_db_foldername}") # the last after "/" is the name used for filenames of db
        else:
            print(f"[INFO] Foldseek DB already exists for '{args.own_db}'. Using existing DB.")
            custom_db = f"{args.out_dir}/db/custom_db/{custom_db_foldername}/{custom_db_foldername}"
        
        # Perform comparison against the custom database
        custom_results_df, custom_structure_ids = compare_structures.foldseek_search_db(mut_pdb, custom_db, k=args.top_k) # TODO: AGA intercept the result to the report



    # ================= SUPPORTED DATABASE COMPARISON =================
    if not args.exclude_pdb:
        # Check if Foldseek DB already exists for the given folder, if not create it
        pdb_db_path = f"{args.out_dir}/db/supported_db/pdb_db"

        if not foldseek.foldseek_db_exists(pdb_db_path):
            print(f"[INFO] Creating Foldseek DB from PDB structures...")
            pdb_db = foldseek.create_supported_foldseek_db("PDB", f"{args.out_dir}/db/supported_db/pdb_db", "pdb_db") # the name of the files here are separete because of tmp folder
        else:
            print(f"[INFO] Foldseek DB already exists for PDB. Using existing DB.")
            pdb_db = f"{args.out_dir}/db/supported_db/pdb_db/pdb_db" # the name of the files here are separete because of tmp folder

        # Perform comparison against the PDB database
        pdb_results_df, pdb_structure_ids = compare_structures.foldseek_search_db(mut_pdb, pdb_db, k=args.top_k) # TODO: AGA intercept the result to the report
        pdb_db_structures = fetch.download_pdb_structures(pdb_structure_ids, f"{args.out_dir}/pdb_db_structures") # TODO: AGA intercept the result to the report

    if not args.exclude_af:
        # Check if Foldseek DB already exists for the given folder, if not create it
        af_db_path = f"{args.out_dir}/db/supported_db/af_db"

        if not foldseek.foldseek_db_exists(af_db_path):
            print(f"[INFO] Creating Foldseek DB from Alphafold/Swiss-Prot structures...")
            af_db = foldseek.create_supported_foldseek_db("Alphafold/Swiss-Prot", f"{args.out_dir}/db/supported_db/af_db", "af_db") # the name of the files here are separete because of tmp folder
        else:
            print(f"[INFO] Foldseek DB already exists for Alphafold/Swiss-Prot. Using existing DB.")
            af_db = f"{args.out_dir}/db/supported_db/af_db/af_db" # the name of the files here are separete because of tmp folder
        
        
        af_results_df, af_structure_ids = compare_structures.foldseek_search_db(mut_pdb, af_db, k=args.top_k) # TODO: AGA intercept the result to the report
        af_db_structures = fetch.download_af_structures(af_structure_ids, f"{args.out_dir}/af_db_structures") # TODO: AGA intercept the result to the report



    # ================= STATISTICS, REPORT GENERATION =================




if __name__ == "__main__":
    main()