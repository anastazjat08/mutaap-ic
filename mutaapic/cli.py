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
import mutaapic.utils.filesystem as filesystem
import mutaapic.utils.foldseek as foldseek
import mutaapic.utils.fetch as fetch
import mutaapic.reporting.report as report
from mutaapic.function.predict_function import predict_function, compare_function

from mutaapic.utils.read_files import read_fasta, read_txt
from mutaapic.orf.validate_sequence import isit_orf, validate_no_internal_stop
from mutaapic.orf.alignment import automatic_alignment, extract_nt_mismatches, merge_adjacent_changes
from mutaapic.analysis.aa_sequence_analysis import compare_proteins, modified_analysis, merge_adjacent_aa_changes
from mutaapic.analysis.input_summary import generate_input_summary


def main():
    parser = ap.ArgumentParser(description="MutAAP-IC: A tool for predicting and analyzing the effects of mutations on protein structure.")

    parser.add_argument("original_fasta", help="FASTA file containing the original nucleotide sequence to analyze.")
    parser.add_argument("modified_fasta", help="FASTA file containing the mutated nucleotide sequence to analyze.")

    parser.add_argument("changes_file", nargs="?", default=None,
                        help="Text file with user-provided ORF modifications. Expected format:" \
                        "one change in one line, eg.:\n" \
                        "8\n" \
                        "28-32\n" \
                        "40\tA\tG\tmismatch\n" \
                        "53-55\tAGGCTT\t------\tdeletion")
    
    # parser.add_argument("--auto-alignment", action="store_true", help="do automatic alignment") # DO WYWALENIA
    # parser.add_argument("--output", default="report", help="output directory") # DO WYWALENIAA

    parser.add_argument("--table", type=int, default=1, help="NCBI translation code table, default NCBI Table 1.")

    parser.add_argument("--db_path", type=str, default="./mutaap_db", help="Path where the default Foldseek databases are stored and downloaded if missing (default: mutaap_db).")
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
    orig_dna = read_fasta(args.original_fasta)
    mut_dna = read_fasta(args.modified_fasta)

    orig_sequence = isit_orf(orig_dna, args.table)

    mut_sequence, frameshift = modified_analysis(orig_dna, mut_dna, args.table)
    validate_no_internal_stop(mut_sequence)


    # if args.auto_alignment:
    #     nt_changes = normalize_changes(automatic_alignment(orig_dna, mut_dna))

    if args.changes_file is None:
        print("[INFO] No changes file provided. Running automatic mutation detection.")
        nt_changes = merge_adjacent_changes(automatic_alignment(orig_dna, mut_dna))
    elif args.changes_file:
        nt_changes = read_txt(args.changes_file)
    mutations = merge_adjacent_aa_changes(compare_proteins(orig_sequence, mut_sequence))

    print(orig_sequence)
    print(mut_sequence)
    
    generate_input_summary(
        output_dir=args.out_dir,
        mutations=mutations,
        nt_changes=nt_changes,
        original_protein=orig_sequence,
        modified_protein=mut_sequence,
        frameshift=frameshift
    )

    def clean_for_esm(protein):
        """
        removes stop codon from the end
        """
        return protein[:-1] if protein.endswith("*") else protein
    
    orig_for_esm = clean_for_esm(orig_sequence)
    mut_for_esm = clean_for_esm(mut_sequence)
    

    # ================= STRUCTURE PREDICTION =================
    orig_pdb = predict_structure.predictESM('orig', orig_for_esm, args.out_dir)
    time.sleep(5)
    mut_pdb = predict_structure.predictESM('mut', mut_for_esm, args.out_dir)

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