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

from mutaapic.utils.read_files import read_fasta, read_txt
from mutaapic.orf.validation import isit_orf, validate_no_internal_stop
from mutaapic.orf.alignment import automatic_alignment, normalize_changes
from mutaapic.orf.protein_analysis import compare_proteins, modified_analysis
from mutaapic.analysis.inputs_summary import generate_input_summary


def main():
    parser = ap.ArgumentParser(description="MutAAP-IC: A tool for predicting and analyzing the effects of mutations on protein structure.")
    # parser.add_argument("fasta_orig", type=str, help="The FASTA file containing the original nucleotide sequence to analyze.")
    parser.add_argument("original_file", 
                        help="FASTA file with original ORF - DNA sequence")
    # parser.add_argument("fasta_mut", type=str, help="The FASTA file containing the mutated nucleotide sequence to analyze.")
    parser.add_argument("modified_file", 
                        help="FASTA file with modified ORF - DNA sequence")
    parser.add_argument("changes_file", nargs="?", default=None,
                        help="user-provided modifications in given ORF, format should be:" \
                        "one change in one line, eg.:\n" \
                        "8\n" \
                        "28-32")    
    parser.add_argument("--auto-alignment", action="store_true",
                        help="do automatic alignment")
    parser.add_argument("--output", default="report", 
                        help="output directory")
    parser.add_argument("--table", type=int, default=11, 
                        help="NCBI translation code table, default is bacterial")
    parser.add_argument("--db_path", type=str, default="./mutaap_db", help="Path where the default Foldseek databases are stored and downloaded if missing.")
    parser.add_argument("--custom_db", type=str, help="Path to a custom database of structures for comparison (optional).")
    parser.add_argument("--exclude_pdb", action="store_true", help="Whether to exclude comparison with structures from PDB DB supported by Foldseek (default: False).")
    parser.add_argument("--exclude_af", action="store_true", help="Whether to exclude comparison with structures from Alphafold DB supported by Foldseek (default: False).")
    parser.add_argument("--top_k", type=int, default=10, help="Number of top structures found by Foldseek in specific DB to take into statistical analysis (default: 10).")
    parser.add_argument("--out_dir", type=str, default="./mutaap_results", help="Directory to save results (default: mutaap_results).")

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

    def _resolve_db_root(db_root: str):
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


    orig_dna = read_fasta(args.original_file)
    mut_dna = read_fasta(args.modified_file)

    orig_sequence = isit_orf(orig_dna, args.table)

    mut_sequence, frameshift = (
        modified_analysis(
            mut_dna, orig_dna, mut_dna, args.table))
    
    validate_no_internal_stop(mut_sequence)


    if args.auto_alignment:
        nt_changes = normalize_changes(
            automatic_alignment(
                orig_dna, mut_dna))
        

    elif args.changes_file:
        nt_changes = read_txt(args.changes_file)
    
    else:
        raise ValueError("Provide --auto-alignment for automatic alignment or your own file")

    mutations = compare_proteins(orig_sequence, mut_sequence)

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



    # ================= STATISTICS, REPORT GENERATION =================

    # Generate a report HTML summarizing results
    try:
        report_path = report.generate_report(
            args.out_dir,
            orig_pdb,
            mut_pdb,
            comparison_results,
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