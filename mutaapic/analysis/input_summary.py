import os
from typing import Any

'''
This module provides functions for generating a TSV summary report
of protein and nucleotide changes from given mutation data.
'''

def generate_input_summary(
        output_dir: str,
        mutations: list[dict[str, Any]],
        original_protein: str,
        modified_protein: str,
        frameshift: bool,
        nt_changes: list[str],
        ) -> str:
    """
    if directory doesn't exist, creates one: report 
    name of file with report: changes_report.tsv
    """
    os.makedirs(output_dir, exist_ok=True)

    report_path = os.path.join(
        output_dir, "changes_report.tsv"
    )

    with open(report_path, "w") as report:
        
        report.write("SUMMARY\n")

        report.write(f"original_protein_length:\t"
                     f"{len(original_protein)}\n")
        
        report.write(f"modified_protein_length:\t"
                     f"{len(modified_protein)}\n")

        report.write(f"frameshift:\t"
                     f"{frameshift}\n\n")
        

        # aminoacid mutations
        report.write("AMINO_ACID_CHANGES\n")

        report.write("aa_position_in_mut\t"
                     "original_aa\t"
                     "modified_aa\t"
                     "mutation_type\n"
                     )

        for mutation in mutations:
            report.write(
                f"{mutation['position']}\t"
                f"{mutation['original']}\t"
                f"{mutation['modified']}\t"
                f"{mutation['type']}\n\n"
            )

        # nucleotide changes
        report.write("\nNUCLEOTIDE_CHANGES\n")
        report.write("nt_position_in_mut\t"
                     "original_nt\t"
                     "modified_nt\t"
                     "change_type\n"
                     )
        for change in nt_changes:
            report.write(
                f"{change}\n")

    return report_path
  