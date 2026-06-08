<p align="center">
  <img src="./images/logo.png" alt="MutAAP-IC logo" width=50%/>
</p>

# MutAAP-IC 
MutAAP-IC (_Mutational Amino Acid Position - Impact Check_) is a Python package designed to evaluate whether a nucleotide mutation affects the encoded amino acid and, consequently, the resulting protein. The tool performs a structured impact check by mapping nucleotide substitutions to codons, determining their effect on amino acids, and identifying potential consequences for the protein structure.

# Table of Contents
- [MutAAP-IC](#mutaap-ic)
- [Installation](#installation)
  - [Requirements](#requirements)
  - [Steps](#step-1---clone-the-repository)
- [Input format](#input-format)
- [Quick usage example](#quick-usage-example)
- [Repository structure](#repository-structure)

## Installation
### Requirements
- Linux or WSL2
- Conda (Miniconda/Anaconda)
- Foldseek installed automatically via environment
- Python >= 3.10

### Step 1 - Clone the repository

```
git clone https://github.com/anastazjat08/mutaap-ic.git
cd mutaap-ic
```

### Step 2 - Create the conda environment

```
conda env create -f environment.yaml
conda activate mutaap_env
```

### Step 3 - Install the MutAAP-IC package
```
pip install .
```
This installs the CLI command:
`mutaap`

### Step 4 - Test the installation

```
mutaap --help
```

## Input format
The input nucleotide sequence should be provided in a FASTA file.
Needs to be an ORF:
- start with START codon
- end with STOP codon
- be divisible by 3
- does not have STOP codons inside
- will give maximum of 400 amino acid length 



## Quick usage examples

```
mutaap orig_fasta_path mut_fasta_path


mutaap orig_fasta_path mut_fasta_path \
      --custom_db custom_structures_folder \
      --exclude_af \
      --top_k 20 \
      --table 11
      --out_dir all_results_path



```

## Repository structure

```

├── .gitignore
├── LICENSE                           # License information
├── Mutaap-ic.pdf                     # Project presentation
├── README.md                         # Project documentation
├── environment.yaml                  # Conda environment specification
├── examples                          # Example input and output files
├── images
│   └── logo.png                      # project logo
├── mutaapic                          # Main Python package
│   ├── __init__.py
│   ├── analysis                      # Structural analysis modules
│   │   ├── __init__.py
│   │   ├── aa_sequence_analysis.py
│   │   ├── compare_structures.py
│   │   └── input_summary.py
│   ├── cli.py                        # Command-line interface (`mutaap` entry point)
│   ├── function                      # Protein function module
|   |   ├── __init__.py
│   │   └── predict_function.py       
│   ├── orf                           # Module for validation input sequences
│   │   ├── __init__.py
│   │   ├── alignment.py
│   │   └── validate_sequence.py
│   ├── reporting                     # Module for reporting
│   │   ├── __init.py
│   │   └── report.py
│   ├── structure                     # Structure prediction module
│   │   ├── __init__.py
│   │   └── predict_structure.py
│   └── utils                         # General-purpose helper functions
│       ├── fetch.py                  # Downloading/fetching external structures
│       ├── filesystem.py             # Path handling, directory management
│       ├── foldseek.py               # Foldseek database creation and search wrappers
│       └── read_files.py             # Module for reading FASTA files
├── pyproject.toml                    # Build system, dependencies, package metadata
└── tests                             # Unit tests for package modules
    ├── conftest.py
    ├── test_alignment.py
    ├── test_cli.py
    ├── test_compare_structures.py
    ├── test_files.py
    ├── test_filesystem.py
    ├── test_foldseek.py
    ├── test_foldseek_search_db.py
    ├── test_input_summary.py
    ├── test_parser.py
    ├── test_predict_structure.py
    └── test_report.py


```

