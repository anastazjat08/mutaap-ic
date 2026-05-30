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

### Step 4 - Test he installatio

```
mutaap --help
```

## Quick usage example

```
mutaap orig_fasta_path mut_fasta_path \
      --own_db own_structures_folder \
      --exclude_af \
      --top_k 20 \
      --out_dir all_results_path
```

## Repository structure
TBU

```
mutaap-ic/
└── mutaap-ic/
    ├── examples/                     # Example input files and usage workflows
    │   └── .gitkeep
    ├── images/                       # Project images (e.g., logo for README)
    │   └── logo.png
    ├── mutaapic/                     # Main Python package
    │   ├── analysis/                 # Structural analysis modules
    │   │   ├── __init__.py
    │   │   └── compare_structures.py # Foldseek/TM-align logic for structure comparison
    │   ├── io/                       # Input/output utilities (future expansion???)
    │   │   ├── __init__.py
    │   │   └── .gitkeep
    │   ├── orf/                      # ORF detection and sequence-level mutation logic
    │   │   ├── __init__.py
    │   │   └── gitkeep
    │   ├── reporting/                # Report generation (tables, summaries, visual outputs)
    │   │   ├── __init.py
    │   │   └── .gitkeep
    │   ├── structure/                # Structure prediction and handling
    │   │   ├── __init__.py
    │   │   └── predict_structure.py  # ESMFold prediction wrapper
    │   ├── utils/                    # General-purpose helper functions
    │   │   ├── __init__.py
    │   │   ├── fetch.py              # Downloading/fetching external structures
    │   │   ├── filesystem.py         # Path handling, directory management
    │   │   └── foldseek.py           # Foldseek database creation and search wrappers
    │   ├── validation/               # Input validation and sanity checks
    │   │   ├── __init__.py
    │   │   └── validate_sequence.py  # Sequence validation logic
    │   ├── __init__.py               
    │   └── cli.py                    # Command-line interface (`mutaap` entry point)
    ├── tests/                        # Unit tests for package modules
    │   └── .gitkeep
    ├── environment.yaml              # Conda environment specification
    ├── LICENSE                       # License information
    ├── pyproject.toml                # Build system, dependencies, package metadata
    └── README.md                     # Project documentation


```

