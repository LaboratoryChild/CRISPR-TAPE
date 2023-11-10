# CRISPR-TAPE
A protein-centric CRISPR gRNA design tool for TArgeted Protein Engineering

# Motivation
Existing CRIPSR gRNA design tools target protein-coding regions within genomic loci and non-specifically target the entire input region of DNA. Current tools fail to consider proteomic-based applications, so CRISPR-TAPE has been developed to reduce the substantial time burden associated with manual curation of gRNA libraries and empower the proteomics community.

# Version
CRISPR-TAPE version 1.0.6

# LICENSE
Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License

# Quick Start
To install CRISPR-TAPE from source, run:
```
git clone https://github.com/LaboratoryChild/CRISPR-TAPE.git
```
To install the CRISPR-TAPE dependencies, run:
```
cd CRISPR-TAPE
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
```
Alternatively, CRISPR-TAPE is available via PyPI:
```
pip install CRISPR-TAPE
```
To use the general amino acid targeting function, add the following to the top of your Python script:
```
from CRISPR_TAPE import general_function
```
To use the specific amino acid targeting function, add the following to the top of your Python script:
```
from CRISPR_TAPE import specific_function
```

# Linux users
Linux users must install python3-tk via apt-get prior to running CRISPR-TAPE for the first time by running:
```
sudo apt-get install python3-tk
```

## Executable
To build the CRISPR-TAPE executable, run:
```
python3 build_executable.py
```
The executable can be found in the generated ```dist``` directory. The genome assembly file used to check for mismatches must be in the same directory as the executable.

## Command line
CRISPSR-TAPE may be run through the GUI or run from the command line. To use the CRISPR-TAPE GUI, run:
```
python3 tape_main-runner.py
```
To run the CRISPR-TAPE specific targeting function through the command line:
```
python3 tape_specific-runner.py --loci <genomic loci file> --cds <coding sequence file> --genome <organism genomic sequence file> --spec-amino <residue position> --motif <protospacer adjacent motif> --distance <maximum distance of cut site from residue> --output <output filename> --up <additional bases upstream of loci> --down <additional bases downstream of loci>
```
To run the CRISPR-TAPE general targeting function through the command line:
```
python3 tape_general-runner.py --loci <genomic loci file> --cds <coding sequence file> --genome <organism genomic sequence file> --aa <amino acide short letter code> --motif <protospacer adjacent motif> --output <output filename> --up <additional bases upstream of loci> --down <additional bases downstream of loci>
```
For more detailed documentation, see the [docs](https://github.com/LaboratoryChild/CRISPR-TAPE/tree/devel/docs)

# Testing
CRIPSR-TAPE tests can be conducted by running:
```
cd test
python3 run_test.py
```

# Troubleshooting
• All inputs are case sensitive and it is important to adhere to the requirements specified adjacent to each entry box

• gRNAs will be incorrect if exonic bases are not capitalised and intronic/untranslated bases are not lowercase in the genomic loci input

• The programme currently only recognises “A”, “T”, “C” and “G” bases and will filter out any other characters

• If no protospacer adjacent motif is specified, the programme will not run and no guide RNAs will be outputted

• The programme will stop running and a pop-up prompt will open if the inputted protein coding sequence and the exon sequences from the inputted genomic loci do not match

• Off target counts will display as “-1” if gRNAs are not found within the forward or reverse complement of the organism genome. If this occurs, ensure the genomic loci and genome originate from the same organism

• If some amino acids are missing from the output, this is because they cannot be target using the current information. If this is the case, add more upstream and downstream bases to the input requesting 100 bases upstream and downstream of the genomic loci

• Please submit GitHub issues for further assistance or to report bugs.

# Technologies
The CRISPR-TAPE python scripts are compatible with Python 3.6 and later versions. Outside of the standard library, it makes use of the following packages: numpy v1.24.4, pandas v2.0.3, pyinstaller v6.1.0 and tqdm v4.66.1.

# Authors
CRISPR-TAPE was developed by Daniel Anderson, Dr Henry Benns and Dr Matthew Child

# Copyright
Copyright (c) 2020 The Child Lab
