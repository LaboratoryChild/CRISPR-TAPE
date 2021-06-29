# CRISPR-TAPE
A protein-centric CRISPR gRNA design tool for TArgeted Protein Engineering

# Motivation
Existing CRIPSR gRNA design tools target protein-coding regions within genomic loci and non-specifically target the entire input region of DNA. Current tools fail to consider proteomic-based applications, so CRISPR-TAPE has been developed to reduce the substantial time burden associated with manual curation of gRNA libraries and empower the proteomics community.

# Version
CRISPR-TAPE version 1.0.3

# LICENSE
Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License

# Quick Start

    pip install CRISPR-TAPE

# Technologies
The CRISPR-TAPE python scripts are compatible with Python 3.6 and later versions. Outside of the standard library, it makes use of the following packages: numpy version 1.18.0, pandas version 0.25.3.

# Troubleshooting
• All inputs are case sensitive and it is important to adhere to the requirements specified adjacent to each entry box

• gRNAs will be incorrect if exonic bases are not capitalised and intronic/untranslated bases are not lowercase in the genomic loci input

• The programme currently only recognises “A”, “T”, “C” and “G” bases and will filter out any other characters

• If no protospacer adjacent motif is specified, the programme will not run and no guide RNAs will be outputted

• The programme will stop running and a pop-up prompt will open if the inputted protein coding sequence and the exon sequences from the inputted genomic loci do not match

• Off target counts will display as “-1” if gRNAs are not found within the forward or reverse complement of the organism genome. If this occurs, ensure the genomic loci and genome originate from the same organism

• If some amino acids are missing from the output, this is because they cannot be target using the current information. If this is the case, add more upstream and downstream bases to the input requesting 100 bases upstream and downstream of the genomic loci

# Authors
CRISPR-TAPE was developed by Daniel Anderson, Henry Benns and Dr Matthew Child

# Copyright
Copyright (c) 2020 The Child Lab
