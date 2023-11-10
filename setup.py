#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# read the contents of your requirements file
with open("requirements.txt", "r") as fh:
    requirements = fh.read().splitlines()

setuptools.setup(
    name="CRISPR-TAPE",
    version="1.0.6",
    author="Daniel Anderson",
    author_email="danielanderson1@hotmail.com",
    license= "Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License",
    description="A protein-centric CRISPR gRNA design tool for TArgeted Protein Engineering",
    long_description="Existing CRIPSR gRNA design tools target protein-coding regions within genomic loci and non-specifically target the entire input region of DNA. Current tools fail to consider proteomic-based applications, so CRISPR-TAPE has been developed to reduce the substantial time burden associated with manual curation of gRNA libraries and empower the proteomics community." ,
    long_description_content_type="text/markdown",
    url="https://github.com/LaboratoryChild/CRISPR-TAPE",
    install_requires=requirements,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',  # Update this to the minimum supported Python version
)
