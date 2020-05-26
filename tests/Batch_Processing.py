#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 19:08:09 2020

@author: danielanderson
"""


from Specific_function import Specific_function
from General_function import General_function
from tkinter import *
import pandas as pd
import numpy as np

data = pd.read_csv("Template.csv")
o = open("INSERT_ORGANISM_GENOME_HERE.txt", 'r')
data = data.replace(np.nan, '', regex=True)
orgen = str(o.read())

for x in range(len(data)):
    filename = str(data['Gene Identifier'][x])
    cds = str(data['Genomic Loci Sequence'][x])
    dna = str(data['Coding Sequence'][x])
    distance = int(data['Distance'][x])
    motif = str(data['Motif'][x])
    spec_amino = data['Residue Position'][x]
    if spec_amino == '':
        pass
    else:
        spec_amino = int(spec_amino)
    aa = str(data['Amino Acid Type (Short Code)'][x])
    if aa == '':
        pass
    else:
        aa = str(aa)
    if distance == '0' or distance == "" or distance == " ":
        distance = 1000000
    else:
        distance = int(distance)
    hundredup = str(data["5' UTR"][x])
    hundreddown = str(data["3' UTR"][x])
    if aa == '':
        guides = Specific_function(spec_amino, distance, motif, cds, dna, hundredup, hundreddown, orgen)
    if spec_amino == '':
        guides = General_function(aa, motif, cds, dna, hundredup, hundreddown, orgen)
    csvfilename = str(filename + ".csv") #Add '.csv' to the specified filename 
    guides.to_csv(csvfilename, index=False)

print("done")
