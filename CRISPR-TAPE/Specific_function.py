#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identify guides surrounding a specific residue position within a maximum distance range
"""

import pandas as pd
import numpy as np
from tqdm import tqdm

from shared_functions import clean_inputs, get_codon_index, PAMposition, analyse_text, notes, pamcolumn, remove_pam, correct_distance

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
    
def off_target(search, reverse_search, orgen):
     if not search == "":
         term = [search,reverse_search]
         count = sum(orgen.count(i) for i in term)- 1
     else:
         count = ""
     return count

def Specific_function(spec_amino, distance, motif, dna, cds, hundredup, hundreddown, orgen):    
    
    dna_exon, dna, dna_edited, cds, orgen = clean_inputs(dna, cds, orgen, hundredup, hundreddown)
    
    if dna_exon == cds: #Confirm that the inputted CDS sequence matches the CDS identified by the tool. If not the programme stops running.
        print("Inputted CDS and concatenated exons match")
    else:
       print("INPUTTED CDS AND CONCATENATED EXONS DO NOT MATCH")
      
    protein_dict = get_codon_index(dna, cds) 
    
    for key, value in protein_dict[spec_amino - 1].items():
        if key == "base_2":   
            middle = value[0]
        if key == "Amino Acid":
            selectedaa = value
    
    selectionmade = selectedaa + '-' + str(spec_amino)
    
    print('The amino acid you have selected is ' + selectionmade) #Confirm the selection is the correct amino acid
    
   # generate lists of all the potential guides within the genomic loci and a list of all guides on the reverse complement.
    print("Searching for guides...")
    guide_positons, gRNA_list, guide_strands = PAMposition(dna_edited, motif)#Identify PAMs in the inputted genomic loci
    
    # convert the lists into dataframes
    gRNA = pd.DataFrame(gRNA_list, columns= ['full gRNA Sequence']) #Convert the gRNA list into a dataframe
    gRNA['Position'] = guide_positons #Add a column of guide RNA cut site positions to the dataframe
    gRNA['Strand'] = guide_strands #Add a column to specify these guides are on the forward strand relative to the inputted genomic loci.
    gRNA = gRNA.sort_values(by=['Position']).reset_index(drop=True) #Sort the guide RNAs by their position in the genomic loci

    gRNA["Distance from Amino Acid (bp)"] = np.array(gRNA["Position"]) - middle #Determine guide distance from base at the centre of the codon
    
    upperguides = "No guides within distance range specified"
    downerguides = "No guides within distance range specified"
    
    if not len(gRNA) == 0:
        print(str(len(gRNA)) + " guides found")
        slice_index = ""
        for dis in range(len(gRNA)-1):
            if gRNA["Distance from Amino Acid (bp)"][dis] < 0 and gRNA["Distance from Amino Acid (bp)"][dis+1] > 0:
                slice_index = dis
        if not slice_index == "":
            upperguides = gRNA.iloc[:slice_index+1,]
            downerguides = gRNA.iloc[slice_index+1:,]
        elif gRNA["Distance from Amino Acid (bp)"][dis] < 0:
            upperguides = gRNA 
        else:
            downerguides = gRNA 
    else:
        raise AttributeError("No guide RNAs identified in sequence")
    
    tqdm.pandas()
    
    noguides = pd.DataFrame(columns=["full gRNA Sequence", "Position", "Strand", "Distance from Amino Acid (bp)", "Reverse complement", "PAM", "gRNA Sequence", "G/C Content (%)", "Notes", 'Off Target Count']) #Generate a new amino acid for the amino acid target information
    noguides = noguides.append({"full gRNA Sequence": "", "Position" : "", "Strand":"","Distance from Amino Acid (bp)":"","Reverse complement": "", "PAM": "", "gRNA Sequence": "No guides within distance range specified", "G/C Content (%)":"","Notes":"", 'Off Target Count': ""}, ignore_index=True) #Generate a new amino acid for the amino acid target information

    if not isinstance(upperguides, str):
        upperguides = upperguides[upperguides["Distance from Amino Acid (bp)"] >= (-distance)]  #Remove guides over the inputted maximum guide distance
        upperguides = upperguides.sort_values(by = ['Distance from Amino Acid (bp)']) #Arrange guide RNAs by their distance from the amino acid
        upperguides['Distance from Amino Acid (bp)'] = upperguides.apply(lambda row: correct_distance(int(row['Distance from Amino Acid (bp)']), row['Strand']), axis =1) + 1
        upperguides["Reverse complement"] = upperguides.apply(lambda row: reverse_complement(row['full gRNA Sequence']), axis =1)
        upperguides["PAM"] = upperguides.apply(lambda row: pamcolumn(row['full gRNA Sequence'], row['Strand'], row["Reverse complement"], motif), axis =1) 
        upperguides["gRNA Sequence"] = upperguides.apply(lambda row: remove_pam(row['full gRNA Sequence'], row['Strand'], row["Reverse complement"], motif), axis =1)
        upperguides['G/C Content (%)'] = upperguides.apply(lambda row: analyse_text(row['full gRNA Sequence']), axis =1) #Calculate GC percentage
        upperguides["Notes"] = upperguides.apply(lambda row: notes(row["full gRNA Sequence"], row["G/C Content (%)"]), axis=1) #Output notes of key guide RNA characterstics to a new column
        upperguides = upperguides.reset_index(drop=True)
    else:
        upperguides = noguides
    
    amino_acid = pd.DataFrame(columns=["full gRNA Sequence", "Position", "Strand", "Distance from Amino Acid (bp)", "Reverse complement", "PAM", "gRNA Sequence", "G/C Content (%)", "Notes", 'Off Target Count']) #Generate a new amino acid for the amino acid target information
    amino_acid = amino_acid.append({"full gRNA Sequence": "", "Position" : "", "Strand":"","Distance from Amino Acid (bp)":"","Reverse complement": "", "PAM": "", "gRNA Sequence": selectionmade, "G/C Content (%)":"","Notes":"", 'Off Target Count': ""}, ignore_index=True) #Generate a new amino acid for the amino acid target information

    if not isinstance(downerguides, str):
        downerguides = downerguides[downerguides["Distance from Amino Acid (bp)"] <= distance]
        downerguides = downerguides.sort_values(by=['Distance from Amino Acid (bp)']) #Arrange guide RNAs by their distance from the amino acid
        downerguides['Distance from Amino Acid (bp)'] = downerguides.apply(lambda row: correct_distance(int(row['Distance from Amino Acid (bp)']),row['Strand']), axis =1) - 1
        downerguides["Reverse complement"] = downerguides.apply(lambda row: reverse_complement(row['full gRNA Sequence']), axis =1)
        downerguides["PAM"] = downerguides.apply(lambda row: pamcolumn(row['full gRNA Sequence'], row['Strand'], row["Reverse complement"], motif), axis =1) 
        downerguides["gRNA Sequence"] = downerguides.apply(lambda row: remove_pam(row['full gRNA Sequence'], row['Strand'], row["Reverse complement"], motif), axis =1)
        downerguides['G/C Content (%)'] = downerguides.apply(lambda row: analyse_text(row['full gRNA Sequence']), axis =1) #Calculate GC percentage
        downerguides["Notes"] = downerguides.apply(lambda row: notes(row["full gRNA Sequence"], row["G/C Content (%)"]), axis=1) #Output notes of key guide RNA characterstics to a new column
        downerguides = downerguides.reset_index(drop=True)
    else:
        downerguides = noguides
    
    guides = pd.concat([upperguides, amino_acid, downerguides])

    print("Counting off targets...")

    guides['Off Target Count'] = guides.progress_apply(lambda row: off_target(row['full gRNA Sequence'], row["Reverse complement"], orgen), axis=1)
    
    guides = guides[['Distance from Amino Acid (bp)', 'gRNA Sequence', 'PAM', 'Strand', 'G/C Content (%)', 'Off Target Count', 'Notes']] #Reorganise the guide RNA dataframe
    
    return guides
