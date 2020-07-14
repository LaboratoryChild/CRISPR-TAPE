#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 23:38:29 2020

@author: danielanderson
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter 

from shared_functions import clean_inputs, get_codon_index, PAMposition, analyse_text, notes, pamcolumn, remove_pam, correct_distance

def context(search, aas): #Return the amino acid being target (marked by '*') and the amino acids immediately surrounding it
    x = int(search)
    if x == 0:
        surrounding= aas[x] + '*' + aas[x+1] + aas[x+2] + aas[x+3] +aas[x+4]
    else:
        if x == 1:
            surrounding= aas[x-1] + aas[x] + '*' + aas[x+1] + aas[x+2] +aas[x+3]
        else:
            if x == len(aas) - 2:
                surrounding= aas[x-4] + aas[x-3] +aas[x-2] + aas[x-1] + aas[x] + '*'
            else:
                if x == len(aas) - 3:
                    surrounding= aas[x-3] + aas[x-2] +aas[x-1] + aas[x] + '*' + aas[x+1]
                else:
                    surrounding= aas[x-2] + aas[x-1] +aas[x] + '*' + aas[x+1] + aas[x+2]
    return surrounding

def get_gRNA_sequence(amino_pos, guide_dict):
    guide_seq = ""
    for x in range(len(guide_dict)):
        if x == amino_pos:
            guide_seq = guide_dict[x]["Sequence"]
    return guide_seq

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    if 'A' in dna or 'T' in dna or 'C' in dna or 'G' in dna:
        comp = ''.join([complement[base] for base in dna[::-1]])
    else:
        comp = dna
    return comp

def get_strand(amino_pos, guide_dict):
    guide_seq = ""
    for x in range(len(guide_dict)):
        if x == amino_pos:
            guide_seq = guide_dict[x]["Strand"]
    return guide_seq

def get_distance(amino_pos, guide_dict):
    guide_seq = ""
    for x in range(len(guide_dict)):
        if x == amino_pos:
            guide_seq = guide_dict[x]["Distance"]
    return guide_seq

def get_gc(amino_pos, guide_dict):
    guide_seq = ""
    for x in range(len(guide_dict)):
        if x == amino_pos:
            guide_seq = guide_dict[x]['G/C Content (%)']
    return guide_seq

def list_search(the_list, orgen):
    the_count = []
    the_set = set(the_list)
    for x in tqdm(range(len(orgen))):
        if orgen[x:x+23] in the_set:
            the_count.append(orgen[x:x+23])
        else:
            pass
    return the_count

def get_count(fw, rv, counter):
    count = -1
    for key, value in counter.items():
        if fw == key or rv == key:
            count += int(value)
    return count



#def gRNA_base(start_aa, stop_aa, protein_dict):
    
    
def General_function(aa, motif, dna, cds, hundredup, hundreddown, orgen):

    dna_exon, dna, dna_edited, cds, orgen = clean_inputs(dna, cds, orgen, hundredup, hundreddown)

    if dna_exon == cds: #Confirm that the inputted CDS sequence matches the CDS identified by the tool. If not the programme stops running.
        print("Inputted CDS and concatenated exons match")
    else:
       print("INPUTTED CDS AND CONCATENATED EXONS DO NOT MATCH")
     
    #if len(aa) == 1:
    protein_dict = get_codon_index(dna, cds)    
    amino_position = []
    amino_acids = []
    
    for key, value in protein_dict.items():
        if value["Amino Acid"] == aa:
            amino_position.append(key)
        amino_acids.append(value["Amino Acid"])
        
    if len(amino_position) == 0:
        raise AttributeError("The amino acid ", aa," is not present in this sequence.")   
        
    print("Searching for guides...")
    
    guide_positons, gRNA_list, guide_strands = PAMposition(dna_edited, motif)#Identify guides in the inputted genomic loci
    gRNA = pd.DataFrame(gRNA_list, columns= ['full gRNA Sequence']) #Convert the gRNA list into a dataframe
    gRNA['Position'] = guide_positons #Add a column of guide RNA cut site positions to the dataframe
    gRNA['Strand'] = guide_strands #Add a column to specify these guides are on the forward strand relative to the inputted genomic loci.
    gRNA["G/C Content (%)"] = gRNA['full gRNA Sequence'].apply(analyse_text)
    gRNA = gRNA[gRNA['G/C Content (%)'] <= 75] #Remove guide RNAs with a G/C percentage above 75%
    gRNA = gRNA.sort_values(by=['Position']).reset_index(drop=True) #Sort the guide RNAs by their position in the genomic loci
    
    start_base= []
    stop_base = []
    for key, value in protein_dict.items():
        if key in amino_position:
            start_base.append(value["base_1"][0])
            stop_base.append(value["base_3"][0])
            
    guide_array = np.array(gRNA['Position'])
    start_base_array = np.array(start_base)
    stop_base_array = np.array(stop_base)
    
    start_distances = np.subtract(guide_array,start_base_array[:,None])
    stop_distances = np.subtract(guide_array,stop_base_array[:,None])
# =============================================================================
#     elif len(aa)>1:
#         first_base_array = np.array(first)
#         first_distances = np.subtract(guide_array,first_base_array[:,None])
#         first_distances = np.where(first_distances > 1, 2, first_distances)
#         
#         final_base_array = np.array(final)
#         final_distances = np.subtract(guide_array,final_base_array[:,None])
#         final_distances = np.where(first_distances < -1, -2, final_distances)
#     else:
#         raise ValueError("Invalid aa entry")
# =============================================================================
        
    upstream_sequences = []
    downstream_sequences = []
    
    for x in start_distances:
        index = ''  
        for y in range(len(x)-1):
            if x[y]<=0 and x[y+1]>0:
                index = y
                up_index = y
                continue       
        if not index == '':
            upstream_sequences.append({"Position": up_index, "Sequence": gRNA['full gRNA Sequence'][up_index], "Distance": x[up_index], "Strand": gRNA['Strand'][up_index], "G/C Content (%)": gRNA["G/C Content (%)"][up_index]})
        elif index == '' and x[len(x)-1]<=0:
            upstream_sequences.append({"Position": len(x)-1, "Sequence": gRNA['full gRNA Sequence'][len(x)-1], "Distance": x[len(x)-1], "Strand": gRNA['Strand'][len(x)-1], "G/C Content (%)": gRNA["G/C Content (%)"][len(x)-1]})
        elif index == '' and x[len(x)-1]>0:
            upstream_sequences.append({"Position": "", "Sequence" : "No 5' guide could be identified", "Distance": "", "Strand": "", "G/C Content (%)": 0})
        else:
            raise AttributeError("No guide RNAs identified in sequence")
    
    for x in stop_distances:
        index = ''  
        for y in range(len(x)-1):
            if x[y]<=0 and x[y+1]>0:
                index = y
                down_index = y+1
                continue       
        if not index == '':
            downstream_sequences.append({"Position": down_index, "Sequence": gRNA['full gRNA Sequence'][down_index], "Distance": x[down_index], "Strand": gRNA['Strand'][down_index], "G/C Content (%)": gRNA["G/C Content (%)"][down_index]})
        elif index == '' and x[len(x)-1]<=0:
            downstream_sequences.append({"Position": "", "Sequence": "No 3' guide could be identified", "Distance": "", "Strand": "", "G/C Content (%)": 0})
        elif index == '' and x[len(x)-1]>0:
            downstream_sequences.append({"Position": 0, "Sequence" : gRNA['full gRNA Sequence'][0], "Distance": x[0], "Strand": gRNA['Strand'][0], "G/C Content (%)": gRNA["G/C Content (%)"][0]})
        else:
            raise AttributeError("No guide RNAs identified in sequence")
    
# =============================================================================
#     for x in stop_distances:
#         index = ''  
#         for y in range(len(x)-1):
#             if x[y]<=0 and x[y+1]>0:
#                 index = y
#                 up_index = y
#                 down_index = y+1
#                 continue       
#         if not index == '':
#             upstream_sequences.append({"Position": up_index, "Sequence": gRNA['full gRNA Sequence'][up_index], "Distance": x[up_index], "Strand": gRNA['Strand'][up_index], "G/C Content (%)": gRNA["G/C Content (%)"][up_index]})
#             downstream_sequences.append({"Position": down_index, "Sequence": gRNA['full gRNA Sequence'][down_index], "Distance": x[down_index], "Strand": gRNA['Strand'][down_index], "G/C Content (%)": gRNA["G/C Content (%)"][down_index]})
#         elif index == '' and x[len(x)-1]<=0:
#             upstream_sequences.append({"Position": len(x)-1, "Sequence": gRNA['full gRNA Sequence'][len(x)-1], "Distance": x[len(x)-1], "Strand": gRNA['Strand'][len(x)-1], "G/C Content (%)": gRNA["G/C Content (%)"][len(x)-1]})
#             downstream_sequences.append({"Position": "", "Sequence": "No 3' guide could be identified", "Distance": "", "Strand": "", "G/C Content (%)": 0})
#         elif index == '' and x[len(x)-1]>0:
#             upstream_sequences.append({"Position": "", "Sequence" : "No 5' guide could be identified", "Distance": "", "Strand": "", "G/C Content (%)": 0})
#             downstream_sequences.append({"Position": 0, "Sequence" : gRNA['full gRNA Sequence'][0], "Distance": x[0], "Strand": gRNA['Strand'][0], "G/C Content (%)": gRNA["G/C Content (%)"][0]})
#         else:
#             raise AttributeError("No guide RNAs identified in sequence")
#     
#     
# =============================================================================
    
    guides = pd.DataFrame(amino_position, columns= ['Amino Acid Position']) 
    guides["Context"] = guides.apply(lambda row: context(row["Amino Acid Position"], amino_acids), axis = 1)     
    guides["Amino Acid Position"] = guides["Amino Acid Position"] + 1
    guides["5' gRNA Sequence"] = guides.apply(lambda row: get_gRNA_sequence(row.name, upstream_sequences), axis = 1)  
    guides["5' gRNA Sequence RC"] = guides["5' gRNA Sequence"].apply(reverse_complement)
    guides["3' gRNA Sequence"] = guides.apply(lambda row: get_gRNA_sequence(row.name, downstream_sequences), axis = 1)
    guides["3' gRNA Sequence RC"] = guides["3' gRNA Sequence"].apply(reverse_complement)
    guides["5' gRNA Strand"] = guides.apply(lambda row: get_strand(row.name, upstream_sequences), axis = 1)   
    guides["3' gRNA Strand"] = guides.apply(lambda row: get_strand(row.name, downstream_sequences), axis = 1)      
    guides["Distance of 5' Cut Site from Amino Acid (bp)"] = guides.apply(lambda row: get_distance(row.name, upstream_sequences), axis = 1)  
    guides["Distance of 5' Cut Site from Amino Acid (bp)"] = guides.apply(lambda row: correct_distance(row["Distance of 5' Cut Site from Amino Acid (bp)"], row["5' gRNA Strand"]), axis = 1)
    guides["Distance of 3' Cut Site from Amino Acid (bp)"] = guides.apply(lambda row: get_distance(row.name, downstream_sequences), axis = 1)   
    guides["Distance of 3' Cut Site from Amino Acid (bp)"] = guides.apply(lambda row: correct_distance(row["Distance of 3' Cut Site from Amino Acid (bp)"], row["3' gRNA Strand"]), axis = 1)      
    guides["5' gRNA G/C Content (%)"] =  guides.apply(lambda row: get_gc(row.name, upstream_sequences), axis = 1)  
    guides["3' gRNA G/C Content (%)"] =  guides.apply(lambda row: get_gc(row.name, downstream_sequences), axis = 1)  
    guides["5' Notes"] = guides.apply(lambda row: notes(row["5' gRNA Sequence"], row["5' gRNA G/C Content (%)"]), axis = 1)  
    guides["3' Notes"] = guides.apply(lambda row: notes(row["3' gRNA Sequence"], row["3' gRNA G/C Content (%)"]), axis = 1)  

    print("Counting off targets...")
    
    ls = list(guides["5' gRNA Sequence"]) + list(guides["5' gRNA Sequence RC"]) + list(guides["3' gRNA Sequence"]) + list(guides["3' gRNA Sequence RC"])
    count_list = list_search(ls, orgen)
    counter = Counter(count_list)
        
    guides["5' gRNA Off Target Count"] = guides.apply(lambda row: get_count(row["5' gRNA Sequence"], row["5' gRNA Sequence RC"], counter), axis=1)
    guides["3' gRNA Off Target Count"] = guides.apply(lambda row: get_count(row["3' gRNA Sequence"], row["3' gRNA Sequence RC"], counter), axis=1)
    
    guides["5' PAM"] = guides.apply(lambda row: pamcolumn(row["5' gRNA Sequence"], row["5' gRNA Strand"], row["5' gRNA Sequence RC"], motif), axis =1) 
    guides["3' PAM"] = guides.apply(lambda row: pamcolumn(row["3' gRNA Sequence"], row["3' gRNA Strand"], row["3' gRNA Sequence RC"], motif), axis =1) 
    guides["5' gRNA Sequence"] = guides.apply(lambda row: remove_pam(row["5' gRNA Sequence"], row["5' gRNA Strand"], row["5' gRNA Sequence RC"], motif), axis =1)
    guides["3' gRNA Sequence"] = guides.apply(lambda row: remove_pam(row["3' gRNA Sequence"], row["3' gRNA Strand"], row["3' gRNA Sequence RC"], motif), axis =1)

    guides[' '] = '' #Empty column in dataframe
     
    guides = guides[["Amino Acid Position", "Context", ' ', "5' gRNA Sequence", "5' PAM", "5' gRNA Strand", "5' gRNA G/C Content (%)", "Distance of 5' Cut Site from Amino Acid (bp)", "5' Notes", "5' gRNA Off Target Count", ' ', "3' gRNA Sequence", "3' PAM", "3' gRNA Strand", "3' gRNA G/C Content (%)", "Distance of 3' Cut Site from Amino Acid (bp)", "3' Notes", "3' gRNA Off Target Count"]] #Reorganise guide dataframe columns for output
    guides.columns = ["Amino Acid Position", "Adjacent amino acids", '                        ', "gRNA Sequence 5' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count", '                        ', "gRNA Sequence 3' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count"] #Rename guide dataframe columns for output
    
    return guides

with open("INSERT_ORGANISM_GENOME_HERE.txt", "r") as f:
    orgen = f.read()

with open("Gene_CDS.txt", "r") as g:
    cds = g.read()

with open("Genomic_loci.txt", "r") as h:
    dna = h.read()  

import time

s = time.time()
new_guides = General_function("M", "NGG", dna, cds, "", "", orgen)
e = time.time()
t1 = str(e-s) + " seconds"
print(t1)