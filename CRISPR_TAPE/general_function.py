#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identify gRNAs immediately adjacent to all amino acids of a specific type
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
import os
import sys

from CRISPR_TAPE.shared_functions import clean_inputs, get_codon_index, PAMposition, analyse_text, notes, pamcolumn, remove_pam, correct_distance

def context(search, aas, motif): #Return the amino acid being target (marked by '*') and the amino acids immediately surrounding it
    x = int(search)
    if len(motif) == 1:
        if x == 0:
            surrounding= aas[x] + '*' + aas[x+1] + aas[x+2] + aas[x+3] +aas[x+4]
        else:
            if x == 1:
                surrounding= aas[x-1] + aas[x] + '*' + aas[x+1] + aas[x+2] +aas[x+3]
            else:
                if x == len(aas) - 2 or x == len(aas) - 1:
                    surrounding= aas[x-4] + aas[x-3] +aas[x-2] + aas[x-1] + aas[x] + '*'
                else:
                    if x == len(aas) - 3:
                        surrounding= aas[x-3] + aas[x-2] +aas[x-1] + aas[x] + '*' + aas[x+1]
                    else:
                        surrounding= aas[x-2] + aas[x-1] +aas[x] + '*' + aas[x+1] + aas[x+2]
    else:
        amino_acids = "".join(aas)
        length = len(motif)%2

        if not len(motif) == 2:
            if length == 0:
                half = int(len(motif)/2)
            else:
                half = int((len(motif)+1)/2)
            surrounding = amino_acids[x - half + 1 : x + half]
        else:
            surrounding = amino_acids[x-1:x+1]

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

def general_function(aa,
                    motif,
                    dna,
                    cds,
                    hundredup,
                    hundreddown,
                    orgen):

    dna_exon, dna, dna_edited, cds, orgen = clean_inputs(dna, cds, orgen, hundredup, hundreddown)

    if len(aa) > 1:
        mid_aa = aa[round((len(aa) + 1)/2)-1]
    else:
        mid_aa = aa

    if dna_exon == cds: #Confirm that the inputted CDS sequence matches the CDS identified by the tool. If not the programme stops running.
        sys.stderr.write("\nInputted CDS and concatenated exons match\n")
    else:
       sys.stderr.write("\nINPUTTED CDS AND CONCATENATED EXONS DO NOT MATCH\n")

    protein_dict = get_codon_index(dna, cds)
    amino_position = []
    amino_acids = []

    for key, value in protein_dict.items():
        if value["Amino Acid"] == mid_aa:
            amino_position.append(key)
        amino_acids.append(value["Amino Acid"])

    if len(amino_position) == 0 and len(aa) == 1:
        raise AttributeError("The amino acid " + aa +" is not present in this sequence.")
    elif len(amino_position) == 0 and len(aa) > 1:
        raise AttributeError("The motif " + aa + " is not present in this sequence.")

    sys.stderr.write("\nSearching for gRNAs in loci...\n")

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

    guides = pd.DataFrame(amino_position, columns= ['Amino Acid Position'])
    guides["Context"] = guides.apply(lambda row: context(row["Amino Acid Position"], amino_acids, aa), axis = 1)

    if len(aa) > 1:
        indexNames = guides[guides['Context'] != aa].index
        guides.drop(indexNames , inplace=True)

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

    try:
        guides["5' Notes"] = guides.apply(lambda row: notes(row["5' gRNA Sequence"], row["5' gRNA G/C Content (%)"]), axis = 1)
        guides["3' Notes"] = guides.apply(lambda row: notes(row["3' gRNA Sequence"], row["3' gRNA G/C Content (%)"]), axis = 1)
    except:
        raise AttributeError("No " + aa + " motifs identified in sequence")

    sys.stderr.write("\nCounting off target matches...\n")

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
    if len(aa) == 1:
        guides.columns = ["Amino Acid Position", "Adjacent amino acids", '                        ', "gRNA Sequence 5' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count", '                        ', "gRNA Sequence 3' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count"] #Rename guide dataframe columns for output
    else:
        guides.columns = ["Amino Acid Position", "Targeted motif", '                        ', "gRNA Sequence 5' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count", '                        ', "gRNA Sequence 3' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count"] #Rename guide dataframe columns for output
    return guides

def get_options():

    import argparse

    parser = argparse.ArgumentParser(description='Target all amino acids of a specific type',
                                     prog='CRISPR-TAPE_specific')

    # input options
    iGroup = parser.add_argument_group('Input files')
    iGroup.add_argument('--loci', required=True, help='File containing genomic loci')
    iGroup.add_argument('--cds', required=True, help='File containing coding sequence')
    iGroup.add_argument('--genome', required=True, help='File containing genomic sequence')

    # target options
    tGroup = parser.add_argument_group('Targeting options')
    tGroup.add_argument('--aa', required=True, type = str, help='Amino acid short letter code')
    tGroup.add_argument('--motif', choices=['NGG', 'YG', 'TTTN'], default='NGG', type=str, help='Cas9 motif')

    # output options
    oGroup = parser.add_argument_group('Output options')
    oGroup.add_argument('--output', required=True, help='Directory and prefix for output guides')

    # other options
    otGroup = parser.add_argument_group('Other options')
    otGroup.add_argument('--up', default=' ', type=str, help='Additional bases upstream of loci')
    otGroup.add_argument('--down', default=' ', type=str, help='Additional bases downstream of loci')

    # combine
    args = parser.parse_args()

    # remove trailing forward slashes
    for arg in [args.loci, args.cds, args.genome]:
        arg = arg.replace(' ', '')

    return args

def main():

    args = get_options()
    # make output directory if absent
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # import input files
    with open(args.loci, 'r') as l:
        dna = l.read()
    with open(args.cds, 'r') as c:
        cds = c.read()
    with open(args.genome, 'r') as g:
        genome = g.read()

    # run function
    guides = general_function(args.aa,
                            args.motif,
                            dna,
                            cds,
                            args.up,
                            args.down,
                            genome)

    # save output
    guides.to_csv(os.path.join(args.output, args.output + '.csv'))
    sys.stderr.write("\nDone\n")

if __name__ == '__main__':
    main()

    sys.exit(0)
