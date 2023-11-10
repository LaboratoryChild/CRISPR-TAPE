#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identify guides surrounding a specific residue position within a maximum distance range
"""
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import os

from CRISPR_TAPE.shared_functions import clean_inputs, get_codon_index, PAMposition, analyse_text, notes, pamcolumn, remove_pam, correct_distance, list_search, get_count, reverse_complement

def off_target(search, reverse_search, reference_genome):
    if not search == "":
        term = [search,reverse_search]
        count = sum(reference_genome.count(i) for i in term)- 1
    else:
        count = ""
    return count

def specific_function(spec_amino,
                    distance,
                    motif,
                    dna,
                    cds,
                    hundredup,
                    hundreddown,
                    reference_genome):

    dna_exon, dna, dna_edited, cds, reference_genome = clean_inputs(dna, cds, reference_genome, hundredup, hundreddown)

    if dna_exon == cds: #Confirm that the inputted CDS sequence matches the CDS identified by the tool. If not the programme stops running.
        sys.stderr.write("Inputted CDS and concatenated exons match")
    else:
        sys.stderr.write("INPUTTED CDS AND CONCATENATED EXONS DO NOT MATCH")

    protein_dict = get_codon_index(dna, cds)
    for key, value in protein_dict[spec_amino - 1].items():
        if key == "base_2":
            middle = value[0]
        if key == "Amino Acid":
            selectedaa = value

    selectionmade = selectedaa + '-' + str(spec_amino)

    sys.stderr.write('\nThe amino acid you have selected is ' + selectionmade + '\n') #Confirm the selection is the correct amino acid
    sys.stderr.write('\nThe motif you have selected is ' + motif + '\n')
    # generate lists of all the potential guides within the genomic loci and a list of all guides on the reverse complement.
    sys.stderr.write("\nSearching for gRNAs in loci...\n")
    guide_positons, gRNA_list, guide_strands = PAMposition(dna_edited, motif)#Identify PAMs in the inputted genomic loci
    # convert the lists into dataframes
    gRNA = pd.DataFrame(gRNA_list, columns= ['full gRNA Sequence']) #Convert the gRNA list into a dataframe
    gRNA['Position'] = guide_positons #Add a column of guide RNA cut site positions to the dataframe
    gRNA['Strand'] = guide_strands #Add a column to specify these guides are on the forward strand relative to the inputted genomic loci.
    if "" in gRNA_list:
        gRNA = gRNA.drop(gRNA[gRNA['full gRNA Sequence'] == ""].index)
    gRNA = gRNA.sort_values(by=['Position']).reset_index(drop=True) #Sort the guide RNAs by their position in the genomic loci
    gRNA["Distance from Amino Acid (bp)"] = np.array(gRNA["Position"]) - middle #Determine guide distance from base at the centre of the codon

    upperguides = "No guides within distance range specified"
    downerguides = "No guides within distance range specified"

    if not len(gRNA) == 0:
        sys.stderr.write('\n' + str(len(gRNA)) + " guides found\n")
        slice_index = " "
        for dis in range(len(gRNA)-1):
            if gRNA["Distance from Amino Acid (bp)"][dis] < 0 and gRNA["Distance from Amino Acid (bp)"][dis+1] > 0:
                slice_index = dis
        if not slice_index == " ":
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
        if len(upperguides) == 0:
            downerguides = noguides
        else:
            upperguides = upperguides.sort_values(by = ['Distance from Amino Acid (bp)']) #Arrange guide RNAs by their distance from the amino acid
            upperguides['Distance from Amino Acid (bp)'] = upperguides.apply(lambda row: correct_distance(int(row['Distance from Amino Acid (bp)']), row['Strand']), axis =1) + 1
            upperguides["Reverse complement"] = upperguides.apply(lambda row: reverse_complement(row['full gRNA Sequence']), axis =1)
            upperguides["PAM"] = upperguides.apply(lambda row: pamcolumn(row['full gRNA Sequence'], row['Strand'], motif), axis =1)
            upperguides["gRNA Sequence"] = upperguides.apply(lambda row: remove_pam(row['full gRNA Sequence'], row['Strand'], motif), axis =1)
            upperguides['G/C Content (%)'] = upperguides.apply(lambda row: analyse_text(row['full gRNA Sequence']), axis =1) #Calculate GC percentage
            upperguides["Notes"] = upperguides.apply(lambda row: notes(row["full gRNA Sequence"], row["G/C Content (%)"]), axis=1) #Output notes of key guide RNA characterstics to a new column
            upperguides = upperguides.reset_index(drop=True)
    else:
        upperguides = noguides

    amino_acid = pd.DataFrame(columns=["full gRNA Sequence", "Position", "Strand", "Distance from Amino Acid (bp)", "Reverse complement", "PAM", "gRNA Sequence", "G/C Content (%)", "Notes", 'Off Target Count']) #Generate a new amino acid for the amino acid target information
    amino_acid = amino_acid.append({"full gRNA Sequence": "", "Position" : "", "Strand":"","Distance from Amino Acid (bp)":"","Reverse complement": "", "PAM": "", "gRNA Sequence": selectionmade, "G/C Content (%)":"","Notes":"", 'Off Target Count': ""}, ignore_index=True) #Generate a new amino acid for the amino acid target information

    if not isinstance(downerguides, str):
        downerguides = downerguides[downerguides["Distance from Amino Acid (bp)"] <= distance]
        if len(downerguides) == 0:
            downerguides = noguides
        else:
            downerguides = downerguides.sort_values(by=['Distance from Amino Acid (bp)']) #Arrange guide RNAs by their distance from the amino acid
            downerguides['Distance from Amino Acid (bp)'] = downerguides.apply(lambda row: correct_distance(int(row['Distance from Amino Acid (bp)']),row['Strand']), axis =1) - 1
            downerguides["Reverse complement"] = downerguides.apply(lambda row: reverse_complement(row['full gRNA Sequence']), axis =1)
            downerguides["PAM"] = downerguides.apply(lambda row: pamcolumn(row['full gRNA Sequence'], row["Strand"], motif), axis =1)
            downerguides["gRNA Sequence"] = downerguides.apply(lambda row: remove_pam(row['full gRNA Sequence'], row["Strand"], motif), axis =1)
            downerguides['G/C Content (%)'] = downerguides.apply(lambda row: analyse_text(row['full gRNA Sequence']), axis =1) #Calculate GC percentage
            downerguides["Notes"] = downerguides.apply(lambda row: notes(row["full gRNA Sequence"], row["G/C Content (%)"]), axis=1) #Output notes of key guide RNA characterstics to a new column
            downerguides = downerguides.reset_index(drop=True)
    else:
        downerguides = noguides

    guides = pd.concat([upperguides, amino_acid, downerguides])
    sys.stderr.write("\nCounting off target matches...\n")
    ls = list(guides['full gRNA Sequence'] + guides["Reverse complement"])
    counter = list_search(ls, reference_genome)
    guides['Off Target Count'] = guides.progress_apply(lambda row: get_count(row["full gRNA Sequence"], row["Reverse complement"], counter), axis=1)
    guides = guides[['Distance from Amino Acid (bp)', 'gRNA Sequence', 'PAM', 'Strand', 'G/C Content (%)', 'Off Target Count', 'Notes']] #Reorganise the guide RNA dataframe

    return guides

def get_options():

    import argparse

    parser = argparse.ArgumentParser(description='Target all amino acids of a specific type',
                                    prog='CRISPR-TAPE_specific')

    # input options
    iGroup = parser.add_argument_group('Input files')
    iGroup.add_argument('--loci', required=True, help='File containing genomic loci')
    iGroup.add_argument('--cds', required=True, help='File containing coding sequence')
    iGroup.add_argument('--reference-genome', required=True, help='File containing genomic sequence')

    # target options
    tGroup = parser.add_argument_group('Targeting options')
    tGroup.add_argument('--spec-amino', required=True, type = int, help='Amino acid/residue position')
    tGroup.add_argument('--PAM', choices=['NGG'], default='NGG', type=str, help='Cas9 motif')
    tGroup.add_argument('--distance', default=10000, type=int, help='Maximum distance from target (base pairs)')

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
    for arg in [args.loci, args.cds, args.reference_genome]:
        arg = arg.rstrip('\\')

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
    with open(args.reference_genome, 'r') as g:
        reference_genome = g.read()

    # run function
    guides = specific_function(args.spec_amino,
                            args.distance,
                            args.PAM,
                            dna,
                            cds,
                            args.up,
                            args.down,
                            reference_genome)

    # save output
    guides.to_csv(os.path.join(args.output, args.output + '.csv'))
    sys.stderr.write("\nDone\n")

if __name__ == '__main__':
    main()

    sys.exit(0)
