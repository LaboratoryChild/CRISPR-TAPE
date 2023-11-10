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

from CRISPR_TAPE.shared_functions import clean_inputs, get_codon_index, PAMposition, analyse_text, notes, pamcolumn, remove_pam, list_search, get_count, reverse_complement, find_sign_change

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
    # convert the CDS to a protein dictionary
    protein_dict = get_codon_index(dna, cds)
    # display information about the selected amino acid
    selectionmade = protein_dict[spec_amino - 1]["Amino Acid"] + '-' + str(spec_amino)
    sys.stderr.write('\nThe amino acid you have selected is ' + selectionmade + '\n') #Confirm the selection is the correct amino acid
    sys.stderr.write('The PAM you have selected is ' + motif + '\n')
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
    assert not len(gRNA) == 0, "No guide RNAs identified in sequence"
    sys.stderr.write('\n' + str(len(gRNA)) + " guides found\n")
    # get the distance of the gRNAs from the start and end base of the amino acid
    start_base = protein_dict[spec_amino - 1]["base_1"][0]
    end_base = protein_dict[spec_amino - 1]["base_3"][0]
    gRNA["Distance from start of Amino Acid (bp)"] = np.array(gRNA["Position"]) - start_base + 1
    gRNA["Distance from end of Amino Acid (bp)"] = np.array(gRNA["Position"]) - end_base
    # define a dataframe for if there are no guides within the distance range
    noguides = pd.DataFrame(columns=["full gRNA Sequence", "Position", "Strand", "Distance from Amino Acid (bp)", "Reverse complement", "PAM", "gRNA Sequence", "G/C Content (%)", "Notes", 'Off Target Count']) #Generate a new amino acid for the amino acid target information
    noguides.loc[0] = {"full gRNA Sequence": "", "Position" : "", "Strand":"","Distance from Amino Acid (bp)":"","Reverse complement": "", "PAM": "", "gRNA Sequence": "No guides within distance range specified", "G/C Content (%)":"","Notes":"", 'Off Target Count': ""} #Generate a new amino acid for the amino acid target information
    # get the upstream guide RNA dataframe
    upstream_guides = gRNA.query("`Distance from start of Amino Acid (bp)` < 3").query(f"`Distance from start of Amino Acid (bp)` >= {-distance}")
    upstream_guides = upstream_guides.rename(columns={"Distance from start of Amino Acid (bp)": "Distance from Amino Acid (bp)"})
    upstream_guides = upstream_guides.drop("Distance from end of Amino Acid (bp)", axis=1)
    if len(upstream_guides) == 0:
        upstream_guides = noguides
    # get the downstream guide RNA dataframe
    downstream_guides = gRNA.query("`Distance from end of Amino Acid (bp)` > -3").query(f"`Distance from end of Amino Acid (bp)` <= {distance}")
    downstream_guides = downstream_guides.rename(columns={"Distance from end of Amino Acid (bp)": "Distance from Amino Acid (bp)"})
    downstream_guides = downstream_guides.drop("Distance from start of Amino Acid (bp)", axis=1)
    if len(downstream_guides) == 0:
        downstream_guides = noguides
    # make a row for the amino acid
    amino_acid = pd.DataFrame(columns=["full gRNA Sequence", "Position", "Strand", "Distance from Amino Acid (bp)", "Reverse complement", "PAM", "gRNA Sequence", "G/C Content (%)", "Notes", 'Off Target Count']) #Generate a new amino acid for the amino acid target information
    amino_acid.loc[0] = {"full gRNA Sequence": "", "Position" : "", "Strand":"","Distance from Amino Acid (bp)":"","Reverse complement": "", "PAM": "", "gRNA Sequence": selectionmade, "G/C Content (%)":"","Notes":"", 'Off Target Count': ""}#Generate a new amino acid for the amino acid target information
    # join the dataframes together
    guides = pd.concat([upstream_guides, amino_acid, downstream_guides])
    # add a progress bar
    tqdm.pandas()
    # get the reverse complement of the guides
    guides["Reverse complement"] = guides.apply(lambda row: reverse_complement(row['full gRNA Sequence']) if row["full gRNA Sequence"] != "" else "", axis = 1)
    # get the PAM of the guides
    guides["PAM"] = guides.apply(lambda row: pamcolumn(row['full gRNA Sequence'], row['Strand'], motif) if row["full gRNA Sequence"] != "" else "", axis = 1)
    # get the guide without the PAM
    guides["gRNA Sequence"] = guides.apply(lambda row: remove_pam(row['full gRNA Sequence'], row['Strand'], motif) if row["full gRNA Sequence"] != "" else "", axis = 1)
    # get the G/C content
    guides['G/C Content (%)'] = guides.apply(lambda row: analyse_text(row['full gRNA Sequence']) if row["full gRNA Sequence"] != "" else "", axis = 1)
    # get the notes
    guides["Notes"] = guides.apply(lambda row: notes(row["full gRNA Sequence"], row["G/C Content (%)"]) if row["full gRNA Sequence"] != "" else "", axis=1) #Output notes of key guide RNA characterstics to a new column
    # count the off target matches
    sys.stderr.write("\nCounting off target matches...\n")
    ls = list(guides['full gRNA Sequence']) + list(guides["Reverse complement"])
    counter = list_search(ls, reference_genome)
    guides['Off Target Count'] = guides.progress_apply(lambda row: get_count(row["full gRNA Sequence"], row["Reverse complement"], counter) if row["full gRNA Sequence"] != "" else "", axis=1)
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
    if not ".csv" in args.output:
        args.output = args.output + '.csv'
    guides.to_csv(args.output)
    sys.stderr.write("\nDone\n")

if __name__ == '__main__':
    main()

    sys.exit(0)
