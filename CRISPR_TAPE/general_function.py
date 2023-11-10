#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identify gRNAs immediately adjacent to all amino acids of a specific type
"""

import numpy as np
import pandas as pd
import os
import sys

from CRISPR_TAPE.shared_functions import clean_inputs, get_codon_index, PAMposition, analyse_text, notes, \
        pamcolumn, remove_pam, correct_distance, list_search, get_count, reverse_complement, find_sign_change

def context(start_index, char_list, substring):
    # Ensure the substring consists of characters from char_list
    assert all(char in char_list for char in substring if not char == "*"), "Invalid substring"
    # Convert char_list to a string for easier manipulation
    char_str = ''.join(char_list)
    # Determine the end index of the motif
    end_index = start_index + len(substring)
    # Check that the motif index is correct
    assert check_if_strings_equal(char_str[start_index: end_index], substring), "Index of amino acid motif is incorrect"
    #assert char_str[start_index: end_index] == substring, "Index of amino acid motif is incorrect"
    # Construct the context string
    left_index = start_index - 2
    if left_index < 0:
        left_index = 0
    right_index = end_index + 2
    if right_index > len(char_list):
        right_index = len(char_list)
    # get the context string
    context_str = char_str[left_index: start_index] + '(' + char_str[start_index: end_index] + ')' + char_str[end_index: right_index]
    return context_str

def process_upstream_distances(distances, gRNA, motif_length):
    # Initialise the list where we are storing the selected gRNA information
    sequences = []
    selected_positions = []
    for amino_acid_row in distances:
        # get the index of the distance immediately before a sign change
        index = find_sign_change(amino_acid_row)
        pos = None
        if index is not None and index + 1 < len(amino_acid_row):
            pos = index if amino_acid_row[index + 1] < -1 * (motif_length - 1) else min(index, index + 1, key=lambda x: abs(amino_acid_row[x]))
        elif amino_acid_row[0] >= 0:
            pos = len(gRNA) - 1
        elif (-1 * motif_length) < amino_acid_row[0] <= 0:
            pos = 0
        if pos is not None:
            sequence_info = {
                "Position": pos,
                "Sequence": gRNA['full gRNA Sequence'][pos],
                "Distance": amino_acid_row[pos],
                "Strand": gRNA['Strand'][pos],
                "G/C Content (%)": gRNA["G/C Content (%)"][pos]
            }
        else:
            sequence_info = {
                "Position": "",
                "Sequence": "No 5' guide could be identified",
                "Distance": "",
                "Strand": "",
                "G/C Content (%)": 0
            }
        # Keep track of the gRNAs we have chosen so that we do not choose the same gRNA again
        selected_positions.append(pos if pos is not None else "")
        sequences.append(sequence_info)
    return sequences, selected_positions

def process_downstream_distances(distances, gRNA, selected_positions, motif_length):
    sequences = []
    for i in range(len(distances)):
        # make a copy of the gRNA dataframe
        row_gRNA = gRNA.copy()
        # make a copy of the distances in the row
        row_distances = distances[i][:]
        # if we have chosen a 5' gRNA, then remove that from the gRNA options for this amino acid
        if not selected_positions[i] == "":
            row_gRNA = row_gRNA.drop([selected_positions[i]]).reset_index(drop=True)
            del row_distances[selected_positions[i]]
        # get the index of the distance immediately after a sign change
        index = find_sign_change(row_distances)
        pos = None
        if index is not None:
            index = index + 1
            pos = index if row_distances[index] < -2 else min([index - 1, index], key=lambda x: abs(row_distances[x]))
        elif (-1 * motif_length) < row_distances[0] < 0:
            pos = len(row_gRNA) - 1
        elif row_distances[0] >= 0:
            pos = 0
        if pos is not None:
            sequence_info = {
                "Position": pos,
                "Sequence": row_gRNA['full gRNA Sequence'][pos],
                "Distance": row_distances[pos],
                "Strand": row_gRNA['Strand'][pos],
                "G/C Content (%)": row_gRNA["G/C Content (%)"][pos]
            }
        else:
            sequence_info = {
                "Position": "",
                "Sequence": "No 3' guide could be identified",
                "Distance": "",
                "Strand": "",
                "G/C Content (%)": 0
            }
        sequences.append(sequence_info)
    return sequences

def check_if_strings_equal(query_motif, reference_motif):
    """ Returns true or false if the reference motifs are the same, allowing for * to denote any amino acid in a motif """
    assert len(query_motif) == len(reference_motif), "Motif lengths are different"
    motifs_identical = True
    for i in range(len(reference_motif)):
        if not reference_motif[i] == "*":
            if not reference_motif[i] == query_motif[i]:
                motifs_identical = False
    return motifs_identical


def find_motifs_in_amino_acid_sequence(motif, protein_dict):
    # iterate through the amino acids to get the indices of the amino acid motifs we are looking for
    amino_position = []
    amino_acids = []
    motif_length = len(motif)
    for amino_acid_index in protein_dict:
        if protein_dict[amino_acid_index]["Amino Acid"] == motif[0]:
            # determine if we have found the motif of interest
            indices_of_interest = [i for i in range(amino_acid_index, amino_acid_index + motif_length)]
            assert len(indices_of_interest) == motif_length
            if check_if_strings_equal("".join([protein_dict[i]["Amino Acid"] for i in indices_of_interest]), motif):
                amino_position.append(amino_acid_index)
        amino_acids.append(protein_dict[amino_acid_index]["Amino Acid"])
    return amino_position, amino_acids

def make_distance_matrix(first_list, second_list):
    # convert the positions to an array
    first_array = np.array(first_list)
    second_array = np.array(second_list)
    # get the distance of each cut site from each start and stop base
    return np.subtract(second_array, first_array[:,None])

def general_function(motif,
                    PAM,
                    dna,
                    cds,
                    hundredup,
                    hundreddown,
                    reference_genome):
    # clean the inputs
    dna_exon, dna, dna_edited, cds, reference_genome = clean_inputs(dna, cds, reference_genome, hundredup, hundreddown)
    # confirm that the inputted CDS sequence matches the CDS identified by the tool. If not the programme stops running.
    if dna_exon == cds:
        sys.stderr.write("\nInputted CDS and concatenated exons match\n")
    else:
        raise AttributeError("\nINPUTTED CDS AND CONCATENATED EXONS DO NOT MATCH\n")
    # convert the coding sequences to an amino acid dictionary
    protein_dict = get_codon_index(dna, cds)
    amino_position = set()
    amino_acids = []
    # iterate through the amino acids to get the indices of the amino acid motifs we are looking for
    amino_position, amino_acids = find_motifs_in_amino_acid_sequence(motif, protein_dict)
    # exit if the motif is not in the sequence
    motif_length = len(motif)
    assert not (len(amino_position) == 0 and motif_length == 1), "The amino acid " + motif +" is not present in this sequence."
    assert not (len(amino_position) == 0 and motif_length > 1), "The motif " + motif + " is not present in this sequence."
    sys.stderr.write("\nSearching for gRNAs in loci...\n")
    # find all guide RNAs in the genomic loci
    guide_positions, gRNA_list, guide_strands = PAMposition(dna_edited, PAM)
    # build the gRNA dataframe
    gRNA = pd.DataFrame({
        'full gRNA Sequence': gRNA_list,
        'Position': guide_positions,
        'Strand': guide_strands,
        'G/C Content (%)': [analyse_text(seq) for seq in gRNA_list]
    }).query("`G/C Content (%)` <= 75").sort_values(by='Position').reset_index(drop=True)
    # get the start and stop bases of all amino acids of interest
    start_base = []
    stop_base = []
    motif_aa_length = len(motif)
    for amino_acid_index in protein_dict:
        if amino_acid_index in amino_position:
            start_base.append(protein_dict[amino_acid_index]["base_1"][0])
            if not motif_aa_length > 1:
                stop_base.append(protein_dict[amino_acid_index]["base_3"][0])
            else:
                stop_base.append(protein_dict[amino_acid_index + motif_aa_length - 1]["base_3"][0])
    # make a matrix of the distance of all gRNA cut sites from all start bases.
    # -1 in this matrix means that the cut is to the right of the start base, relative to the strand of the inputted genomic loci
    start_distances = -1 * (make_distance_matrix(start_base, gRNA['Position']) + 1)
    # make a matrix of the distance of all gRNA cut sites from all stop bases
    # -1 in this matrix means that the cut is to the left of the start base, relative to the strand of the inputted genomic loci
    stop_distances = make_distance_matrix(stop_base, gRNA['Position'])
    # get the length of the motif in base pairs
    motif_bp_length = len(motif) * 3
    # choose the closest upstream gRNA for each amino acid of interest
    upstream_sequences, selected_positions = process_upstream_distances(start_distances.tolist(), gRNA, motif_bp_length)
    # choose the closest downstream gRNA for each amino acid of interest
    downstream_sequences = process_downstream_distances(stop_distances.tolist(), gRNA, selected_positions, motif_bp_length)
    # initialise the guide dataframe
    guides = pd.DataFrame(amino_position, columns=['Amino Acid Position'])
    # get the context of each amino acid
    guides["Context"] = guides.apply(lambda row: context(row["Amino Acid Position"], amino_acids, motif), axis = 1)
    # extract the 5' guide information
    guides["5' gRNA Sequence"] = [row["Sequence"] for row in upstream_sequences]
    guides["5' gRNA Sequence RC"] = guides["5' gRNA Sequence"].apply(reverse_complement)
    guides["Distance of 5' Cut Site from Amino Acid (bp)"] = [row["Distance"] for row in upstream_sequences]
    guides["5' gRNA Strand"] = [row["Strand"] for row in upstream_sequences]
    guides["5' gRNA G/C Content (%)"] = [row["G/C Content (%)"] for row in upstream_sequences]
    # extract the 3' guide information
    guides["3' gRNA Sequence"] = [row["Sequence"] for row in downstream_sequences]
    guides["3' gRNA Sequence RC"] = guides["3' gRNA Sequence"].apply(reverse_complement)
    guides["Distance of 3' Cut Site from Amino Acid (bp)"] = [row["Distance"] for row in downstream_sequences]
    guides["3' gRNA Strand"] = [row["Strand"] for row in downstream_sequences]
    guides["3' gRNA G/C Content (%)"] = [row["G/C Content (%)"] for row in downstream_sequences]
    assert len(guides) != 0, f"No {motif} motifs identified in sequence"
    # check for issues with the guide RNAs
    guides["5' Notes"] = guides.apply(lambda row: notes(row["5' gRNA Sequence"], row["5' gRNA G/C Content (%)"]), axis = 1)
    guides["3' Notes"] = guides.apply(lambda row: notes(row["3' gRNA Sequence"], row["3' gRNA G/C Content (%)"]), axis = 1)
    # count the number of exact off target matches
    # (NOTE: this can be made to work with non-exact matches by mapping the guide RNAs back to the reference using e.g. Minimap2)
    sys.stderr.write("\nCounting off target exact matches...\n")
    ls = list(guides["5' gRNA Sequence"]) + list(guides["5' gRNA Sequence RC"]) + list(guides["3' gRNA Sequence"]) + list(guides["3' gRNA Sequence RC"])
    counter = list_search(ls, reference_genome)
    guides["5' gRNA Off Target Count"] = guides.apply(lambda row: get_count(row["5' gRNA Sequence"], row["5' gRNA Sequence RC"], counter), axis=1)
    guides["3' gRNA Off Target Count"] = guides.apply(lambda row: get_count(row["3' gRNA Sequence"], row["3' gRNA Sequence RC"], counter), axis=1)
    # make a separate column for the PAM sequence
    guides["5' PAM"] = guides.apply(lambda row: pamcolumn(row["5' gRNA Sequence"], row["5' gRNA Strand"], PAM), axis =1)
    guides["3' PAM"] = guides.apply(lambda row: pamcolumn(row["3' gRNA Sequence"], row["3' gRNA Strand"], PAM), axis =1)
    # remove the PAM from the gRNA sequence
    guides["5' gRNA Sequence"] = guides.apply(lambda row: remove_pam(row["5' gRNA Sequence"], row["5' gRNA Strand"], PAM), axis =1)
    guides["3' gRNA Sequence"] = guides.apply(lambda row: remove_pam(row["3' gRNA Sequence"], row["3' gRNA Strand"], PAM), axis =1)
    # make the amino acid positions 1-based
    guides["Amino Acid Position"] = [p + 1 for p in list(guides["Amino Acid Position"])]
    # sort the rows by ascending amino acid position
    guides = guides.sort_values(by=["Amino Acid Position"], ascending=True).reset_index(drop=True)
    # reorganise guide dataframe columns for output
    guides = guides[["Amino Acid Position",
                    "Context",
                    "5' gRNA Sequence",
                    "5' PAM",
                    "5' gRNA Strand",
                    "5' gRNA G/C Content (%)",
                    "Distance of 5' Cut Site from Amino Acid (bp)",
                    "5' Notes",
                    "5' gRNA Off Target Count",
                    "3' gRNA Sequence",
                    "3' PAM",
                    "3' gRNA Strand",
                    "3' gRNA G/C Content (%)",
                    "Distance of 3' Cut Site from Amino Acid (bp)",
                    "3' Notes",
                    "3' gRNA Off Target Count"]]
    return guides

def get_options():

    import argparse

    parser = argparse.ArgumentParser(description='Target all amino acids or multi-amino acid motifs in a sequence',
                                    prog='CRISPR-TAPE_general')

    # input options
    iGroup = parser.add_argument_group('Input files')
    iGroup.add_argument('--loci', required=True, help='File containing genomic loci')
    iGroup.add_argument('--cds', required=True, help='File containing coding sequence')
    iGroup.add_argument('--reference-genome', required=True, help='File containing genomic sequence')

    # target options
    tGroup = parser.add_argument_group('Targeting options')
    tGroup.add_argument('--motif', required=True, type = str, help='Amino acid short letter code or a run of adjacent amino acids to target. "*" can be used to denote an ambiguous base position.')
    tGroup.add_argument('--PAM', choices=['NGG'], default='NGG', type=str, help='Cas9 motif')

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
    with open(args.reference_genome, 'r') as g:
        genome = g.read()

    # run function
    guides = general_function(args.motif,
                            args.PAM,
                            dna,
                            cds,
                            args.up,
                            args.down,
                            genome)

    # save output
    if not ".csv" in args.output:
        args.output = args.output + ".csv"
    guides.to_csv(os.path.join(args.output, args.output))
    sys.stderr.write("\nDone\n")

if __name__ == '__main__':
    main()

    sys.exit(0)
