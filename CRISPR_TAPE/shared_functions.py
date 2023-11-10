#!/usr/bin/env python3
"""
Functions shared by specific and general functions
"""
import re
import pandas as pd
from tqdm import tqdm

def translate(seq):
    """ Translate the exon sequence of the gene into its respective amino acid codes using a dictionary """
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
        'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }
    # initialise the protein sequence
    protein = ""
    # if the length of the coding sequence is divisible by 3
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            # Defining a codon as 3 bases
            codon = seq[i:i + 3]
            # Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
            protein += table[codon]
    return protein

def clean_inputs(loci,
            coding_sequence,
            reference_genome,
            hundred_up,
            hundred_down):
    """ clean and reformat user inputs """
    # Remove any characters that aren't a letter
    remove_lower = lambda text: re.sub('[a-z]', '', text)
    # Remove any lower case characters. These correspond to intron sequences so the output is the exon sequence of the gene.
    loci_exon = remove_lower(loci).replace(" ", "").replace("\n", "").replace("\r", "")
    # Append the 100 bases up and downstream to the genomic loci
    loci = hundred_up.lower() + loci + hundred_down.lower()
    # Remove spaces
    loci = loci.replace(" ", "").replace("\n", "").replace("\r", "")
    # Uppercase the loci
    loci_edited = loci.upper()
    # Remove any characters that are not ACTG
    reference_genome = re.sub(r'[^ACTG]', '', reference_genome.upper())
    coding_sequence = coding_sequence.replace("\n", "").replace("\r", "").replace(" ", "").upper()
    return loci_exon, loci, loci_edited, coding_sequence, reference_genome

def PAMposition(string, motif):
    """ get index position of all PAMs within the genomic loci """
    # initialise lists to store gRNA information
    position = []
    entry = []
    strand = []
    # check what PAM we are looking for
    if motif == "NGG":
        for n in tqdm(range(len(string) - 1)):
            if string[n] == 'G' and string[n+1] == 'G' and n-21 >= 0: #If two Gs in a row and Gs not at the start of the genomic loci
                position.append(n-5) #Append cut site
                entry.append(string[n-21:n+2].upper())
                strand.append("forward")
            if string[n-1] == "C" and string[n] == "C" and n+21 <= len(string):
                position.append(n+4) #Append cut site
                entry.append(reverse_complement(string[n-1:n+22].upper()))
                strand.append("reverse")
    elif motif == "YG":
        for n in tqdm(range(len(string) - 1)):
            if (string[n] == 'C' or string[n] == 'T') and string[n + 1] == 'G' and n-21 >= 0: #If a C or T is followed by G
                position.append(n-9) #Append cut site
                entry.append(string[n-21:n+2].upper())
                strand.append("forward")
            if string[n-1] == "C" and (string[n+1] == "A" or string[n+1] == "G") and n+21 <= len(string):
                position.append(n+8) #Append cut site
                entry.append(reverse_complement(string[n-1:n+22].upper()))
                strand.append("reverse")
    elif motif == "TTTN":
        for n in tqdm(range(len(string) - 7)):
            if string[n] == 'T' and string[n+1] == 'T' and string[n+2] == 'T': #If 3 Ts in a row
                if string[n+3] == 'A' or string[n+3] == 'C' or string[n+3] == 'G': #If 4th base is not T
                    position.append(n+23) #Append cut site
                    entry.append(string[n:n+31])
                    strand.append("forward")
            if string[n-3] == 'T' or string[n-3] == 'G' or string[n-3] == 'C':
                if string[n-2] == 'A' and string[n-1] == 'A' and string[n] == 'A':
                    position.append(n-24) #Append cut site
                    entry.append(reverse_complement(string[n-30:n+1]))
                    strand.append("reverse")
    else:
        raise ValueError("PAM not recognised")
    return position, entry, strand

def get_codon_index(dna, cds):
    """ return a dictionary of each amino acid position and the positions of the bases of the codon in the genomic loci """
    assert len(cds) % 3 == 0, "The length of the coding sequence is not a multiple of 3"
    Base_df = pd.DataFrame(list(dna), columns = ["Base"]) #List of all the bases in the genomic loci
    Base_df["Position"] = Base_df.index
    Base_df = Base_df[Base_df['Base'].str.istitle()].reset_index(drop=True) #Remove lowercase bases so dataframe only codes for exon

    pos1 = Base_df[Base_df.index % 3 == 0].reset_index(drop=True) #Generate a list of every third base from the first base
    pos2 = Base_df[Base_df.index % 3 == 1].reset_index(drop=True) #Generate a list of every third base from the second base
    pos3 = Base_df[Base_df.index % 3 == 2].reset_index(drop=True) #Generate a list of every third base from the third base
    aas = list(translate(cds))
    protein_dict = {}
    for aa in range(len(aas)):
        protein_dict[aa] = {"Amino Acid": aas[aa],
                        "base_1": (pos1["Position"][aa],
                                pos1["Base"][aa]),
                        "base_2": (pos2["Position"][aa],
                                pos2["Base"][aa]),
                        "base_3": (pos3["Position"][aa],
                                pos3["Base"][aa])}

    return protein_dict

def analyse_text(text):
    """ Return the G/C content of the guide RNA as a percentage."""
    count = 0
    letter_count = 0
    for char in text:
        if char.isalpha(): #Count the length of the guide RNA
            count += 1
        if char == "C" or char =="G":
            letter_count += 1 #Count the number of Gs or Cs in the guide RNA
    perc = float(letter_count)/float(count) * 100 #Calculate a percentage of Gs and Cs in the length of the guide RNA
    perc = round(perc, 2) #Round to two decimal places
    return perc

def notes(string, content):
    """ Return key information on the generated guide RNA """
    if string == "No 5' guide could be identified" or string == "No 3' guide could be identified":
        return ""
    note = []
    for n in range(len(string)-3):
        if string[n:n + 4] == 'TTTT': #Check for four thymines in a row
            note.append('PolyT present.')
    if string[0] != 'G': #Check the guide RNA starts with a 'G' at the most 5' position
        note.append('No leading G.')
    if content >= 75:
        note.append('G/C content over 75%.') #Check if the G/C content of the guide is more than or equal to 75%
    return " ".join(note)

def pamcolumn(guide, strand, motif):
    """ Return a string of the pam adjacent to the guide RNA generated."""
    if guide == "No 5' guide could be identified" or guide == "No 3' guide could be identified":
        return ""
    # check that the gRNA is the correct length
    if motif == "NGG" or motif == "YG":
        assert len(guide) == 23, f"gRNA: {guide} is the wrong length"
    # extract the PAM for this gRNA depending on the motif
    assert not strand == ""
    # PAMs are at the 3' end of the guide
    if motif == 'NGG' or motif == 'YG':
        pam = guide[20:]
    if motif == 'YG':
        pam = guide[21:]
    # PAMs are at the 5' end of the guide
    if motif == 'TTTN':
        pam = guide[0:4]
    assert not pam == "", f"No {motif} was found in the guide: {guide}"
    return pam

def remove_pam(guide, strand, PAM):
    """ Return a string of the PAM-less guide RNA """
    if strand == "":
        return guide
    PAM_removed = False
    if PAM == 'NGG':
        gRNA = guide[:20]
        removed_PAM = guide[20:]
        PAM_removed = True
        assert any(removed_PAM == p for p in ["AGG", "TGG", "GGG", "CGG"])
    if PAM == 'YG':
        gRNA = guide[0:21]
        removed_PAM = guide[21:]
        PAM_removed = True
        assert any(removed_PAM == p for p in ["CG", "TG"])
    if PAM == 'TTTN':
        gRNA = guide[4:]
        removed_PAM = guide[:4]
        PAM_removed = True
        assert any(removed_PAM == p for p in ["TTTT", "TTTA", "TTTG", "TTTC"])
    assert PAM_removed, f"{PAM} could not removed from {guide}"
    return gRNA

def correct_distance(distance, strand):
    if not distance == "":
        if distance > 0 and strand == "reverse":
            distance = distance - 1
        if distance <= 0 and strand == "forward":
            distance = distance + 1
    else:
        distance = distance
    return distance

def list_search(guides_of_interest,
            reference_genome):
    """ Return a dictionary counting the number of times each guide of interest is seen in the reference genome. """
    found_guides = {}
    guide_length = len(guides_of_interest[0])
    guides_of_interest = set(guides_of_interest)  # Convert list to set for O(1) lookups
    # Loop through the reference genome and check for presence of each guide
    for x in tqdm(range(len(reference_genome) - guide_length + 1)):
        current_substring = reference_genome[x:x + guide_length]
        if current_substring in guides_of_interest:
            found_guides[current_substring] = found_guides.get(current_substring, 0) + 1
    return found_guides

def get_count(fw_query, rv_query, counter):
    """ return an integer of the number of times - 1 that a string occurs in a counter"""
    count = -1
    for gRNA in counter:
        if fw_query == gRNA or rv_query == gRNA:
            count += int(counter[gRNA])
    return count

def reverse_complement(dna):
    """ Return a string of the reverse complement sequence of a given nucleotide sequence """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    comp = ''.join([complement[base] if base in complement else base for base in dna[::-1]])
    return comp

def find_sign_change(row):
    return next((i for i, val in enumerate(row[:-1])
            if val * row[i + 1] <= 0), None)