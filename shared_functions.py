#!/usr/bin/env python3
"""
Functions shared by specific and general functions
"""
import re 
import pandas as pd
from tqdm import tqdm

def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3] #Defining a codon as 3 bases
            protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
    return protein

def clean_inputs(loci, coding_sequence, organism_genome, hundredup, hundreddown):
    """ clean and reformat user inputs """
    
    remove_lower = lambda text: re.sub('[a-z]', '', text) #Remove any characters that aren't a letter
    loci_exon = remove_lower(loci).replace(" ", "").replace("\n", "").replace("\r", "")  #Remove any lower case characters. These correspond to intron sequences so the output is the exon sequence of the gene.
    loci = hundredup.lower() + loci + hundreddown.lower() #Append the 100 bases up and downstream to the genomic loci
    loci = loci.replace(" ", "").replace("\n", "").replace("\r", "") #Remove spaces
    loci_edited = loci.upper()
    organism_genome = re.sub(r'[^ACTG]', '', organism_genome) #Remove any characters that are not ACTG
    coding_sequence = coding_sequence.replace("\n", "").replace("\r", "").replace(" ", "").upper()
    
    return loci_exon, loci, loci_edited, coding_sequence, organism_genome

def PAMposition(string, motif):
    """ get index position of all PAMs within the genomic loci """
    
    position = [] #Empty list to store identified gRNAs
    entry = []
    strand = []
    
    if motif == "NGG":
        for n in tqdm(range(len(string) - 1)):
            if string[n] == 'G' and string[n+1] == 'G' and n-21 >= 0: #If two Gs in a row and Gs not at the start of the genomic loci
                position.append(n-5) #Append cut site
                entry.append(string[n-21:n+2].upper())
                strand.append("forward")
            if string[n-1] == "C" and string[n] == "C" and n+21 <= len(string):
                position.append(n+5) #Append cut site
                entry.append(string[n-1:n+22].upper())
                strand.append("reverse")
                
    elif motif == "YG":
        for n in tqdm(range(len(string) - 1)):
            if string[n] == 'C' or string[n] == 'T' and string[n + 1] == 'G' and n-21 >= 0: #If a C or T is followed by G
                position.append(n-9) #Append cut site
                entry.append(string[n-21:n+2].upper())
                strand.append("forward")
            if string[n-1] == "C" and string[n+1] == "A" or string[n+1] == "G" and n+21 <= len(string):
                position.append(n+9) #Append cut site
                entry.append(string[n-1:n+22].upper())
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
                    position.append(n-23) #Append cut site
                    entry.append(string[n-30:n+1])
                    strand.append("reverse")
    
    else:
    raise ValueError("PAM not recognised")
    return position, entry, strand

def get_codon_index(dna, cds):
    
    Base_df = pd.DataFrame(list(dna), columns = ["Base"]) #List of all the bases in the genomic loci
    Base_df["Position"] = Base_df.index
    Base_df = Base_df[Base_df['Base'].str.istitle()].reset_index(drop=True) #Remove lowercase bases so dataframe only codes for exon
    
    pos1 = Base_df[Base_df.index % 3 == 0].reset_index(drop=True) #Generate a list of every third base from the first base
    pos2 = Base_df[Base_df.index % 3 == 1].reset_index(drop=True) #Generate a list of every third base from the second base
    pos3 = Base_df[Base_df.index % 3 == 2].reset_index(drop=True) #Generate a list of every third base from the third base
    aas = list(translate(cds))
    
    protein_dict = {}
    for aa in range(len(aas)):
        protein_dict.update({aa: {"Amino Acid": aas[aa], "base_1": (pos1["Position"][aa], pos1["Base"][aa]), "base_2": (pos2["Position"][aa], pos2["Base"][aa]), "base_3": (pos3["Position"][aa], pos3["Base"][aa])}})

    return protein_dict

def analyse_text(text): #Analyse the G/C content of the guide RNA as a percentage.
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
 
def notes(string, content): #Return key information on the generated guide RNA
    polyt = ''
    for n in range(len(string)-3):
        if string[n:n + 4] == 'TTTT': #Check for four thymines in a row
            polyt = 'PolyT present. '
    if string[0] != 'G': #Check the guide RNA starts with a 'G' at the most 5' position
        polyt += 'No leading G. '
    if content >= 75:
        polyt += 'G/C content over 75%. ' #Check if the G/C content of the guide is more than or equal to 75%
    return polyt

def pamcolumn(entry, strand, reverse_entry, motif): #Add a column to the guide dataframe specifying the pam adjacent to the guide RNA generated.
    if not strand == "":
        if strand == "forward":
            guide = entry
        else:
            guide = reverse_entry
        if motif == 'NGG' or motif == 'YG': #PAMs at 3' of the guide RNA
                pam = guide[20:]
        if motif == 'TTTN': #PAM at the 5' of the guide RNA
                  pam = guide[0:4]
    else:
        pam = ""
    return pam

def remove_pam(entry, strand, reverse_entry, motif): #Remove the PAM from the guide RNA column
    if not strand == "":
        if strand == "forward":
            guide = entry
        else:
            guide = reverse_entry
        if motif == 'NGG' or motif == 'YG':
                gRNA = guide[0:20]
        if motif == 'TTTN':
                gRNA = guide[4:]
    else:
        gRNA = entry
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
