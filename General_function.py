#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 00:28:49 2019

@author: danielanderson
"""

import numpy as np
import pandas as pd
import re

def General_function(aa, motif, cds, dna, hundredup, hundreddown, orgen):
    
    def reverse_complement(dna):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([complement[base] for base in dna[::-1]])

    psn = []
    def baseposition(code): #Get the position of each base within the genomic loci and output to a list
        for  n in range(len(code)):
            i = [n + 1]
            psn.append(i)
        return psn

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
            'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
            'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
        }
        protein =""
        if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3] #Defining a codon as 3 bases
                protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
        return protein
        
    def charposition(string, char):
        pos = [] #list to store positions for each amino acid in the amino acid sequence
        for n in range(len(string)):
            if string[n] == char:
                pos.append(n + 1)
        return pos

    ################################ ADD FUNCTION IF ADDITTIONAL PAMS ARE TO BE INCLUDED ################################

    # get index position of all PAMs within the genomic loci
    def PAMposition(string): #Identify the position of NGGs within the genomic loci
        pos = [] #Empty list to store identified gRNAs
        for n in range(len(string) - 1):
            if string[n] == 'G' and string[n+1] == 'G' and n-21 >= 0: #If two Gs in a row and Gs not at the start of the genomic loci
                pos.append(n-4) #Append the position of the base 5' of the cut site
        return pos

    def YGposition(string): #Identify the position of YGs in the genomic loci (Y is a pyrimidine- C or T)
        pos = []
        for n in range(len(string) - 1):
            if string[n] == 'C' or string[n] == 'T': #If a C or T is followed by G
                if string[n + 1] == 'G':
                    pos.append(n-9)  #Append the position of the base 5' of the cut site
        return pos

    def TTTNposition(string): #Identify the position of TTTN motifs in the genomic loci
        pos = []
        for n in range(len(string) - 7):
            if string[n] == 'T' and string[n+1] == 'T' and string[n+2] == 'T': #If 3 Ts in a row
                if string[n+3] == 'A' or string[n+3] == 'C' or string[n+3] == 'G': #If 4th base is not T
                    pos.append(n + 23) #Append the position of the base 5' of the cut site
        return pos

    ###############################################################################################################

    #Determine closest downstream gRNA to specified aa
               
    def closest_downstream(posa): #Determine which guide RNA is closest downstream of the 3' base coding for the amino acid.
        small = 1000000000
        position = 0
        for x in gRNA.iloc[:, 1]:
            posg = x
            diff = int(posa) - int(posg)
            if diff > -2 and diff < small: #iterate through the list of guide RNAs and their positions to determine the cut site resistance that is the smallest distance away from the 3' base of the amino acid.
                small = diff
                position = posg
        return position #Returns the cut site position of the closest 3' guide  RNA.

    #Determine closest upstream gRNA to specified aa

    def closest_upstream(posa): #Determine which guide RNA is closest upstream of the 5' base coding for the amino acid.
        small = 1000000000
        position = 0
        for x in gRNA.iloc[:, 1]:
            posg = x
            diff = int(posg) - int(posa)
            if diff >-2 and diff < small: #iterate through the list of guide RNAs and their positions to determine the cut site resistance that is the smallest distance away from the 3' base of the amino acid.
                small = diff
                position = posg
        return position #Returns the cut site position of the closest 3' guide  RNA.

    def position_finder(x):
        mask = np.where(gRNA['Position'] == x) #Determine which guide RNA sequence corresponds to the position of the closest guide RNA
        indexes = mask[0]
        target = gRNA['gRNA Sequence'][indexes]
        try:
            target = target.iloc[0]
        except:
            target = 0
        return target

    def strand_finder(x): #Determine on which strand the closest guide RNA identified is located
        mask = np.where(gRNA['Position'] == x)
        indexes = mask[0]
        target = gRNA['Strand'][indexes]
        try:
            target = target.iloc[0]
        except:
            target = 0
        return target

    #Calculate GC percentage

    def analyse_text(text): #Analyse the G/C content of the guide RNA as a percentage.
        count = 0
        letter_count = 0
        try:
            for char in text:
                if char.isalpha(): #Count the length of the guide RNA
                    count += 1
                if char == "C" or char =="G":
                    letter_count += 1 #Count the number of Gs or Cs in the guide RNA
            perc = float(letter_count)/float(count) * 100 #Calculate a percentage of Gs and Cs in the length of the guide RNA
            perc = round(perc, 2) #Round to two decimal places
        except:
            perc = 0
        return perc

    def intconverter(entry): #Convert the entry into an integer
        integer = int(entry)
        return integer

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

    def off_target(search):
        count = 0
        fwcount = orgen.count(search) #Count the occurrence of the guide in the organism genome
        rvcount = revgen.count(search) #Count the occurrence of the guide in the reverse complement of the organism genome
        count = fwcount + rvcount - 1 #Add the counts in both genomes together and subtract 1 to measure only off targets
        return count

    def upstream_real_distance(strand, base_distance): #Determine the actual distance between the 5' base of the amino acid and the upstream guide RNA cut site
        distance = 0
        if strand == "forward":
            distance = base_distance + 1
        if strand == "reverse":
            distance = base_distance + 1
        return distance

    def downstream_real_distance(strand, base_distance): #Determine the actual distance between the 3' base of the amino acid and the downstream guide RNA cut site
        distance = 0
        if strand == "forward":
            distance = base_distance
        if strand == "reverse":
            distance = base_distance
        return distance

    def context(search): #Return the amino acid being target (marked by '*') and the amino acids immediately surrounding it
        x = int(search) - 1
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

    ################################ MODIFY IF ADDITTIONAL PAMS ARE TO BE INCLUDED ################################

    def pamcolumn(entry): #Add a column to the guide dataframe specifying the pam adjacent to the guide RNA generated.
        pam = ''
        if motif == 'NGG' or motif == 'YG': #PAMs at 3' of the guide RNA
            if entry == fiveprime or entry == threeprime:
                pass #Pass this function if there are no guides outputted
            else:
                pam = entry[20:]
        if motif == 'TTTN': #PAM at the 5' of the guide RNA
            if entry == fiveprime or entry == threeprime:
                pass #Pass this function if there are no guides outputted
            else:
                    pam = entry[0:4]
        return pam
        
    def pamidentifier(entry): #Remove the PAM from the guide RNA column
        gRNA = ''
        if motif == 'NGG' or motif == 'YG':
            if entry == fiveprime or entry == threeprime:
                gRNA = entry #Do not apply function if there are no guides.
            else:
                gRNA = entry[0:20]
        if motif == 'TTTN':
            if entry == fiveprime or entry == threeprime:
                gRNA = entry #Do not apply function if there are no guides.
            else:
                gRNA = entry[4:]
        return gRNA

    ###############################################################################################################
    cds = cds.replace("\n", "")
    cds = cds.replace("\r", "") #Remove new lines
    cds = cds.replace(" ", "") #Remove spaces
    remove_lower = lambda text: re.sub('[a-z]', '', text) #Remove any characters that aren't a letter
    cds_exon = remove_lower(cds) #Remove any lower case characters. These correspond to intron sequences so the output is the exon sequence of the gene.

    # Make genomic loci sequence uppercase

    hundredup = hundredup.lower()
    hundreddown = hundreddown.lower()

    cds = hundredup + cds + hundreddown #Append the 100 bases up and downstream to the genomic loci
    cds = cds.replace(" ", "") #Remove spaces
    cds = cds.replace("\n", "")
    cds = cds.replace("\r", "")
    cds_edited = cds.upper()

    orgen = re.sub(r'[^ACTG]', '', orgen) #Remove any characters that are not ACTG
    revgen = reverse_complement(orgen) #Reverse complement of the whole organism genome

    # Reverse translate to get guides on other strand

    cds_reverse = reverse_complement(cds_edited) #Reverse coomplement the genomic loci to identify guides on both strands
      
    # Get positions of exon coding bases

    epos = list(cds)#List of all the bases in the genomic loci
    epospos = pd.DataFrame(epos) #Convert the base list to a dataframe
    epospos = epospos.rename(columns={0: "Base"}) #Rename column to Base
                          
    baseposition(cds) #Get the position of each base in the genomic loci

    epospos['Position'] = psn #Append the base position to the dataframe

    epospos = epospos[epospos['Base'].str.istitle()] #Remove lowercase bases so dataframe only codes for exon

    epospos = epospos.reset_index(drop=True) #Reset the index of the dataframe

    pos1 = epospos[epospos.index % 3 == 0] #Generate a dataframe of every third base from the first base
    pos2 = epospos[epospos.index % 3 == 1] #Generate a dataframe of every third base from the second base
    pos3 = epospos[epospos.index % 3 == 2] #Generate a dataframe of every third base from the third base

    pos1 = list(pos1["Position"]) #Convert the dataframes to a list of positions. These are the positions of the exon coding bases within the context of the genomic loci
    pos2 = list(pos2["Position"])
    pos3 = list(pos3["Position"])

    dna = dna.replace("\n", "")
    dna = dna.replace("\r", "")
    dna = dna.replace(" ", "") #Remove new lines and spcaes from the inputted CDS
    dna = dna.upper()

    # Confirm concatenated exon = gene CDS

    if cds_exon == dna: #Confirm that the inputted CDS sequence matches the CDS identified by the tool. If not the programme stops running.
        print("Inputted CDS and concatenated exons match")
    else:
        print("INPUTTED CDS AND CONCATENATED EXONS DO NOT MATCH")

      
    p = translate(dna) #Translate the CDS identified by the tool into an amino acid sequence using the dictionary of codons.
        
        
        #           GENERAL FUNCTIONS

        # Determine position of specific amino aicds
            
    if aa in p:
        amino_position = charposition(p, aa) #Generate a list of the positions of all the types of one amino acid
    else:
        print ('The amino acid ', aa,' is not present in this sequence.') #Print this string if the amino acid is not present in the amino acid sequence
        
    # Create a 2D array for each amino acid and the respective codon and position of each base

    start = 0
    end = 3
    codon =[]
    while len(dna[start:end])>0: #Split the coding sequence into its respective codons and output them in order to a list called codons
        codon.append(dna[start:end])
        start+=3
        end+=3

    codes = []
    aas = list(p) #List all the amino acids in a protein sequence
    codes = np.array(aas) #Convert the amino acid list into an array
    triple = np.array(codon) #Convert the codon list into an array
    translation = np.column_stack((codes, triple, pos1, pos2, pos3)) #Generate a matrix of the amino acid, the respective codon and the position of each base within the codon

    ################################ ADD ADDITIONAL FUNCTION IF NEW PAMS ARE TO BE INCLUDED. ################################
    # These functions check for PAMs in the inputted genomic loci and its reverse complement.
    #They generate lists of all the potential guides within the genomic loci and a list of all guides on the reverse complement.

    gRNA_list= []#Empty list for gRNAs
    gRNA_list_reverse= []#Empty list for gRNAs

    if motif == "NGG":
        pos = PAMposition(cds_edited)#Identify PAMs in the inputted genomic loci
        pos_reversed = PAMposition(cds_reverse)#Identify PAMs in the reversed genomic loci
        for x in pos:
            entry = cds_edited[x-17:x+6]#Return bases surrounding the cut site
            gRNA_list.append(entry)
        for x in pos_reversed:
            entry = cds_reverse[x-17:x+6]#Return bases surrounding the cut site
            gRNA_list_reverse.append(entry)
            
    if motif == 'YG':
        pos = YGposition(cds_edited)#Identify PAMs in the inputted genomic loci
        pos_reversed = YGposition(cds_reverse)#Identify PAMs in the reversed genomic loci
        for x in pos:
            entry = cds_edited[x-11:x+11]#Return bases surrounding the cut site
            gRNA_list.append(entry)
        for x in pos_reversed:
            entry = cds_reverse[x-11:x+11]#Return bases surrounding the cut site
            gRNA_list_reverse.append(entry)
            
    if motif == 'TTTN':
        pos = TTTNposition(cds_edited)#Identify PAMs in the inputted genomic loci
        pos_reversed = TTTNposition(cds_reverse)#Identify PAMs in the reversed genomic loci
        for x in pos:
            entry = cds_edited[x - 23:x + 8]#Return bases surrounding the cut site
            gRNA_list.append(entry)
        for x in pos_reversed:
            entry = cds_reverse[x - 23:x + 8]#Return bases surrounding the cut site
            gRNA_list_reverse.append(entry)
        
    ###############################################################################################################
    
    # convert the list into a dataframe
        
    gRNA = pd.DataFrame(gRNA_list, columns= ['gRNA Sequence']) #Convert the gRNA list into a dataframe

    gRNA['Position'] = pos #Add a column of guide RNA cut site positions to the dataframe
    gRNA['Strand'] = 'forward' #Add a column to specify these guides are on the forward strand relative to the inputted genomic loci.

    gRNA_reverse = pd.DataFrame(gRNA_list_reverse, columns= ['gRNA Sequence'] ) #Convert  the list of guide RNAs in the reverse complement of the genomic loci into a dataframe
    cds_len = int(len(cds_reverse)) #Store the length of the reverse translated CDS
    gRNA_position = pd.DataFrame(pos_reversed, columns= ['Position'])
    gRNA_reverse['Position'] = cds_len - gRNA_position['Position'] #Put the position of the guide RNA on the reverse strand into the context of the forward strand
    gRNA_reverse['Strand'] = 'reverse'


    gRNA = pd.concat([gRNA, gRNA_reverse]) #Concatenate the forward and reverse guide RNAs
    gRNA['gRNA G/C Content (%)'] = gRNA['gRNA Sequence'].apply(analyse_text) #Calculate the G/C percentage of each guide RNAs
    gRNA = gRNA[gRNA['gRNA G/C Content (%)'] <= 75] #Remove guide RNAs with a G/C percentage above 75%
    gRNA= gRNA.sort_values(by=['Position']) #Sort the guide RNAs by their position in the genomic loci
    gRNA = gRNA.reset_index(drop=True) #Reset the guide RNA dataframe indexes
    del gRNA["gRNA G/C Content (%)"]  #Delete the G/C content column
    
    fiveprime = "No 5' guide could be identified"
    threeprime = "No 3' guide could be identified"
    new_row = pd.DataFrame({'gRNA Sequence':fiveprime, 'Position': 0, 'Strand':''}, index =[0]) 
    gRNA = pd.concat([new_row, gRNA[:]]).reset_index(drop = True) 
    gRNA = gRNA.append({'gRNA Sequence':threeprime, 'Position':len(cds), 'Strand':''}, ignore_index = True)
    #Determine closest gRNAs to specified aa

    base1 = []
    base3 = []
    for x in range(len(translation)):
        if translation[x,0] == aa: #If the amino acid within the translation matrix is the amino acid of interest
            b1 = translation[x,2] #b1 is the position of the most 5' base of the amino acid within the genomic loci
            b3 = translation[x,4] #b3 is the position of the most 3' base of the amino acid within the genomic loci
            base1.append(b1) #Append b1 positions to a list called base1
            base3.append(b3) #Append b3 positions to a list called base3

    #Output gRNA positions to dataframe
    guides = pd.DataFrame(amino_position, columns= ['Amino Acid Position']) #Generate a dataframe called guides where the column is the position of each amino acid of interest in the aminoa cid sequence
    guides["5' Base"] = base1 #New row in the dataframe for the position of the most 5' base in the amino acid codon
    guides["3' Base"] = base3 #New row in the dataframe for the position of the most 3' base in the amino acid codon

    guides["5' gRNA Base Position"] = guides["5' Base"].apply(closest_downstream) #Determine the closest guide RNA cut site to the 5' base of the amino acid codon

    guides["3' gRNA Base Position"] = guides["3' Base"].apply(closest_upstream) #Determine the closest guide RNA cut site to the3' base of the amino acid codon

    guides["5' gRNA Sequence"] = guides["5' gRNA Base Position"].apply(position_finder) #Determine the guide RNA sequence corresponding to the guide RNA position determined by the closest_upstream function
    guides["3' gRNA Sequence"] = guides["3' gRNA Base Position"].apply(position_finder) #Determine the guide RNA sequence corresponding to the guide RNA position determined by the closest_downstream function

    guides["5' PAM"] = guides["5' gRNA Sequence"].apply(pamcolumn) #Generate a new column for the 5' PAM sequence adjacent to the gRNA
    guides["3' PAM"] = guides["3' gRNA Sequence"].apply(pamcolumn) #Generate a new column for the 3' PAM sequence adjacent to the gRNA

    guides["5' gRNA Strand"] = guides["5' gRNA Base Position"].apply(strand_finder) #Determine the strand of the 5' guide RNA
    guides["3' gRNA Strand"] = guides["3' gRNA Base Position"].apply(strand_finder) #Determine the strand of the 3' guide RNA

    guides["5' gRNA Base Position"] = guides["5' gRNA Base Position"].apply(intconverter) #Convert the base position to an integer
    guides["3' gRNA Base Position"] = guides["3' gRNA Base Position"].apply(intconverter) #Convert the base position to an integer
    guides["5' Base"] = guides["5' Base"].apply(intconverter) #Convert the base position to an integer
    guides["3' Base"] = guides["3' Base"].apply(intconverter) #Convert the base position to an integer

    guides["Distance of 5' Cut Site from Amino Acid (bp)"] = guides["5' gRNA Base Position"] - guides["5' Base"] #Determine the crude distance between the 5' guide RNA cut site and the most 5' base of the amino acid codon
    guides["Distance of 3' Cut Site from Amino Acid (bp)"] = guides["3' Base"] - guides["3' gRNA Base Position"] #Determine the crude distance between the 3' guide RNA cut site and the most 3' base of the amino acid codon

    guides["Distance of 3' Cut Site from Amino Acid (bp)"] = - guides["Distance of 3' Cut Site from Amino Acid (bp)"] #Inverse sign of the distance between the amino acid and cut site

    guides["5' gRNA Sequence"] = guides["5' gRNA Sequence"].apply(pamidentifier) #Remove the PAM from the 5' guide RNA column
    guides["3' gRNA Sequence"] = guides["3' gRNA Sequence"].apply(pamidentifier) #Remove the PAM from the 3' guide RNA column

     # Off target search

    guides["5' gRNA Off Target Count"] = guides.apply(lambda row: off_target(row["5' gRNA Sequence"]), axis = 1) #Count the number of guide RNA off targets
    guides["3' gRNA Off Target Count"] = guides.apply(lambda row: off_target(row["3' gRNA Sequence"]), axis = 1) #Count the number of guide RNA off targets

    guides["5' gRNA G/C Content (%)"] = guides["5' gRNA Sequence"].apply(analyse_text) #Calculate the G/C percentage of the 5' guide RNA
    guides["3' gRNA G/C Content (%)"] = guides["3' gRNA Sequence"].apply(analyse_text) #Calculate the G/C percentage of the 3' guide RNA

    guides["5' Notes"] = guides.apply(lambda row: notes(row["5' gRNA Sequence"], row["5' gRNA G/C Content (%)"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
    guides["3' Notes"] = guides.apply(lambda row: notes(row["3' gRNA Sequence"], row["3' gRNA G/C Content (%)"]), axis=1) #Output notes of key 3' guide RNA characterstics to a new column

    del guides["5' Base"] #Delete the column specifying the 5' base position of the amino acid codon
    del guides["3' Base"] #Delete the column specifying the 3' base position of the amino acid codon
        
    guides["Distance of 5' Cut Site from Amino Acid (bp)"] = guides.apply(lambda row: upstream_real_distance(row["5' gRNA Strand"], row["Distance of 5' Cut Site from Amino Acid (bp)"]), axis=1) #Adjust the 5' guide RNA distance to account for strand
    guides["Distance of 3' Cut Site from Amino Acid (bp)"] = guides.apply(lambda row: downstream_real_distance(row["3' gRNA Strand"], row["Distance of 3' Cut Site from Amino Acid (bp)"]), axis=1) #Adjust the 3' guide RNA distance to account for strand

    guides["Context"] = guides["Amino Acid Position"].apply(context) #Generate new column  of the targeted amino acid and its surrounding amino acids

    del guides["5' gRNA Base Position"] #Delete the column specifying the 5' guide RNA cut site position
    del guides["3' gRNA Base Position"] #Delete the column specifying the 3' guide RNA cut site position

    guides[' '] = '' #Empty column in dataframe
     
    guides = guides[["Amino Acid Position", "Context", ' ', "5' gRNA Sequence", "5' PAM", "5' gRNA Strand", "5' gRNA G/C Content (%)", "Distance of 5' Cut Site from Amino Acid (bp)", "5' Notes", "5' gRNA Off Target Count", ' ', "3' gRNA Sequence", "3' PAM", "3' gRNA Strand", "3' gRNA G/C Content (%)", "Distance of 3' Cut Site from Amino Acid (bp)", "3' Notes", "3' gRNA Off Target Count"]] #Reorganise guide dataframe columns for output
    guides.columns = ["Amino Acid Position", "Adjacent amino acids", '                        ', "gRNA Sequence 5' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count", '                        ', "gRNA Sequence 3' of Amino Acid", "PAM", "gRNA Strand", "gRNA G/C Content (%)", "Distance of Cut Site from Amino Acid (bp)", "Notes", "gRNA Off Target Count"] #Rename guide dataframe columns for output
    
    return guides

