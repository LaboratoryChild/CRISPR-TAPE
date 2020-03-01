#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 00:32:03 2019

@author: danielanderson
"""

import pandas as pd
import numpy as np
import re

def Specific_function(spec_amino, distance, motif, cds, dna, hundredup, hundreddown, orgen):
    
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
       

################################ ADD FUNCTION IF ADDITTIONAL PAMS ARE TO BE INCLUDED ################################

    def PAMposition(string): #Identify the position of NGGs within the genomic loci. The position returned corresponds to the base 5' of the cut site of the Cas targetting this motif.
        pos = []
        for n in range(len(string) - 1):
            if string[n] == 'G' and string[n+1] == 'G' and n-21 >= 0:
                pos.append(n-4)
        return pos
    
    def YGposition(string): #Identify the position of YGs in the genomic loci (Y is a pyrimidine- C or T). The position returned corresponds to the base 5' of the cut site of the Cas targetting this motif.
        pos = []
        for n in range(len(string) - 1):
            if string[n] == 'C' or string[n] == 'T':
                if string[n + 1] == 'G':
                    pos.append(n-9)
        return pos
    
    def TTTNposition(string): #Identify the position of TTTN motifs in the genomic loci. The position returned corresponds to the base 5' of the cut site of the Cas targetting this motif.
        pos = []
        for n in range(len(string) - 7):
            if string[n] == 'T' and string[n+1] == 'T' and string[n+2] == 'T':
                if string[n+3] == 'A' or string[n+3] == 'C' or string[n+3] == 'G':
                    pos.append(n + 23)
        return pos

###############################################################################################################

    def analyse_text(text): #Analyse the G/C content of the guide RNA as a percentage.
        count = 0
        letter_count = 0
        if text == str(selectionmade):
            perc = ''
        else:
            try:
                for char in text:
                    if char.isalpha(): #Count the length of the guide RNA
                        count += 1
                    if char == "C" or char =="G":
                        letter_count += 1 #Count the number of Gs or Cs in the guide RNA
                perc = float(letter_count)/float(count) * 100 #Calculate a percentage of Gs and Cs in the length of the guide RNA
                perc = round(perc, 2) #Round to two decimal places
            except:
                perc = ''
        return perc
    
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
            if diff > -2 and diff < small: #iterate through the list of guide RNAs and their positions to determine the cut site resistance that is the smallest distance away from the 3' base of the amino acid.
                small = diff
                position = posg
        return position #Returns the cut site position of the closest 3' guide  RNA.
    
    def distanceupfrom(reference): #Return the crude distance between the 5' guide RNA cut site and the base 5' of the amino acid
       difference = up - reference
       return difference
    
    def distancedownfrom(reference): #Return the crude distance between the 3' guide RNA cut site and the base 5' of the amino acid
        difference = reference - down
        return difference
    
    def notes(string, content): #Return key information on the generated guide RNA
        polyt = ''
        if string == str(selectionmade) or content == '':
            polyt = ''
        else:
            for n in range(len(string)-3):
                if string[n:n + 4] == 'TTTT': #Check for four thymines in a row
                    polyt = 'PolyT present. '
            if string[0] != 'G': #Check the guide RNA starts with a 'G' at the most 5' position
                polyt += 'No leading G. '
            if content >= 75:
                polyt += 'G/C content over 75%. ' #Check if the G/C content of the guide is more than or equal to 75%
        return polyt
    
    def real_distance(strand, sign, base_distance): #Adjust the calculated distance between the guide RNA cut site and bases coding for the amino acid depending on which strand the guide is located and whether the guide is 5' or 3' of the amino acid
        distance = 0
        if strand == "forward" and sign >= 0: #If guide is on the sense strand and downstream of the amino acid
            distance = base_distance - 2
        if strand == "forward" and sign <= 0: #If guide is on the sense strand and upstream of the amino acid
            distance = base_distance + 3
        if strand == "reverse" and sign >= 0: #If guide is on the antisense strand and upstream of the amino acid
            distance = base_distance - 2
        if strand == "reverse" and sign <= 0: #If guide is on the antisense strand and downstream of the amino acid
            distance = base_distance + 3
        return distance
    
    def off_target(search):
        fwcount = orgen.count(search) #Count the occurrence of the guide in the organism genome
        rvcount = revgen.count(search) #Count the occurrence of the guide in the reverse complement of the organism genome
        count = fwcount + rvcount - 1 #Add the counts in both genomes together and subtract 1 to measure only off targets
        if search == selectionmade: #If the string in the guide column is the selected amino acid information do not count
            count = ""
        return count

################################ MODIFY IF ADDITTIONAL PAMS ARE TO BE INCLUDED ################################
    def pamcolumn(entry): #Add a column to the guide dataframe specifying the pam adjacent to the guide RNA generated.
        pam = ''
        if entry == selectionmade:
            pam = '' #If guide is the amino acid information then pass
        else:
            if motif == 'NGG' or motif == 'YG': #PAMs at 3' of the guide RNA
                pam = entry[20:]
            if motif == 'TTTN': #PAM at the 5' of the guide RNA
                  if entry == 'No guides within distance range specified':
                      pass #Pass this function if there are no guides outputted
                  else:
                      pam = entry[0:4]
        return pam
    
    def pamidentifier(entry): #Remove the PAM from the guide RNA column
        gRNA = ''
        if entry == str(selectionmade):
            gRNA = selectionmade #Pass if guide RNA is the amino acid selection.
        else:
            if motif == 'NGG' or motif == 'YG':
                if entry == 'No guides within distance range specified':
                      gRNA = noguidestring #Do not apply function if there are no guides.
                else:
                    gRNA = entry[0:20]
            if motif == 'TTTN':
                if entry == 'No guides within distance range specified':
                      gRNA = noguidestring #Do not apply function if there are no guides.
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
    # Create a 2D array for each amino acid and the respective codon and position of each base
    
    start = 0
    end = 3
    codon =[]
    while len(dna[start:end])>0:
        codon.append(dna[start:end]) #Group 3 consecutive bases together to identify the amino acid condons.
        start+=3
        end+=3
    
    codes = []
    aas = list(p) #list the amino acids in the protein sequence
    codes = np.array(aas) #Convert this list to an array
    triple = np.array(codon) #Convert the codon list to an array
    translation = np.column_stack((codes, triple, pos1, pos2, pos3)) #Concatenate these lists together to generate a matrix of the amino acid, the respective condon and each of the 3 base positions in the amino acid.
    
    
    # get base positions of the specific amino acid
    
    selected = translation[spec_amino - 1] #Determine which row of the matrix corresponds to the amino acid of interest
    down = selected[2] #lower base of the codon
    up = selected[4] #upper base of the codon
    selectedaa = selected[0] #The selected amino acid position in the matrix
    selectedpos = spec_amino #Inputted amino acid position
    
    selectionmade = selectedaa + '-' + str(selectedpos)
    print('The amino acid you have selected is ' + selectionmade) #Confirm the selection is the correct amino acid
    
    #           PAM SPECIFIC FUNCTIONS
################################ ADD ADDITIONAL FUNCTION IF NEW PAMS ARE TO BE INCLUDED. ################################
    # These functions check for PAMs in the inputted genomic loci and its reverse complement. They generate lists of all the potential guides within the genomic loci and a list of all guides on the reverse complement.
    gRNA_list= []
    gRNA_list_reverse= []

    if motif == "NGG":
        pos = PAMposition(cds_edited)
        pos_reversed = PAMposition(cds_reverse)
        for x in pos:
            entry = cds_edited[x-17:x+6]
            gRNA_list.append(entry)
        for x in pos_reversed:
            entry = cds_reverse[x-17:x+6]
            gRNA_list_reverse.append(entry)
            
    if motif == 'YG':
        pos = YGposition(cds_edited)
        pos_reversed = YGposition(cds_reverse)
        for x in pos:
            entry = cds_edited[x-11:x+11]
            gRNA_list.append(entry)
        for x in pos_reversed:
            entry = cds_reverse[x-11:x+11]
            gRNA_list_reverse.append(entry)
            
    if motif == 'TTTN':
        pos = TTTNposition(cds_edited)
        pos_reversed = TTTNposition(cds_reverse)
        for x in pos:
            entry = cds_edited[x - 23:x + 8]
            gRNA_list.append(entry)
        for x in pos_reversed:
            entry = cds_reverse[x - 23:x + 8]
            gRNA_list_reverse.append(entry)

###############################################################################################################

    #           GENERAL FUNCTIONS
            
    # convert the lists into dataframes
    
    gRNA = pd.DataFrame(gRNA_list, columns= ['gRNA Sequence']) #Convert the gRNA list into a dataframe
    gRNA['Position'] = pos #Add a column of guide RNA cut site positions to the dataframe
    gRNA['Strand'] = 'forward' #Add a column to specify these guides are on the forward strand relative to the inputted genomic loci.
    
    gRNA_reverse = pd.DataFrame(gRNA_list_reverse, columns= ['gRNA Sequence'] ) #Convert  the list of guide RNAs in the reverse complement of the genomic loci into a dataframe
    cds_len = int(len(cds_reverse)) #Store the length of the reverse translated CDS
    gRNA_position = pd.DataFrame(pos_reversed, columns= ['Position'])
    gRNA_reverse['Position'] = cds_len - gRNA_position['Position'] #Put the position of the guide RNA on the reverse strand into the context of the forward strand
    gRNA_reverse['Strand'] = 'reverse'
    
    gRNA = pd.concat([gRNA, gRNA_reverse]) #Concatenate the forward and reverse guide RNAs
    gRNA= gRNA.sort_values(by=['Position']) #Sort the guide RNAs by their position in the genomic loci
    
    gRNA = gRNA.reset_index(drop=True) #Reset gRNA dataframe indexes
    
     # Guides around specified position
        
    upposition = closest_upstream(up) #The position correspdonding to the closest 5' guide RNA
    downposition = closest_downstream(down)#The position correspdonding to the closest 3' guide RNA
    
    downstream = np.where(gRNA['Position'] == upposition) #The index of the guide RNA in the list with the corresponding guide RNA position
    d = downstream[0]
    d = int(d[0])
    downerguides = gRNA.iloc[d:-1, :] #Split the guide RNA dataframe into two depending on the closest downstream guides
    
    upstream = np.where(gRNA['Position'] == downposition) #The index of the guide RNA in the list with the corresponding guide RNA position
    u = upstream[0]
    u = int(u[0])
    upperguides = gRNA.iloc[0:u+1, :] #Split the guide RNA dataframe into two depending on the closest upstream guides
    
    up = int(up)
    down = int(down)
    
    upperguides["Distance from Amino Acid (bp)"] = upperguides['Position'].apply(distanceupfrom) #The crude distance between the guide RNA cut site and the most 5' base of the amino acid
    
    downerguides["Distance from Amino Acid (bp)"] = downerguides['Position'].apply(distancedownfrom) #The crude distance between the guide RNA cut site and the most 3' base of the amino acid
    
    downerguides = downerguides[downerguides['gRNA Sequence'] != ''] #Delete guide RNA elements that are empty
    upperguides = upperguides[upperguides['gRNA Sequence'] != ''] #Delete guide RNA elements that are empty
   
    #           GENERAL FUNCTIONS
    
    noguidestring = "No guides within distance range specified"

    downerguides = downerguides[downerguides["Distance from Amino Acid (bp)"] <= distance] #Remove guides over the inputted maximum guide distance
    upperguides = upperguides[upperguides["Distance from Amino Acid (bp)"] <= distance]  #Remove guides over the inputted maximum guide distance
    
    noguides = pd.DataFrame(columns=["gRNA Sequence", "Strand", "G/C Content (%)", "Distance from Amino Acid (bp)", "Notes"]) #Generate a new dataframe for cases where there are no guides within the specified distance
    noguides = noguides.append({"gRNA Sequence": noguidestring, "Strand":"","G/C Content (%)":"","Distance from Amino Acid (bp)":"","Notes":""}, ignore_index=True) #Generate a new dataframe for cases where there are no guides within the specified distance
    
    upperguides = upperguides.reset_index(drop=True) #Reset upstream guide RNA indexes
    downerguides = downerguides.reset_index(drop=True) #Reset downstream guide RNA indexes
    
    if len(upperguides) == 0:
        pass
    else:
        upperguides["Distance from Amino Acid (bp)"] = - upperguides["Distance from Amino Acid (bp)"] #Make distance of upstream guides negative
    
    #           SPECIFIC FUNCTIONS
    
    if len(downerguides) == 0:
        pass
    else:
        downerguides["Distance from Amino Acid (bp)"] = downerguides.apply(lambda row: real_distance(row["Strand"], row["Distance from Amino Acid (bp)"], row["Distance from Amino Acid (bp)"]), axis=1) #Adjust guide RNA distance to account for strand and whether upstream or downstream
    
    if len(upperguides) == 0:
        pass
    else:
        upperguides["Distance from Amino Acid (bp)"] = upperguides.apply(lambda row: real_distance(row["Strand"], row["Distance from Amino Acid (bp)"], row["Distance from Amino Acid (bp)"]), axis=1) #Adjust guide RNA distance to account for strand and whether upstream or downstream
     
    # Get names of indexes for which column has value greater than distance
    
    downerguides = downerguides[(downerguides['Distance from Amino Acid (bp)'] >= -1)]  #Remove guides if they are less than 0 base pairs away from the amino acid
    upperguides = upperguides[(upperguides['Distance from Amino Acid (bp)'] <= 1)]  #Remove guides if they are more than 0 base pairs away from the amino acid
    
    downerguides = downerguides.sort_values(by=['Distance from Amino Acid (bp)']) #Arrange guide RNAs by their distance from the amino acid
    upperguides = upperguides.sort_values(by = ['Distance from Amino Acid (bp)']) #Arrange guide RNAs by their distance from the amino acid
    
    if len(downerguides) == 0:
        downerguides = noguides
  
    if len(upperguides) == 0:
        upperguides = noguides
        
    amino_acid = pd.DataFrame(columns=["gRNA Sequence", "Strand", "G/C Content (%)", "Distance from Amino Acid (bp)", "Notes"]) #Generate a new amino acid for the amino acid target information
    amino_acid = amino_acid.append({"gRNA Sequence": selectedaa + '-' + str(selectedpos), "Strand":"","G/C Content (%)":"","Distance from Amino Acid (bp)":"","Notes":""}, ignore_index=True) #Generate a new amino acid for the amino acid target information
    guides = pd.concat([upperguides, amino_acid, downerguides]) #Concatenate upstream, amino acid information and downstream dataframes
        
    guides["PAM"] = guides["gRNA Sequence"].apply(pamcolumn) #Create new column of the specific PAM for each guide RNA
    guides["gRNA Sequence"] = guides["gRNA Sequence"].apply(pamidentifier) #Remove the PAM from the guide RNA column
    
    guides['Off Target Count'] = guides.apply(lambda row: off_target(row['gRNA Sequence']), axis =1) #Count the number of guide RNA off targets

    #Calculate GC percentage

    guides['G/C Content (%)'] = guides['gRNA Sequence'].apply(analyse_text) #Add a new column to the gRNA dataframe with G/C percentages
    guides["Notes"] = guides.apply(lambda row: notes(row["gRNA Sequence"], row["G/C Content (%)"]), axis=1) #Output notes of key guide RNA characterstics to a new column

    guides = guides[['Distance from Amino Acid (bp)', 'gRNA Sequence', 'PAM', 'Strand', 'G/C Content (%)', 'Off Target Count', 'Notes']] #Reorganise the guide RNA dataframe
    
    
    return guides
