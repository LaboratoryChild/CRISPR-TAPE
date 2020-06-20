#!/usr/bin/env python3
"""
Functions shared by specific and general functions
"""

def reverse_complement(dna):
       complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
       return ''.join([complement[base] for base in dna[::-1]])

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
