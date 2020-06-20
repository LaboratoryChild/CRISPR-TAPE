#!/usr/bin/env python3
"""
PAM-specific functions
"""

################################ ADD FUNCTION IF ADDITTIONAL PAMS ARE TO BE INCLUDED ################################

# get index position of all PAMs within the genomic loci
def PAMposition(stringm motif): #Identify the position of NGGs within the genomic loci
    position = [] #Empty list to store identified gRNAs
    entry = []
    if motif == "NGG":
        for n in range(len(string) - 1):
            if string[n] == 'G' and string[n+1] == 'G' and n-21 >= 0: #If two Gs in a row and Gs not at the start of the genomic loci
                pos = n - 4
                position.append(pos)
                entry.append(cds[pos-17:pos+6])
    elif motif == "YG":
        for n in range(len(string) - 1):
            if string[n] == 'C' or string[n] == 'T' and string[n + 1] == 'G': #If a C or T is followed by G
                    pos = n-9  #Append the position of the base 5' of the cut site
                    position.append(pos)
                    entry.append(cds[pos-11:pos+11])
    elif motif == "TTTN":
        for n in range(len(string) - 7):
            if string[n] == 'T' and string[n+1] == 'T' and string[n+2] == 'T': #If 3 Ts in a row
                if string[n+3] == 'A' or string[n+3] == 'C' or string[n+3] == 'G': #If 4th base is not T
                    pos = n + 23 #Append the position of the base 5' of the cut site
                    position.append(pos)
                    entry.append(cds[pos-23:pos+8])
    return position, entry
###############################################################################################################

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
