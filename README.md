# CRISPR-TAPE
A protein-centric CRISPR gRNA design tool for TArgeted Protein Engineering

# Motivation 
Existing CRIPSR gRNA design tools target protein-coding regions within genomic loci and non-specifically target the entire input region of DNA. Current tools fail to consider proteomic-based applications, so CRISPR-TAPE has been developed to reduce the substantial time burden associated with manual curation of gRNA libraries and empower the proteomics community. 

# Version
CRISPR-TAPE version 1.0.0

# Quick Start 

    pip install CRISPR-TAPE

General_function:

    from CRISPR_TAPE import General_function
  
    gRNAs = General_function(aa, motif, cds, dna, hundredup, hundreddown, orgen)

Specific_function:

    from CRISPR_TAPE import Specific_function
        
    gRNAs = Specific_function(spec_amino, distance, motif, cds, dna, hundredup, hundreddown, orgen)
    
# Inputs
General_function:

• aa = String of the short letter code of all amino acids to be targeted

• motif = String of the protospacer adjacent motif (PAM) sequence

• cds = String of the coding sequence for the protein of interest

• hundredup = String of the 5' UTR or additional upstream sequences if desired

• hundreddown = String of the 3' UTR or additional downstream sequences if desired

• orgen = String of the forward genomic loci sequence of the organism of interest


Specific_function:
• spec_amino = Integer residue position of the residue in the protein sequence

• distance = Integer distance of the maximum distance for gRNA identification

• motif = String of the protospacer adjacent motif (PAM) sequence

• cds = String of the coding sequence for the protein of interest

• hundredup = String of the 5' UTR or additional upstream sequences if desired

• hundreddown = String of the 3' UTR or additional downstream sequences if desired

• orgen = String of the forward genomic loci sequence of the organism of interest
 
# Output 

General_function:
A dataframe of 5' and 3' gRNAs immediately upstream and downstream of the specified amino acids and some basic properties including the off target count, G/C content and presence of a leading G seqence and polyA sequence. 

Specific_function:
A dataframe of 5' and 3' gRNAs immediately upstream and downstream of the specified amino acids and some basic properties including the off target count, G/C content and presence of a leading G seqence and polyA sequence. 

# Technologies 
The CRISPR-TAPE python scripts are compatible with Python 3.6 and later versions. Outside of the standard library, it makes use of the following packages: numpy version 1.18.0, pandas version 0.25.3.

# General_functions PAM-Specific Functions
CRISPR-TAPE is an open source python programme. The modularity of CRISPR-TAPE allows for easy incorporation of additional PAM sequences by modifying the code shown below:
   
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

    def pamcolumn(entry):   #Generate a new column containing the PAM adjacent to the guide RNA
        if motif == 'NGG' or motif == 'YG':
            pam = entry[20:] #If NGG or YG PAM is at the end of the identified gRNA
        if motif == 'TTTN':
            pam = entry[0:4] #If TTTN PAM is at the beginning of the identified gRNA
        return pam

    #Remove the PAM sequence from the gRNA column
    def pamidentifier(entry):
        if motif == 'NGG' or motif == 'YG':
            gRNA = entry[0:20] #If NGG or YG remove final 3 bases
        if motif == 'TTTN':
            gRNA = entry[4:] #If TTTN remove initial 4 bases
        return gRNA
    
    if motif == "NGG":
        pos = PAMposition(cds_edited)#Identify PAMs in the inputted genomic loci
        pos_reversed = PAMposition(cds_reverse)#Identify PAMs in the reversed genomic loci
        gRNA_list= []#Empty list for gRNAs
        for x in pos:
            entry = cds_edited[x-17:x+6]#Return bases surrounding the cut site
            gRNA_list.append(entry)
        gRNA_list_reverse= []#Empty list for gRNAs
        for x in pos_reversed:
            entry = cds_reverse[x-17:x+6]#Return bases surrounding the cut site
            gRNA_list_reverse.append(entry)
            
    if motif == 'YG':
        pos = YGposition(cds_edited)#Identify PAMs in the inputted genomic loci
        pos_reversed = YGposition(cds_reverse)#Identify PAMs in the reversed genomic loci
        gRNA_list= []#Empty list for gRNAs
        for x in pos:
            entry = cds_edited[x-11:x+11]#Return bases surrounding the cut site
            gRNA_list.append(entry)
        gRNA_list_reverse= []#Empty list for gRNAs
        for x in pos_reversed:
            entry = cds_reverse[x-11:x+11]#Return bases surrounding the cut site
            gRNA_list_reverse.append(entry)
            
    if motif == 'TTTN':
        pos = TTTNposition(cds_edited)#Identify PAMs in the inputted genomic loci
        pos_reversed = TTTNposition(cds_reverse)#Identify PAMs in the reversed genomic loci
        gRNA_list= []#Empty list for gRNAs
        for x in pos:
            entry = cds_edited[x - 23:x + 8]#Return bases surrounding the cut site
            gRNA_list.append(entry)
        gRNA_list_reverse= []#Empty list for gRNAs
        for x in pos_reversed:
            entry = cds_reverse[x - 23:x + 8]#Return bases surrounding the cut site
            gRNA_list_reverse.append(entry)

# Specific_function PAM-Specific Functions
CRISPR-TAPE is an open source python programme. The modularity of CRISPR-TAPE allows for easy incorporation of additional PAM sequences by modifying the code shown below:
   
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
    
     if motif == "NGG":
        pos = PAMposition(cds_edited)
        pos_reversed = PAMposition(cds_reverse)
        gRNA_list= []
        for x in pos:
            entry = cds_edited[x-17:x+6]
            gRNA_list.append(entry)
        gRNA_list_reverse= []
        for x in pos_reversed:
            entry = cds_reverse[x-17:x+6]
            gRNA_list_reverse.append(entry)
            
    if motif == 'YG':
        pos = YGposition(cds_edited)
        pos_reversed = YGposition(cds_reverse)
        gRNA_list= []
        for x in pos:
            entry = cds_edited[x-11:x+11]
            gRNA_list.append(entry)
        gRNA_list_reverse= []
        for x in pos_reversed:
            entry = cds_reverse[x-11:x+11]
            gRNA_list_reverse.append(entry)
            
    if motif == 'TTTN':
        pos = TTTNposition(cds_edited)
        pos_reversed = TTTNposition(cds_reverse)
        gRNA_list= []
        for x in pos:
            entry = cds_edited[x - 23:x + 8]
            gRNA_list.append(entry)
        gRNA_list_reverse= []
        for x in pos_reversed:
            entry = cds_reverse[x - 23:x + 8]
            gRNA_list_reverse.append(entry)
            
# Troubleshooting
• All inputs are case sensitive and it is important to adhere to the requirements specified adjacent to each entry box

• gRNAs will be incorrect if exonic bases are not capitalised and intronic/untranslated bases are not lowercase in the genomic loci input

• The programme currently only recognises “A”, “T”, “C” and “G” bases and will filter out any other characters

• If no protospacer adjacent motif is specified, the programme will not run and no guide RNAs will be outputted

• The programme will stop running and a pop-up prompt will open if the inputted protein coding sequence and the exon sequences from the inputted genomic loci do not match

• Off target counts will display as “-1” if gRNAs are not found within the forward or reverse complement of the organism genome. If this occurs, ensure the genomic loci and genome originate from the same organism

• If some amino acids are missing from the output, this is because they cannot be target using the current information. If this is the case, add more upstream and downstream bases to the input requesting 100 bases upstream and downstream of the genomic loci

# Authors
CRISPR-TAPE was developed by Daniel Anderson, Henry Benns and Dr Matthew Child
