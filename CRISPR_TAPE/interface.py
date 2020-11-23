from tkinter import *
from .specific_function import specific_function
from .general_function import general_function

def main():

    def retrieve_spec(): #Function to get the values of the entry boxes and open the run completed pop-up
        try:
            dna_seq = gen.get()
            dna_seq = str(dna_seq)
            uputr = hup.get()
            downutr = hdn.get()
            uputr = str(uputr)
            downutr = str(downutr)
            motif = thepam.get()
            coding_seq = dn.get()
            coding_seq = str(coding_seq)
            filename = fn.get()
            filename = str(filename)
            spec_amino = a.get()
            spec_amino = int(spec_amino)
            distance = d.get()
            if distance == '0' or distance == "" or distance == " ":
                distance = 1000000
            else:
                distance = int(distance)

            # import organism genome
            with open("INSERT_ORGANISM_GENOME_HERE.txt", "r") as f:
                orgen = f.read()

            guides = specific_function(spec_amino, 
                                distance,
                                motif, 
                                dna_seq, 
                                coding_seq, 
                                uputr, 
                                downutr, 
                                orgen)
            csvfilename = str(filename + ".csv")
            guides.to_csv(csvfilename, index=False)

            done = Tk()
            Label(done, text='Your guide RNAs have successfully been generated.', font='Helvetica 12').grid(row=1)
            Label(done, text='To view them, open the "' + filename + '.csv" file in the CRISPR-TAPE folder.', font='Helvetica 12').grid(row=2)
            done.title("**GUIDES SUCCESSFULLY GENERATED**")
            done.lift()
            done.attributes('-topmost',True)
            done.focus_force()
            done.bind('<FocusIn>', OnFocusIn)
            done.mainloop() 
        except:
            general = Tk()
            Label(general, text='Guides have not been successfully generated.', font='Helvetica 12').grid(row=1)
            Label(general, text='Please check the toubleshooting section of the README file.', font='Helvetica 12').grid(row=2)
            general.title("**ERROR**")
            general.lift()
            general.attributes('-topmost',True)
            general.focus_force()
            general.bind('<FocusIn>', OnFocusIn)
            general.mainloop() 
        
        return 

    def retrieve_gen(): #Function to get the values of the entry boxes and open the run completed pop-up
        try:
            dna_seq = gen.get()
            dna_seq = str(dna_seq)
            uputr = hup.get()
            downutr = hdn.get()
            uputr = str(uputr)
            downutr = str(downutr)
            motif = thepam.get()
            coding_seq = dn.get()
            coding_seq = str(coding_seq)
            filename = fn.get()
            filename = str(filename)
            distance = 10000000
            aa = v.get()
            aa = str(aa)

            # import organism genome
            with open("INSERT_ORGANISM_GENOME_HERE.txt", "r") as f:
                orgen = f.read()

            guides = general_function(aa, 
                                    motif, 
                                    dna_seq, 
                                    coding_seq, 
                                    uputr, 
                                    downutr, 
                                    orgen)
            csvfilename = str(filename + ".csv")
            guides.to_csv(csvfilename, index=False)

            done = Tk()
            Label(done, text='Your guide RNAs have successfully been generated.', font='Helvetica 12').grid(row=1)
            Label(done, text='To view them, open the "' + filename + '.csv" file in the CRISPR-TAPE folder.', font='Helvetica 12').grid(row=2)
            done.title("**GUIDES SUCCESSFULLY GENERATED**")
            done.lift()
            done.attributes('-topmost',True)
            done.focus_force()
            done.bind('<FocusIn>', OnFocusIn)
            done.mainloop() 
        except:
            general = Tk()
            Label(general, text='Guides have not been successfully generated.', font='Helvetica 12').grid(row=1)
            Label(general, text='Please check the toubleshooting section of the README file.', font='Helvetica 12').grid(row=2)
            general.title("**ERROR**")
            general.lift()
            general.attributes('-topmost',True)
            general.focus_force()
            general.bind('<FocusIn>', OnFocusIn)
            general.mainloop() 

        return

    def instructions(): #Function to open the instructions for the programme 
        root = Tk()
        root.title("CRISPR-TAPE INSTRUCTIONS")
        Label(root, text='To use the programme, input the genome sequence of your organism of interest into the "INSERT_ORGANISM_GENOME_HERE.txt"', font='Helvetica 12').grid(row=1, columnspan = 2, sticky=W)
        Label(root, text='file, then copy and paste your genomic loci and coding sequences in the entry boxes below. UTRs are not strictly required', font='Helvetica 12').grid(row=2, columnspan = 2, sticky=W)
        Label(root, text='(see README). Once you have selected your parameters, choose your PAM sequence then choose the CRISPR-TAPE mode', font='Helvetica 12').grid(row=4, columnspan = 2, sticky=W)
        Label(root, text='depending on whether you want to target a specific amino acid ("Run Option 1") or all amino acids of one type ("Run Option 2").', font='Helvetica 12').grid(row=5, columnspan = 2, sticky=W)
        Label(root, text='When the programme has finished running, your guides will be saved to the CRISPR-TAPE folder as a .csv file. If you would ', font='Helvetica 12').grid(row=6, columnspan = 2, sticky=W)
        Label(root, text='like to run the programme again, press "RESET" to clear all fields. For more information, view the README file.', font='Helvetica 12').grid(row=7, columnspan = 2, sticky=W)
        root.lift()
        root.attributes('-topmost',True)
        root.focus_force()
        root.bind('<FocusIn>', OnFocusIn)
        root.mainloop()
        return

    def reset(): #Function to reset the entries in the GUI window 
        gen.set("")
        hup.set("")
        hdn.set("")
        dn.set("")
        fn.set("")
        v.set("")
        a.set("0")
        d.set("0")
        return

    def OnFocusIn(event):
        if type(event.widget).__name__ == 'Tk':
            event.widget.attributes('-topmost', False)
        return

    self = Tk()
    Label(self, text="Welcome to CRISPR-TAPE, a CRISPR gRNA design tool for Targeted Protein Engineering. CRISPR-TAPE is a python-", font='Helvetica 12').grid(row=0, columnspan = 2, sticky=W)
    Label(self, text="based programme that outputs gRNAs designed to target Cas9 to specified residues or amino acid types.", font='Helvetica 12').grid(row=1, columnspan = 2, sticky=W)
    Button(self, text="INSTRUCTIONS", command = lambda: instructions()).grid(row=2, sticky=W)

    Label(self, text="").grid(row=3)
    Label(self, text="Name of guide output file (no spaces):" , font='Helvetica 12 bold').grid(row=4, sticky = W)
    fn = StringVar()
    Entry(self, textvariable = fn).grid(row=4, column=1, sticky=W)    
    Label(self, text="Please input the genomic loci sequence of your protein here:", font='Helvetica 12 bold').grid(row=5, column = 0, sticky=W)
    Label(self, text="(UTRs and introns lowercase, exons uppercase)", font='Helvetica 12').grid(row=6, column = 0, sticky=W)
    Label(self, text="", font='Helvetica 12 bold').grid(row=7, column = 0, sticky=W)
    gen = StringVar()
    Entry(self, textvariable = gen).grid(row=5, column=1, sticky=W)
    Label(self, text="If your genomic loci does not include 5' and 3' UTRs:", font='Helvetica 12 bold').grid(row=8, sticky=W)
    Label(self, text="Please input 5' UTR here:", font='Helvetica 12').grid(row=9, sticky=W)
    Label(self, text="Please input 3' UTR here:", font='Helvetica 12').grid(row=10, sticky=W)
    hup = StringVar()
    hdn = StringVar()
    Entry(self, textvariable = hup).grid(row=9, column=1, sticky=W)
    Entry(self, textvariable = hdn).grid(row=10, column=1, sticky=W)
    hup.set(' ')
    hdn.set(' ')
    thepam = StringVar()
    Label(self, text="").grid(row=11, sticky=W)
    Label(self, text="Please input the coding sequence of your protein here (no UTRs):", font='Helvetica 12 bold').grid(row=12, sticky = W)
    dn = StringVar()
    Entry(self, textvariable = dn).grid(row=12, column=1, sticky=W)
    Label(self, text="").grid(row=13, sticky = W)
    Label(self, text="Specify nuclease protospacer adjacent motif (PAM):", font='Helvetica 12 bold').grid(row=14, sticky=W)
    Radiobutton(self, text="NGG", value="NGG", variable=thepam, font='Helvetica 12').grid(row=15, sticky=W)
    Radiobutton(self, text="YG", value="YG", variable=thepam, font='Helvetica 12').grid(row=16, sticky=W)
    Radiobutton(self, text="TTTN", value="TTTN", variable=thepam, font='Helvetica 12').grid(row=17, sticky=W)
    Label(self, text="").grid(row=18, sticky = W)
    Label(self, text="Now choose to target a specific amino acid or target all amino acids of one type:", font='Helvetica 12').grid(row=19, sticky=W)
    Label(self, text="").grid(row=20, sticky = W)

    Label(self, text="OPTION 1: Target a specific amino acid (e.g. cysteine at position 106 = 106):", font='Helvetica 12 bold').grid(row=21, sticky = W)
    a = IntVar()
    Label(self, text="Input the residue position of this amino acid:", font='Helvetica 12').grid(row=22, sticky=W)
    Entry(self, textvariable = a).grid(row=22, column=1, sticky=W)
    d = StringVar()
    d.set("")
    Label(self, text="Please specify the maximum guide distance (in nucleotides) from the amino acid:", font='Helvetica 12').grid(row=23, sticky=W)
    Entry(self, textvariable = d).grid(row=23, column=1, sticky=W)

    Button(self, text="Run Option 1", command = lambda: retrieve_spec()).grid(row=24, column =1, sticky=W)

    v = StringVar()
    Label(self, text="OPTION 2: Target all amino acids of a certain type (e.g. cysteine = C):", font='Helvetica 12 bold').grid(row = 25, sticky = W)
    Label(self, text="Input your target amino acid single letter code:", font='Helvetica 12').grid(row=26, sticky=W)
    Entry(self, textvariable = v).grid(row=26, column=1, sticky=W)
    Button(self, text="Run Option 2", command = lambda: retrieve_gen()).grid(row=27, column =1, sticky=W)

    Label(self, text="").grid(row=29, column =1, sticky=W)
    Button(self, text="RESET", command= lambda: reset()).grid(row=30, column =1, sticky=W)
    self.title("CRISPR-TAPE")
    self.lift()
    self.attributes('-topmost',True)
    self.focus_force()
    self.bind('<FocusIn>', OnFocusIn)

    self.mainloop() 
    return

if __name__ == '__main__':
    main()

    sys.exit(0)
