from CRISPR_TAPE.general_function import general_function
from CRISPR_TAPE.specific_function import specific_function
import multiprocessing
import os
import sys
from tkinter import *

multiprocessing.freeze_support()

def retrieve_spec(gen, hup, hdn, thepam, dn, fn, a, d, gf): #Function to get the values of the entry boxes and open the run completed pop-up
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
        genome_file = gf.get()
        spec_amino = a.get()
        spec_amino = int(spec_amino)
        distance = d.get()
        if distance == '0' or distance == "" or distance == " ":
            distance = 1000000
        else:
            distance = int(distance)

        with open(os.path.abspath(genome_file), "r") as f:
            reference_genome = f.read()

        guides = specific_function(spec_amino,
                            distance,
                            motif,
                            dna_seq,
                            coding_seq,
                            uputr,
                            downutr,
                            reference_genome)
        if not ".csv" in filename:
            csvfilename = str(filename + ".csv")
        else:
            csvfilename = filename
        output_path = os.path.join(os.getcwd(), csvfilename)
        guides.to_csv(output_path, index=False)

        done = Tk()
        Label(done, text='Your guide RNAs have successfully been generated.', font='Helvetica 12').grid(row=1)
        Label(done, text='To view them, open the "' + output_path + '" file in the CRISPR-TAPE folder.', font='Helvetica 12').grid(row=2)
        done.title("**GUIDES SUCCESSFULLY GENERATED**")
        done.lift()
        done.attributes('-topmost',True)
        done.focus_force()
        done.bind('<FocusIn>', OnFocusIn)
        done.mainloop()
    except Exception as e:
        general = Tk()
        Label(general, text='Guides have not been successfully generated.', font='Helvetica 12').grid(row=1)
        Label(general, text='Please check the toubleshooting section of the README file.', font='Helvetica 12').grid(row=2)
        Label(general, text='Error is:', font='Helvetica 12').grid(row=3)
        Label(general, text=e, font='Helvetica 12').grid(row=4)
        general.title("**ERROR**")
        general.lift()
        general.attributes('-topmost',True)
        general.focus_force()
        general.bind('<FocusIn>', OnFocusIn)
        general.mainloop()

def retrieve_gen(gen, hup, hdn, thepam, dn, fn, v, gf): #Function to get the values of the entry boxes and open the run completed pop-up
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
        genome_file = gf.get()
        aa = v.get()
        aa = str(aa)

        # import organism genome
        with open(os.path.abspath(genome_file), "r") as f:
            reference_genome = f.read()

        guides = general_function(aa,
                                motif,
                                dna_seq,
                                coding_seq,
                                uputr,
                                downutr,
                                reference_genome)
        if not ".csv" in filename:
            csvfilename = str(filename + ".csv")
        else:
            csvfilename = filename
        output_path = os.path.join(os.getcwd(), csvfilename)
        guides.to_csv(output_path, index=False)

        done = Tk()
        Label(done, text='Your guide RNAs have successfully been generated.', font='Helvetica 12').grid(row=1)
        Label(done, text='To view them, open the "' + output_path + '" file in the CRISPR-TAPE folder.', font='Helvetica 12').grid(row=2)
        done.title("**GUIDES SUCCESSFULLY GENERATED**")
        done.lift()
        done.attributes('-topmost',True)
        done.focus_force()
        done.bind('<FocusIn>', OnFocusIn)
        done.mainloop()
    except Exception as e:
        general = Tk()
        Label(general, text='Guides have not been successfully generated.', font='Helvetica 12').grid(row=1)
        Label(general, text='Please check the toubleshooting section of the README file.', font='Helvetica 12').grid(row=2)
        Label(general, text='Error is:', font='Helvetica 12').grid(row=3)
        Label(general, text=e, font='Helvetica 12').grid(row=4)
        general.title("**ERROR**")
        general.lift()
        general.attributes('-topmost',True)
        general.focus_force()
        general.bind('<FocusIn>', OnFocusIn)
        general.mainloop()

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

def reset(gen, hup, hdn, dn, fn, v, a, d, gf): #Function to reset the entries in the GUI window
    gen.set("")
    hup.set("")
    hdn.set("")
    dn.set("")
    fn.set("")
    v.set("")
    a.set("0")
    d.set("0")
    gf.set("")

def OnFocusIn(event):
    if type(event.widget).__name__ == 'Tk':
        event.widget.attributes('-topmost', False)

def main():
    self = Tk()
    Label(self, text="Welcome to CRISPR-TAPE, a CRISPR gRNA design tool for Targeted Protein Engineering. CRISPR-TAPE is a python-", font='Helvetica 12').grid(row=0, columnspan = 2, sticky=W)
    Label(self, text="based programme that outputs gRNAs designed to target Cas9 to specified residues or amino acid types.", font='Helvetica 12').grid(row=1, columnspan = 2, sticky=W)
    Button(self, text="INSTRUCTIONS", command = lambda: instructions()).grid(row=2, sticky=W)

    Label(self, text="").grid(row=3)
    Label(self, text="Name of guide output file (no spaces):" , font='Helvetica 12 bold').grid(row=4, sticky = W)
    Label(self, text="Genome assembly file of the organism of interest:" , font='Helvetica 12 bold').grid(row=5, sticky = W)
    fn = StringVar()
    Entry(self, textvariable = fn).grid(row=4, column=1, sticky=W)
    gf = StringVar()
    Entry(self, textvariable = gf).grid(row=5, column=1, sticky=W)
    Label(self, text="Please input the genomic loci sequence of your protein here:", font='Helvetica 12 bold').grid(row=6, column = 0, sticky=W)
    Label(self, text="(UTRs and introns lowercase, exons uppercase)", font='Helvetica 12').grid(row=7, column = 0, sticky=W)
    Label(self, text="", font='Helvetica 12 bold').grid(row=8, column = 0, sticky=W)
    gen = StringVar()
    Entry(self, textvariable = gen).grid(row=6, column=1, sticky=W)
    Label(self, text="If your genomic loci does not include 5' and 3' UTRs:", font='Helvetica 12 bold').grid(row=9, sticky=W)
    Label(self, text="Please input 5' UTR here:", font='Helvetica 12').grid(row=10, sticky=W)
    Label(self, text="Please input 3' UTR here:", font='Helvetica 12').grid(row=11, sticky=W)
    hup = StringVar()
    hdn = StringVar()
    Entry(self, textvariable = hup).grid(row=10, column=1, sticky=W)
    Entry(self, textvariable = hdn).grid(row=11, column=1, sticky=W)
    hup.set(' ')
    hdn.set(' ')
    thepam = StringVar()
    Label(self, text="").grid(row=12, sticky=W)
    Label(self, text="Please input the coding sequence of your protein here (no UTRs):", font='Helvetica 12 bold').grid(row=13, sticky = W)
    dn = StringVar()
    Entry(self, textvariable = dn).grid(row=13, column=1, sticky=W)
    Label(self, text="").grid(row=14, sticky = W)
    Label(self, text="Specify nuclease protospacer adjacent motif (PAM):", font='Helvetica 12 bold').grid(row=15, sticky=W)
    Radiobutton(self, text="NGG", value="NGG", variable=thepam, font='Helvetica 12').grid(row=16, sticky=W)
    Radiobutton(self, text="YG", value="YG", state=DISABLED, variable=thepam, font='Helvetica 12').grid(row=17, sticky=W)
    Radiobutton(self, text="TTTN", value="TTTN", state=DISABLED, variable=thepam, font='Helvetica 12').grid(row=18, sticky=W)
    Label(self, text="").grid(row=19, sticky = W)
    Label(self, text="Now choose to target a specific amino acid or target all amino acids of one type:", font='Helvetica 12').grid(row=20, sticky=W)
    Label(self, text="").grid(row=21, sticky = W)

    Label(self, text="OPTION 1: Target a specific amino acid (e.g. cysteine at position 106 = 106):", font='Helvetica 12 bold').grid(row=22, sticky = W)
    a = IntVar()
    Label(self, text="Input the residue position of this amino acid:", font='Helvetica 12').grid(row=23, sticky=W)
    Entry(self, textvariable = a).grid(row=23, column=1, sticky=W)
    d = StringVar()
    d.set("")
    Label(self, text="Please specify the maximum guide distance (in nucleotides) from the amino acid:", font='Helvetica 12').grid(row=24, sticky=W)
    Entry(self, textvariable = d).grid(row=24, column=1, sticky=W)

    Button(self, text="Run Option 1", command = lambda: retrieve_spec(gen, hup, hdn, thepam, dn, fn, a, d, gf)).grid(row=25, column =1, sticky=W)

    v = StringVar()
    Label(self, text="OPTION 2: Target all amino acids of a certain type (e.g. cysteine = C):", font='Helvetica 12 bold').grid(row = 26, sticky = W)
    Label(self, text="Input a target amino acid single letter code or motif:", font='Helvetica 12').grid(row=27, sticky=W)
    Entry(self, textvariable = v).grid(row=27, column=1, sticky=W)
    Button(self, text="Run Option 2", command = lambda: retrieve_gen(gen, hup, hdn, thepam, dn, fn, v, gf)).grid(row=28, column =1, sticky=W)

    Label(self, text="").grid(row=30, column =1, sticky=W)
    Button(self, text="RESET", command= lambda: reset(gen, hup, hdn, dn, fn, v, a, d, gf)).grid(row=31, column =1, sticky=W)
    self.title("CRISPR-TAPE")
    self.lift()
    self.attributes('-topmost',True)
    self.focus_force()
    self.bind('<FocusIn>', OnFocusIn)

    self.mainloop()
    sys.exit(0)

if __name__ == '__main__':
    main()
