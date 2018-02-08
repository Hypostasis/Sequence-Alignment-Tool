from tkinter import Tk, Label, Button, Entry, Checkbutton, IntVar
from algorithms import NeedlemanWunsch, align_star
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

class GUI:
    def __init__(self, master):
        self.master = master
        master.title("Sequence Alignment Tool")

        self.label = Label(master, text = "Enter your sequences here:")
        self.label.pack()
        self.field_S1 = Entry()
        self.field_S1.pack()
        self.field_S2 = Entry()
        self.field_S2.pack()

        self.var_biop = IntVar()
        self.biopython_check = Checkbutton(master, text = "Use Biopython implementation",
                                           variable = self.var_biop, onvalue=1, offvalue=0)
        self.biopython_check.pack()

        self.var_fromfile = IntVar()
        self.fromfile_check = Checkbutton(master, text = "Read sequences from file",
                                          variable = self.var_fromfile, onvalue=1, offvalue=0)
        self.fromfile_check.pack()

        self.align_button = Button(master, text="Align global", command=self.align_g)
        self.align_button.pack()
        self.align_button = Button(master, text="Align local", command=self.align_l)
        self.align_button.pack()
        self.align_button = Button(master, text="Dot matrix", command=self.dotMatrix)
        self.align_button.pack()
        self.align_button = Button(master, text="MSA", command=self.MSA)
        self.align_button.pack()

        self.close_button = Button(master, text="Close", command=master.quit)
        self.close_button.pack()
        self.sequences = []



    def align_g(self):

        if self.var_fromfile.get() == 1:
            pass
        else:
            seq1 = self.field_S1.get()
            seq2 = self.field_S2.get()

        if self.var_biop.get() == 1:
            alignments = pairwise2.align.globalxx(seq1, seq2)
            with open('output.txt', 'w') as file:
                file.write("SEQUENCE ALIGNMENT TOOL OUTPUT - GLOBAL ALIGNMENT\n")
                for i, alignment in enumerate(alignments):
                    file.write(format_alignment(*alignments[i]))
                    file.write("\n")
        else:
            NeedlemanWunsch(seq1, seq2)

    def align_l(self):
        if self.var_fromfile.get() == 1:
            pass
        else:
            seq1 = self.field_S1.get()
            seq2 = self.field_S2.get()
        alignments = pairwise2.align.localxx(seq1, seq2)
        with open('output.txt', 'w') as file:
            file.write("SEQUENCE ALIGNMENT TOOL OUTPUT - LOCAL ALIGNMENT\n")
            for i, alignment in enumerate(alignments):
                file.write(format_alignment(*alignments[i]))
                file.write("\n")

    def dotMatrix(self):
        return

    def load_fasta(self):
        """Loads sequences to be aligned from a fasta file"""

        self.sequences = []
        for seq_record in SeqIO.parse("input.fasta", "fasta"):
            self.sequences.append((seq_record.description, seq_record.seq))
        return

    def MSA(self):
        if self.var_fromfile.get() == 1:
            self.load_fasta()
            pass
        else:
            self.sequences = []
            seq1 = self.field_S1.get()
            seq2 = self.field_S2.get()
            self.sequences.append(("S1", seq1))
            self.sequences.append(("S2", seq2))
        align_star(self.sequences)

        return

root = Tk()
gui = GUI(root)
root.mainloop()
