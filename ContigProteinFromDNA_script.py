import os

#for converting RNA to protein
translatetable = {
        "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
        "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
        "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
        "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

#for converting DNA to protein
translatetable2 = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                   'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                   'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
                   'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
                   'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                   'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                   'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                   'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                   'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                   'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                   'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                   'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                   'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                   'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                   'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                   'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

###################################
# FUTURE USE OF THIS PROGRAM JUST NEEDS TO REDIRECT THE FILE PATH
# THIS PROGRAM ITERATES THROUGH A FOLDER CONTAINING RNA SEQUENCES WHICH ARE SEPARATED AND CONTAINS A HEADER
###################################

def TranslateSeq():
    path = '/filepath containing CDS' #shifted path to all the sequences found 
    folder = os.fsencode(path)
    contigprotein = open("filepath to save the text file", "a+") #shifted cspcontigprotein to desktop
    for seqs in os.listdir(folder):
        seq = os.fsdecode(seqs)
        #seq is the str format of the file name in the current iteration
        if seq == '.DS_Store':
            continue
        else:
            protein = "" #add the growing protein sequence here
            protpath = 'output file path'
            proteinoutfilename = protpath
            proteinoutfilename += "Protein"
            proteinoutfilename += seq
            proteinout = open((proteinoutfilename), "a+") #creates a new txt file for the output protein
            opener = path + "/" + seq #directory for the seq file in question
            open_seq = open(opener, "r") #opens the seq file 
            read_seq = open_seq.read() #reads the seq file which has been opened
            if "+" in read_seq:
                startHere = read_seq.index("+") + 1
            elif "-" in read_seq:
                startHere = read_seq.index("-") + 1
            DNAonly = read_seq[startHere:].replace("\n","") #removes the line break that separates sequences
            for i in range(0, len(DNAonly), 3):
                try:
                    protein += translatetable2[DNAonly[i:i+3]]
                except Exception: #in case of errors
                    continue
            print(seq)
            print(protein)
            proteinout.write(read_seq[:startHere] + "\n")
            proteinout.write(protein) #correct
            contigprotein.write(read_seq[:startHere])
            contigprotein.write("\n")
            contigprotein.write(protein)
            contigprotein.write("\n")

