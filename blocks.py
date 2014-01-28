#!/usr/lib/python2.7
# Created by F. Madeira, 2012 
# Biopython - http://biopython.org/wiki/Download
# EMBOSS - http://emboss.sourceforge.net/

import sys
import os
from itertools import combinations
from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC

def printUsage():
    """Prints the usage.
    """
    Usage = \
    """
    Routines to get blocks of conserved/gappy columns in
	Multiple Sequence Alignments (MSA), FASTA format to MSF
	conversion, and scoring with the Sum-of-pairs function.
    
    Usage:
    -c              Creates *.blocks with conserved columns
    -g              Creates *.blocks with gappy columns
    -msf            Converts any *.fa or *.fasta file to *.msf
                    MSF format (needs seqret from EMBOSS)
    -s              Scores the alignment (fasta) with the SP
                    (Sum-of-Pairs) scoring function
    -h or -help     Prints help 
    """
    print Usage        
        
def checkArguments():
    if command != 'c' and command != 'g' and \
       command != 's' and command != 'msf':
        printUsage()

def checkDependencies():
    try: 
        import Bio
        del Bio
    except ImportError:
        print "ERROR: Unable to import Biopython"
        sys.exit()
    
def main():
    args = sys.argv
    global command
    
    if len(args) != 2:
        if args[0] == "-h" or args[0] == "-help":
            printUsage()
        else: 
            printUsage()
        return
    else:
        arg1 = args[1].lstrip("-")
        command = arg1
        checkArguments()
        checkDependencies()
        if command == "c":
            conservedBlocks()
        elif command == "g":
            gappyBlocks()
        elif command == "s":
            alignScore()
        elif command == "msf":
            fasta2MSF()
        else: 
            return


def gappyBlocks():
    """Gappy blocks are defined as regions of consecutive columns,
    with length > 3, on which more than half of its positions are
    gaps (in accordance with the specified threshold).
    """
    path = "./data/"
    for file in os.listdir(path):
        if file.endswith(".fa") or file.endswith(".fasta"):
            alignin = AlignIO.read(path + file, "fasta")
            try:
                filecore = file.rstrip(".fa")
            except:
                filecore = file.rstrip(".fasta")
            fileout = path + filecore + ".blocks"
            
            # constants
            align = []
            gap = []
            border = []
            blocks = []
            
            # specify cut-off of gaps in column (in percentage)
            cut_min = 0.1
            cut_max = 0.9
            
            # alignment
            for pos in range(0,(alignin.get_alignment_length())):
                column=alignin[:,pos]
                align.append(column)
                if "-" in column:
                    col=list(column)
                    gaps=col.count("-")
                    if gaps > (cut_min*len(col)) and gaps < (cut_max*len(col)):
                        gap.append(pos)
                        
            if gap != []:
                border.append(gap[0])
                border.append(gap[len(gap)-1])
                for i in range(0,(len(gap)-1)):
                    if int(gap[i]+1)!=int(gap[i+1]):
                        border.append(gap[i])
                        
                for j in range((len(gap)-1), 0, -1):
                    if int(gap[j]-1)!=int(gap[j-1]):
                        border.append(gap[j])
                # list of positions for the blocks
                order=sorted(border)
                
                # get the blocks and writes to the .blocks file
                o=open(fileout, "w")
                
                for i in range(0,len(order)-1,2):
                    beg=int(order[i])
                    end=int(order[i+1])
                    count = end-beg    
                    block=alignin[:,beg:end]
                    
                    # specify the minimum length of a gap
                    if count < 3:
                        pass
                    else:        
                        blocks.append(block)                    
                        o.write('***Block***'+"\n"+"Start:"+str(beg)+\
                            "\n"+"Count:"+str(count)+"\n")
                        for record in block:
                            o.write(str(record.seq)+"\n")
                o.close()
            else:
                o=open(fileout, "w")
                o.close()
                pass
    return

def conservedBlocks():
    """Conserved blocks are defined as regions of consecutive columns,
    with length > 3, on which none of its positions are
    gaps (in accordance with the specified threshold).
    """
    path = "./data/"
    for file in os.listdir(path):
        if file.endswith(".fa") or file.endswith(".fasta"):
            alignin = AlignIO.read(path + file, "fasta")
            try:
                filecore = file.rstrip(".fa")
            except:
                filecore = file.rstrip(".fasta")
            fileout = path + filecore + ".blocks"
            
            # constants
            align = []
            cons = []
            border = []
            blocks = []
            
            # alignment
            for pos in range(0,(alignin.get_alignment_length())):
                column=alignin[:,pos]
                if "-" not in column:
                        align.append(column)
                        cons.append(pos)
                        
            
            if cons != []:           
                border.append(cons[0])
                border.append(cons[len(cons)-1])
                for i in range(0, len(cons)-1):
                    if int(cons[i]+1)!=int(cons[i+1]):
                        border.append(cons[i])
                        
                for j in range((len(cons)-1), 0, -1):
                    if int(cons[j]-1)!=int(cons[j-1]):
                        border.append(cons[j])    
                    
                # list of positions for the blocks
                order=sorted(border)
                
                # get the blocks and writes to the .blocks file
                o=open(fileout, "w")
                
                for i in range(0,len(order)-1,2):
                    beg=int(order[i])
                    end=int(order[i+1])
                    count = end-beg    
                    block=alignin[:,beg:end]
                    
                    # specify the minimum length of a gap
                    if count < 3:
                        pass
                    else:        
                        blocks.append(block)                    
                        o.write('***Block***'+"\n"+"Start:"+str(beg)+\
                            "\n"+"Count:"+str(count)+"\n")
                        for record in block:
                            o.write(str(record.seq)+"\n")
                o.close()
            else:
                o=open(fileout, "w")
                o.close()
                pass
    return

def fasta2MSF():
    """converts alignments from fasta format to msf format.
    """
    path = "./data/"
    for file in os.listdir(path):
        if file.endswith(".fa") or file.endswith(".fasta"):
            os.chdir(path)
            try:
                filecore = file.rstrip(".fa")
            except:
                filecore = file.rstrip(".fasta")
            fileout = filecore + ".msf2"
            
            seqret = os.system("seqret fasta::" + file + \
                            " msf::" + fileout)
            print seqret
            
            outmsf = filecore + ".msf"
            out = open(outmsf, "w")
            op = open(fileout, "r")
            msf = op.readlines()
            op.close()
            for line in msf:
                if line[0] == "\n":
                    print >> out, line.rstrip("\n")
                elif line[0] != "!" and line[0] != "/" and \
                                             line[0] != "\n":
                    line = line.replace(".", "-")
                    line = line.replace("~", "-")
                    print >> out, line.rstrip("\n")        
                else:
                    print >> out, line.rstrip("\n")
            out.close()
            
            # remove the comment if you want to remove the
            # original file
            #os.remove(file)
            
            os.remove(fileout)
            os.chdir("../")
    return

def alignScore():
    """Computes a score for the MSA inputed.
    Methods implemented:
    Sum-of-pairs (SP) score -  Murata et al, 1985
    as explained in Gonnet et al, 2000. 
    SP is the sum of all possible combinations of
    pairwise scores. 
    """
    matrix = mapMatrix("BLOSUM62")
    
    path = "./data/"
    for file in os.listdir(path):
        if file.endswith(".fa") or file.endswith(".fasta"):
            sequences = []
            input_sequences = SeqIO.parse(path + file, "fasta", \
                                    IUPAC.protein)

            for record in input_sequences:
                seq = str(record.seq)
                sequences.append(seq) 
    
    SumOfPairs = 0
    for pair in combinations(sequences, 2): 
        SumOfPairs += pairwiseScore(pair[0], pair[1], matrix)
    
    print SumOfPairs


def pairwiseScore(seq1, seq2, matrix):
    """
    Method Sum-of-pairs - as on Gonnet et al, 2000.
    s(x,y) = { matchScore(x,y) if x!='-' and y!='-';
               0 if x=='-' and y=='-', because the delection as caused earlier
               gap penalty, depending on the gap length  gap + length * increment}
               
    gap - depends on the scoring matrix (PAM, BLOSUM, etc)
    length - length of the gap
    increment - incremental penalty that depends on the scoring matrix
    
    BLOSUM62, gap = -4, increment = 1 -> increment = length 
    """
    
    gap = -4.0
    incr_top = 0
    incr_bottom = 0
    pairwise_score = 0
    for i,j in zip(range(len(seq1)), range(len(seq2))):
        aa1 = seq1[i]
        aa2 = seq2[j] 
        if aa1=="-" and aa2 =="-" :
            pairwise_score += 0
        elif aa1!="-" and aa2!="-":
            pairwise_score += float(matchScore(aa1, aa2, matrix))
        elif aa1=="-" and aa2!="-":
            try:
                aa11 = seq1[i+1]
                aa22 = seq2[j+1]
                if aa11=="-" and aa22!="-":
                    incr_top += 1
                else: 
                    pairwise_score += gap + incr_top * incr_top
                    incr_top = 0
            except: 
                pairwise_score += gap
                pass
        elif aa1!="-" and aa2=="-":
            try:
                aa11 = seq1[i+1]
                aa22 = seq2[j+1]
                if aa11!="-" and aa22=="-":
                    incr_bottom += 1
                else: 
                    pairwise_score += gap + incr_bottom * incr_bottom
                    incr_bottom = 0
            except: 
                pairwise_score += gap
                pass
        else: pass
        
    return pairwise_score
         
def matchScore(alpha, beta, matrix):
    "Matches scores from a matrix"
    
    alphabet = {}    
    alphabet["A"] = 0
    alphabet["R"] = 1
    alphabet["N"] = 2
    alphabet["D"] = 3
    alphabet["C"] = 4
    alphabet["Q"] = 5
    alphabet["E"] = 6
    alphabet["G"] = 7
    alphabet["H"] = 8
    alphabet["I"] = 9
    alphabet["L"] = 10
    alphabet["K"] = 11
    alphabet["M"] = 12
    alphabet["F"] = 13
    alphabet["P"] = 14
    alphabet["S"] = 15
    alphabet["T"] = 16
    alphabet["W"] = 17
    alphabet["Y"] = 18
    alphabet["V"] = 19
    alphabet["B"] = 20
    alphabet["Z"] = 21
    alphabet["X"] = 22
    alphabet["-"] = 22
    lut_x = alphabet[alpha]
    lut_y = alphabet[beta]
    
    return matrix[lut_x][lut_y]
    
def mapMatrix(matrix):
    "Maps a matrix of floats"
    matrix = matrix.upper()
    
    score_matrix = []
    input = matrix
    input_matrix = open(input, 'r')
    for line in input_matrix.readlines():
        score_matrix.append(map(float, line.split()))
    input_matrix.close()
    
    return score_matrix

if __name__ == "__main__":
    main()