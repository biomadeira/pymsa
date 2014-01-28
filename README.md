## Some Python routines to work with Multiple Sequence Alignments (MSA)
Routines to get blocks of conserved/gappy columns in Multiple Sequence Alignments (MSA). 
Converts fasta format to msf (needs seqret). 
Scores an MSA with the Sum-of-pairs scoring function.

###Dependencies

[Biopython 1.58](http://biopython.org/)

Optional:
[seqret](http://emboss.sourceforge.net/) for converting fasta->msf

###Usage
	-c              Creates *.blocks with conserved columns
	-g              Creates *.blocks with gappy columns
	-msf            Converts any *.fa or *.fasta file to *.msf
	                MSF format (needs seqret from EMBOSS)
    -s              Scores the alignment (fasta) with the SP
	                (Sum-of-Pairs) scoring function	
	-h or -help     Prints help 

Make sure you insert input files in the ./data folder.
Examples in ./data/ folder.

*F. Madeira, 2012*

