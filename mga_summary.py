'''
Created on Oct 29, 2012

@author: astar
'''
from Bio import SeqIO
from pyparsing import *
from csvWriter import UnicodeWriter
import sys
sequences = "/media/Transcend/TAGS/Jura_meta/MET5_6_7/MET5_6_7_roi.fasta"
f = open("/media/Transcend/TAGS/Jura_meta/stats/A1-67_full_plate_roi.MetaGeneMark.txt")
mga = f.read()
f.close()
#############
start = LineStart().suppress() + LineEnd().suppress()
block = start + SkipTo(start)
##############
header = OneOrMore(Combine(LineStart().suppress() + Word(alphanums) + SkipTo(LineEnd().suppress())))("header")
##############
prot_start = LineStart() + Literal("##Protein").suppress() + Word(nums) + LineEnd().suppress()
prot_end = SkipTo("##end-Protein")
protein = OneOrMore(prot_start + prot_end)("protein")
#############
dna_start = LineStart() + Literal("##DNA").suppress() + Word(nums) + LineEnd().suppress()
dna_end = SkipTo("##end-DNA")
dna = OneOrMore(dna_start + dna_end)("dna")
######################
#matches = block.searchString(mga)
#for match in matches:
#    headers = header.searchString(match[0])
#    proteins = protein.searchString(match[0])
#    dnas = dna.searchString(match[0])
#    if headers:
#        genes = [element.header[0].split("GeneMark.hmm")[0].rstrip() for element in headers]
#    if protein:
#        for protein in proteins:
#            print protein.protein
            #sys.exit(1)
count = 0
for seq in SeqIO.parse(sequences, "fasta"):
    count += 1
print count