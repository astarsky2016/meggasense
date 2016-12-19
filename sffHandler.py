'''
Created on Mar 20, 2012

@author: astar
'''
from Bio import SeqIO
import sys
fasta = open("test.fasta", "w")
sffFile = open("test.sff", "wb")
#reads = SeqIO.to_dict(SeqIO.parse("GIQLXTN02.sff", "sff")) brzo rjesenje
reads = SeqIO.index("GIQLXTN02.sff", "sff")
ids = ["GIQLXTN02EL3SK", "GIQLXTN02EGDK6", "GIQLXTN02ENE6X"]
def generateRecords(params):
    for id in params:
        yield reads[id]
for read in generateRecords(ids):
    print read.format("fasta")
    print read.format("qual")
    print "#############"
sff = SeqIO.SffIO.SffWriter(sffFile)    
sff.write_file(generateRecords(ids))  
for rec in SeqIO.parse("test.sff", "sff"):
    print rec.id, len(rec), rec.seq

