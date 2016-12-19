'''
Created on Dec 2, 2011

@author: astar
'''
from Bio import SeqIO
from mySqlite import Database
import os, glob, redis
redis = redis.Redis(host='localhost', port=6379, db=0)
print redis.get("Aaci_0001")
#path = '/media/Transcend/Working_KEGG/seq_pep'
#for infile in glob.glob( os.path.join(path, '*.pep') ):
#    for seq_record in SeqIO.parse(infile, "fasta"):
#        redis.set(seq_record.id.split(":")[1], seq_record.format("fasta"))
