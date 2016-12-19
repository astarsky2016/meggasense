'''
Created on Sep 13, 2012

@author: astar
'''
from Bio import SeqIO, SeqUtils
from sqliteManager import sqliteManager
import os, sys
path = "/media/Transcend/TAGS/Jura_meta/reads"
out = "/media/Transcend/TAGS/Jura_meta/stats"
db = "/media/Transcend/TAGS/Jura_meta/reads/metagenomes.sqlite"
conn = sqliteManager(db)
sff_files = [f for f in os.listdir(path) if f.endswith('.sff')]
for sff in sff_files:
    print sff
    average_untrimmed = 0
    average_trimmed = 0
    input = os.path.join(path, sff)
    output = sff.split(".")[0] + ".jebek"
    lib = sff.split(".")[0] + "_reads"
    outFile = open(os.path.join(out, output), "w")
    count = 0
    query_res = conn.query("SELECT COUNT (DISTINCT sequences.id) FROM sequences, hits WHERE sequences.assembly_type IS NOT 'contig' AND sequences.metaid = '%s' AND sequences.id = hits.seqid" % lib)
    for row in query_res:
        ko_reads = row[0]
    for record in SeqIO.parse(input, 'sff'):
        count += 1
        average_untrimmed += len(record)
        average_trimmed += len([el for el in record.seq if el.isupper()])
    average_trimmed = round(float(average_trimmed)/float(count), 2)
    average_untrimmed = round(float(average_untrimmed)/float(count), 2)    
    outFile.write("total number of reads in library:\t%i\n" % count)
    outFile.write("total number of functional annotation reads in library:\t%i\n" % ko_reads) 
    outFile.write("average lenght of untrimmed reads in library: %.2f\n" % average_untrimmed)
    outFile.write("average lenght of trimmed reads in library: %.2f\n" % average_trimmed)
    outFile.close()

conn.close()
        
        
    