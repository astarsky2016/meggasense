'''
Created on Dec 6, 2011

@author: astar
'''
from Bio import SeqIO
from mySqlite import MySqlite
import os, glob, sys
pathPep = '/home/astar/BioApps/TAGS/seq_pep'
pathNuc = '/home/astar/BioApps/TAGS/seq_nuc'
db = MySqlite("/home/astar/BioApps/TAGS/fastaSequences.db")
table = 'CREATE TABLE keggFasta (Id VARCHAR(50) not null, Organism VARCHAR(50) not null, SeqPep TEXT, SeqNuc TEXT, PRIMARY KEY (Id, Organism))'
db.createTable(table)
for infile in glob.glob(os.path.join(pathPep, '*.pep')):
    for seq_record in SeqIO.parse(infile, "fasta"):
        org = seq_record.id.split(":")[0]
        geneId = seq_record.id.split(":")[1] 
        seqPep = seq_record.format("fasta")
        insert = 'INSERT INTO keggFasta (Id, Organism, SeqPep) VALUES (?, ?, ?)'
        records = (geneId, org, seqPep)
        db.execute(insert, records)
    db.commit()
print "Inserted all peptides!"
for infile in glob.glob(os.path.join(pathNuc, '*.nuc')):
    for seq_record in SeqIO.parse(infile, "fasta"):
        org = seq_record.id.split(":")[0]
        geneId = seq_record.id.split(":")[1] 
        seqNuc = seq_record.format("fasta")
        if db.getRecords('SELECT * FROM keggFasta WHERE Id=? AND Organism=?', (geneId,org)):
            update = 'UPDATE keggFasta SET SeqNuc=? WHERE Id=? AND Organism=?'
            records = (seqNuc, geneId, org)
            db.execute(update, records)
        else:
            insert = 'INSERT INTO keggFasta (Id, Organism, SeqNuc) VALUES (?, ?, ?)'
            records = (geneId, org, seqNuc)
            db.execute(insert, records)
    db.commit()
print "Inserted all nucleotides!"
db.closeHandle()

    
