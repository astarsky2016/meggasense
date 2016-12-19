'''
Created on Dec 6, 2011

@author: astar
'''
from myZODB import MyZODB
from mySqlite import MySqlite
import logging, sys, re, os, shutil
logging.basicConfig()
removeBrackets = re.compile(r'\(.*\)')
db = MyZODB('/media/Transcend/TAGS/myZODB/keggKOs.fs')
KO_dir = '/media/Transcend/TAGS/KO_seq'
KO_nuc = '/media/Transcend/TAGS/KO_seq/KO_nuc'
KO_pep = '/media/Transcend/TAGS/KO_seq/KO_pep'
#ko_empty = open('/media/Transcend/TAGS/KO_seq/empty_KO.txt', 'w')
dbroot = db.dbroot
count = 0
sqlite = MySqlite('/media/Transcend/TAGS/Sqlite3/fastaSequences.db')
ko_file = os.path.join(KO_dir, "K01980_CHECK.fasta")
f = open(ko_file, 'w')
peptides = False
nucleotides = False
for gene in dbroot["K01980"].genes:
    org = gene[0][0]
    for accs in gene[1]:
        gene = removeBrackets.sub('', accs)
        res = sqlite.getRecords('SELECT SeqPep, SeqNuc FROM keggFasta WHERE Id=? AND Organism=? COLLATE NOCASE', (gene, org))
        try:
            seqPep = res[0][0]
            seqNuc = res[0][1]
        except IndexError:
            print res
        if seqPep:
            if not nucleotides: f.write(seqPep + "\n")
            else: 
                print "DNA and protein sequences are not allowed to mix - error at: " + org + "-" + accs
                continue
            peptides = True
            nucleotides = False
            continue
        else:
            if seqNuc:
                if not peptides: f.write(seqNuc + "\n")
                else: 
                    print "DNA and protein sequences are not allowed to mix - error at: " + org + "-" + accs
                    continue
                nucleotides = True
                peptides = False
                continue
            else:
                print "bad entry: " + org + "-" + accs
                continue 
f.close() 
if peptides: shutil.move(ko_file, KO_pep)
else: shutil.move(ko_file, KO_nuc)
    
    
#for key in dbroot.keys():
#    ko = dbroot[key]
#    ko_file = os.path.join(KO_dir, key + ".fasta")
#    f = open(ko_file, 'w')
#    peptides = False
#    nucleotides = False
#    if ko.genes != 'undetermined':
#        for gene in ko.genes:
#            org = gene[0][0]
#            try:
#                for accs in gene[1]:
#                    gene = removeBrackets.sub('', accs)
#                    res = sqlite.getRecords('SELECT SeqPep, SeqNuc FROM keggFasta WHERE Id=? AND Organism=? COLLATE NOCASE', (gene, org))
#                    seqPep = res[0][0]
#                    seqNuc = res[0][1]
#                    if seqPep:
#                        if not nucleotides: f.write(seqPep + "\n")
#                        else: 
#                            print "DNA and protein sequences are not allowed to mix - error at: " + org + "-" + accs
#                            continue
#                        peptides = True
#                        nucleotides = False
#                        continue
#                    else:
#                        if seqNuc:
#                            if not peptides: f.write(seqNuc + "\n")
#                            else: 
#                                print "DNA and protein sequences are not allowed to mix - error at: " + org + "-" + accs
#                                continue
#                            nucleotides = True
#                            peptides = False
#                            continue
#                        else:
#                            print "bad entry: " + org + "-" + accs
#                            continue 
#            except IndexError, ex:
#                print key + " makes problems!"
#        else:
#            ko_empty.write(key + "\n")
#            count += 1
#    f.close() 
#    if peptides: shutil.move(ko_file, KO_pep)
#    else: shutil.move(ko_file, KO_nuc)
#ko_empty.close()                  
sqlite.closeHandle()
db.close()
print count