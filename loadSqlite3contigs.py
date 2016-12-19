'''
Created on Jun 19, 2012

@author: astar
'''
from HMMER.parseHmmerScaffold_tab import HmmerTabParser 
from keggDB.getData import KeggData
from storage.mySqlite import MySqlite
from pyparsing import Keyword
import re, os
from Bio import SeqIO
path_in = '/media/Transcend/TAGS/Jura_meta/MET5_6_7'
path_out = '/media/Transcend/TAGS/Jura_meta/reads'
class SqliteLoader(object):
    def __init__(self, db, hmmerData, seqData, tables):
        self.db = MySqlite(db)
        self.hmmerData = hmmerData
        self.seqData = seqData
        self.kegg = KeggData('ko') 
        self.tables = tables
        self.__load()
    def __loadKO(self, ko):
        name = self.kegg.dbroot[ko].aliases
        definition = self.kegg.dbroot[ko].definition
        class_element = None
        pathway = None
        pubId = None
        if self.kegg.dbroot[ko].class_element != 'undetermined':
            class_element = "\n".join(self.kegg.dbroot[ko].class_element)
        if self.kegg.dbroot[ko].pathways != 'undetermined':
            pathway = ''
            for path in self.kegg.dbroot[ko].pathways:
                pathway += " ".join(path) + "\n"
        if self.kegg.dbroot[ko].references != 'undetermined':
            pubId = ''
            for reference in self.kegg.dbroot[ko].references:
                pubId += str(reference[0]) + "\n"
        records = (ko, name, definition, class_element, pathway, pubId)
        self.db.execute('INSERT INTO keggData (keggko, name, definition, class, pathways, pubids) VALUES (?,?,?,?,?,?)', records)
    def __load(self):
        seq_count = self.db.getRecords('SELECT MAX(id) FROM sequences')[0][0]
        koSet = set([el[0] for el in self.db.getRecords("SELECT keggko FROM keggData")])
#        for table in self.tables:
#            db.createTable(table)
        for hmmer_file in self.hmmerData:
            hmmerFile = os.path.join(path_in, hmmer_file)
            for seq in self.seqData:
                if hmmer_file.split(".")[0] in seq.split('.'): 
                    seqFile = os.path.join(path_in, seq)
                    sequences = SeqIO.to_dict(SeqIO.parse(seqFile, "fasta"))
                    x = HmmerTabParser(hmmerFile, seqFile, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))                    
                    metaid = seq.split('.')[0] + "_contigs"
                    records = (metaid, "contigs")
                    self.db.execute('INSERT INTO metagenomes (id, seq_type) VALUES (?, ?)', records)        
                    for id, item in x.hmmerAnnotation.items():
                        seq_count += 1
                        fasta = sequences[id].format("fasta")
                        records = (seq_count, id, fasta, metaid, None)
                        self.db.execute('INSERT INTO sequences (id, seqid, fasta, metaid, assembly_type) VALUES (?, ?, ?, ?, ?)', records)
                        min_Evalue = 100
                        annotation = None
                        description = 0
                        for el in item:
                            if el[2] not in koSet:
                                koSet.add(el[2])
                                self.__loadKO(el[2])
                            records = (el[2], seq_count, el[3], el[5], description, int(el[0]), int(el[1]))
                            self.db.execute('INSERT INTO hits (keggid, seqid, evalue, orientation, description, start, stop) VALUES (?,?,?,?,?,?,?)', records)
                            if el[3] < min_Evalue:
                                min_Evalue = el[3]
                                annotation = el[2]
                        self.db.execute('UPDATE hits SET description=? WHERE keggid=? AND seqid=?', (1, annotation, seq_count))
        self.db.closeHandle() 
        self.kegg.close()  
if __name__ == "__main__":
    db = os.path.join(path_out, "metagenomes.sqlite")
    table1 = 'CREATE TABLE metagenomes (id VARCHAR(50) PRIMARY KEY, seq_type VARCHAR(10))'
    table2 = 'CREATE TABLE sequences (id INTEGER PRIMARY KEY, seqid VARCHAR(50), fasta TEXT, metaid VARCHAR(50), assembly_type VARCHAR(50), FOREIGN KEY(metaid) REFERENCES metagenomes(id))'
    table3 = 'CREATE TABLE keggData (keggko VARCHAR(8) PRIMARY KEY, name TEXT, definition TEXT, class TEXT, pathways TEXT, pubids VARCHAR(20))'
    table4 = 'CREATE TABLE hits (keggid VARCHAR(8), seqid INTEGER, evalue REAL, orientation VARCHAR(10), description INTEGER, start INTEGER, stop INTEGER\
              FOREIGN KEY(seqid) REFERENCES sequences(id), FOREIGN KEY(keggid) REFERENCES keggData(id))'
    table5 = 'CREATE TABLE assembly (readid INTEGER, contigid INTEGER, kegg VARCHAR(8), FOREIGN KEY(readid) REFERENCES sequences(id), FOREIGN KEY(contigid) REFERENCES sequences(id))'
    tables = [table1, table2, table3, table4, table5]
    hmmer_files = [f for f in os.listdir(path_in) if f.endswith('.out')]
    seq_files = [f for f in os.listdir(path_in) if f.endswith('.fna')]
    SqliteLoader(db, hmmer_files, seq_files, tables)
    print "finished"
            

