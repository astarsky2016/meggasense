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
path_in = '/media/Transcend/TAGS/Jura_meta/contigs'
path_out = '/home/astar/jura_metagenomes'
class SqliteLoader(object):
    def __init__(self, db, hmmerData, seqData, tables):
        self.db = db
        self.hmmerData = hmmerData
        self.seqData = seqData
        self.tables = tables
        self.__load()
    def __load(self):
        db = MySqlite(self.db)
        keggKoAnnotation = KeggData('ko') 
        for table in self.tables:
            db.createTable(table)
        for hmmer_file in self.hmmerData:
            hmmerFile = os.path.join(path_in, hmmer_file)
            for seq in self.seqData:
                if hmmer_file.split(".")[0] in seq.split('.'): 
                    print hmmer_file + " - " + seq
                    seqFile = os.path.join(path_in, seq)
                    sequences = SeqIO.to_dict(SeqIO.parse(seqFile, "fasta"))
                    x = HmmerTabParser(hmmerFile, seqFile, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))
                    records = (seq.split('.')[0], "contigs")
                    db.execute('INSERT INTO metagenomes (id, seq_type) VALUES (?, ?)', records)        
                    for id, item in x.hmmerAnnotation.items():
                        fasta = sequences[id].format("fasta")
                        records = (id, fasta, seq.split('.')[0], 'NULL')
                        db.execute('INSERT INTO sequences (id, fasta, metaid, assembly) VALUES (?, ?, ?, ?)', records)
                        min_Evalue = 100
                        annotation = 'NULL'
                        description = 0
                        for el in item:
                            name = keggKoAnnotation.dbroot[el[2]].aliases
                            definition = keggKoAnnotation.dbroot[el[2]].definition
                            class_element = 'NULL'
                            pathway = 'NULL'
                            pubId = 'NULL'
                            description = 0
                            if keggKoAnnotation.dbroot[el[2]].class_element != 'undetermined':
                                class_element = "\n".join(keggKoAnnotation.dbroot[el[2]].class_element)
                            if keggKoAnnotation.dbroot[el[2]].pathways != 'undetermined':
                                pathway = ''
                                for path in keggKoAnnotation.dbroot[el[2]].pathways:
                                    pathway += " ".join(path) + "\n"
                            if keggKoAnnotation.dbroot[el[2]].references != 'undetermined':
                                pubId = ''
                                for reference in keggKoAnnotation.dbroot[el[2]].references:
                                    pubId += str(reference[0]) + "\n"
                            records = (el[2], el[3], name, definition, class_element, pathway, pubId, el[5], description, id)
                            db.execute('INSERT INTO keggData (keggko, evalue, name, definition, class, pathways, pubids, orientation, description, seqid) VALUES (?,?,?,?,?,?,?,?,?,?)', records)
                            if el[3] < min_Evalue:
                                min_Evalue = el[3]
                                annotation = el[2]
                    print "UPDATE %i %s %s" % (1, annotation, id)
                    db.execute('UPDATE keggData SET description=? WHERE keggko=? AND seqid=?', (1, annotation, id))
        db.closeHandle() 
        keggKoAnnotation.close()  
if __name__ == "__main__":
    db = os.path.join(path_out, "metagenomesDB.sqlite")
    table1 = 'CREATE TABLE metagenomes (id VARCHAR(50) PRIMARY KEY, seq_type VARCHAR(10))'
    table2 = 'CREATE TABLE sequences (id VARCHAR(50), fasta TEXT, metaid VARCHAR(50), assembly VARCHAR(50), FOREIGN KEY(metaid) REFERENCES metagenomes(id), FOREIGN KEY(assembly) REFERENCES sequences(id))'
    table3 = 'CREATE TABLE keggData (keggko VARCHAR(8), evalue REAL, name TEXT, definition TEXT, class TEXT, pathways TEXT, pubids VARCHAR(20), \
              orientation VARCHAR(10), description INTEGER, seqid VARCHAR(50), FOREIGN KEY(seqid) REFERENCES sequences(id))'
    tables = [table1, table2, table3]
    hmmer_files = [f for f in os.listdir(path_in) if f.endswith('.out')]
    seq_files = [f for f in os.listdir(path_in) if f.endswith('.fna')]
    SqliteLoader(db, hmmer_files, seq_files, tables)
    print "finished"
            

