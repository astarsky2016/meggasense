'''
Created on Sep 10, 2012

@author: astar
'''
from parseHmmerScaffold_tab import HmmerTabParser 
from sqliteManager import sqliteManager
from pyparsing import Keyword
import re, os, math
from Bio import SeqIO
path_in = '/media/Transcend/TAGS/Jura_meta/contigs'
path_out = '/media/Transcend/TAGS/Jura_meta/reads'
class SqliteUpdate(object):
    def __init__(self, db, hmmerData, seqData, columns):
        self.conn = sqliteManager(db)
        self.hmmerData = hmmerData
        self.seqData = seqData
        self.columns = columns
        self.seq_count = 0
        self.__update()
        self.__load()
        print "now backing up"
        self.conn.backupToDisc()
    def __update(self):
        self.conn.alterTable(self.columns)
    def __load(self):
        for hmmer_file in self.hmmerData:
            print hmmer_file
            hmmerFile = os.path.join(path_in, hmmer_file)
            for seq in self.seqData:
                if hmmer_file.split(".")[0] in seq.split('.'): 
                    print "updating %s" % hmmer_file
                    seqFile = os.path.join(path_in, seq)
                    sequences = SeqIO.to_dict(SeqIO.parse(seqFile, "fasta"))
                    x = HmmerTabParser(hmmerFile, seqFile, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))                    
                    metaid = seq.split('.')[0] + "_contigs"
                    for id, item in x.hmmerAnnotation.items():
                        for el in item:
                            updates = [int(el[0]), int(el[1]), el[2]]
                            self.conn.addColumnValues(id, metaid, updates)
if __name__ == "__main__":
    db = os.path.join(path_out, "metagenomes.sqlite")
    new_columns = ("start", "stop")
    hmmer_files = [f for f in os.listdir(path_in) if f.endswith('.out')]
    seq_files = [f for f in os.listdir(path_in) if f.endswith('.fna')]
    SqliteUpdate(db, hmmer_files, seq_files, new_columns)
    print "finished"
            

