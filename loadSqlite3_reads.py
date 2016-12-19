'''
Created on Jun 19, 2012

@author: astar
'''
from HMMER.parseHmmer_tab import HmmerTabParser 
from keggDB.getData import KeggData
from storage.mySqlite import MySqlite
from pyparsing import Keyword
import re, os
from Bio import SeqIO
from assembler import Newbler
from collections import defaultdict
db_path = '/media/Transcend/TAGS/Jura_meta/reads'
path = "/media/Transcend/TAGS/Jura_meta/MET5_6_7"
assemblyDir = '/home/astar/assembly'
class SqliteLoader(object):
    def __init__(self, db, hmmerData):
        self.sqlitedb = MySqlite(db)
        self.hmmerData = hmmerData
        self.kegg = KeggData('ko') 
        self.koSet = set([el[0] for el in self.sqlitedb.getRecords("SELECT keggko FROM keggData")])
        self.seq_count = self.sqlitedb.getRecords('SELECT MAX(id) FROM sequences')[0][0]
        self.__assemble()
    def __assemble(self):
        for hmm in self.hmmerData:
            annotationData = defaultdict(list)
            hmmerFile = os.path.join(path, hmm)
            print hmm 
            sffFile = os.path.join(path, hmm.split(".")[0] + ".sff")
            assembly_tmp_dir = os.path.join(assemblyDir, hmm.split(".")[0])
            annotation = HmmerTabParser(hmmerFile, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))                  
            [annotationData[key.split("_")[0]].extend(val) for key, val in annotation.hmmerAnnotation.items()]
            if not os.path.exists(assembly_tmp_dir):
                os.makedirs(assembly_tmp_dir)
                self.assembly = Newbler(annotationData, sffFile, assembly_tmp_dir)
                seqFile = os.path.join(path, hmm.split(".")[0] + ".fasta")
                self.__load(seqFile, annotationData)
            else:
                self.assembly = Newbler(annotationData, sffFile, assembly_tmp_dir)
                #continue
                seqFile = os.path.join(path, hmm.split(".")[0] + ".fasta")
                self.__load(seqFile, annotationData)
        self.sqlitedb.closeHandle() 
        self.kegg.close() 
    def __loadHits(self, read, annotationData):
        annotated = annotationData.get(read, None)  
        min_Evalue = 100
        anno = None
        if annotated != None: 
            for candidate in annotated:
                records = (candidate[2], self.seq_count, candidate[3], None, 0, None, None)
                self.sqlitedb.execute('INSERT INTO hits (keggid, seqid, evalue, orientation, description, start, stop) VALUES (?,?,?,?,?,?,?)', records)
                if candidate[3] < min_Evalue:
                    min_Evalue = candidate[3]
                    anno = candidate[2]
            self.sqlitedb.execute('UPDATE hits SET description=? WHERE keggid=? AND seqid=?', (1, anno, self.seq_count))

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
        self.sqlitedb.execute('INSERT INTO keggData (keggko, name, definition, class, pathways, pubids) VALUES (?,?,?,?,?,?)', records)
        self.koSet.add(ko)
    def __load(self, seqFile, annotationData):        
        assembled_reads = {}                      
        sequences = SeqIO.to_dict(SeqIO.parse(seqFile, "fasta"))
        metaid = os.path.basename(seqFile).split('.')[0] + "_reads"
        records = (metaid, "reads")
        self.sqlitedb.execute('INSERT INTO metagenomes (id, seq_type) VALUES (?, ?)', records)    
        for ko, seqs in self.assembly.ko2contig2Reads.items():
            if ko not in self.koSet:
                self.__loadKO(ko)
            for contId, reads in seqs.items():
                if self.assembly.contig_dict[ko].get(contId, None) != None:
                    self.seq_count += 1
                    assembly_reference = self.seq_count
                    records = (self.seq_count, contId, self.assembly.contig_dict[ko][contId].format("fasta"), metaid, "contig")
                    self.sqlitedb.execute('INSERT INTO sequences (id, seqid, fasta, metaid, assembly_type) VALUES (?, ?, ?, ?, ?)', records)
                    for read in reads:
                        if read not in assembled_reads.keys():
                            self.seq_count += 1
                            assembled_reads[read] = self.seq_count
                            records = (self.seq_count, read, sequences[read].format("fasta"), metaid, "assembled_read")
                            self.sqlitedb.execute('INSERT INTO sequences (id, seqid, fasta, metaid, assembly_type) VALUES (?, ?, ?, ?, ?)', records)
                            records = (self.seq_count, assembly_reference, ko)
                            self.sqlitedb.execute('INSERT INTO assembly (readid, contigid, kegg) VALUES (?, ?, ?)', records)
                            self.__loadHits(read, annotationData)
                        else:
                            records = (assembled_reads[read], assembly_reference, ko)
                            self.sqlitedb.execute('INSERT INTO assembly (readid, contigid, kegg) VALUES (?, ?, ?)', records)
        for singleton in set(sequences.keys()).difference(assembled_reads):
            self.seq_count += 1
            records = (self.seq_count, singleton, sequences[singleton].format("fasta"), metaid, "singleton")
            self.sqlitedb.execute('INSERT INTO sequences (id, seqid, fasta, metaid, assembly_type) VALUES (?, ?, ?, ?, ?)', records)
            self.__loadHits(singleton, annotationData)
                        
                        
if __name__ == "__main__":
    db = os.path.join(db_path, "metagenomes.sqlite")
    hmmer_files = [f for f in os.listdir(path) if f.endswith('.out')]
    SqliteLoader(db, hmmer_files)
    print "finished"
            

