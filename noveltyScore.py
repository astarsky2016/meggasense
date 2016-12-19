'''
Created on Sep 14, 2012

@author: astar
'''
from sqliteManager import sqliteManager
import csv, os, glob
from csvWriter import UnicodeWriter
from BLAST.blastXMLparser import Blast
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
class GOIExtractor(object):
    def __init__(self, genes):
        goi = genes
        out = "/media/Transcend/TAGS/Jura_meta/MET5_6_7"
        db = "/media/Transcend/TAGS/Jura_meta/reads/metagenomes.sqlite"
        self.data = defaultdict(lambda: defaultdict(list))
        geneNames = {}
        conn = sqliteManager(db)
        csvfile  = open(goi, "rb")
        reader = csv.reader(csvfile, delimiter=',')
        reader.next()
        for row in reader:
            geneNames[row[0]] = row[1].split(";")
        csvfile.close()
        SEQUENCES = defaultdict(list)
        counter = defaultdict(int)
        unique_seq = defaultdict(set)
        for row in conn.query("SELECT id FROM metagenomes"):
            partial = 0
            assembled = 0
            func_goi = 0
            cons_goi = 0
            contig_goi = 0
            GOI_sequences = {}
            metagenome = row[0]
            lib_key = "_".join(metagenome.split("_")[0:-1])
            if lib_key == "MET5_6_7":
                for gene, functions in geneNames.items():
                    kos = ", ".join(["'" + el + "'" for el in functions])
                    reads = [row for row in conn.query("SELECT  DISTINCT sequences.fasta, sequences.id FROM sequences, hits WHERE hits.keggid IN (%s) AND sequences.id = hits.seqid AND sequences.metaid = '%s' AND sequences.assembly_type = 'singleton'" % (kos, metagenome))]
                    assembled_reads = [str(row[0]) for row in conn.query("SELECT  DISTINCT sequences.id FROM sequences, hits WHERE hits.keggid IN (%s) AND sequences.id = hits.seqid AND sequences.metaid = '%s' AND sequences.assembly_type = 'assembled_read'" % (kos, metagenome))]
                    assembled_genes = [row for row in conn.query("SELECT DISTINCT sequences.fasta, sequences.id FROM sequences, assembly WHERE assembly.readid IN (%s) AND assembly.kegg IN (%s) AND assembly.contigid = sequences.id" % (", ".join(assembled_reads), kos))]
                    contig_genes = [row for row in conn.query("SELECT DISTINCT sequences.fasta, hits.start, hits.stop, sequences.id FROM sequences, hits WHERE hits.keggid IN (%s) AND hits.seqid = sequences.id AND assembly_type IS NULL AND sequences.metaid = '%s'" % (kos, metagenome))]
                    if contig_genes:
                        contig_goi += len(contig_genes)
                        for record in contig_genes:
                            counter[lib_key] += 1
                            seq = Seq("".join(record[0].split("\n")[1:]))[record[1]:record[2]]
                            desc = "%s | %i bp region extracted from contig" % (gene, len(seq))
                            unique_seq[lib_key].add(record[3])
                            SEQUENCES[lib_key].append(SeqRecord(seq, id = str(record[3]), description = desc))
                            self.data[lib_key][gene].append((record[3], record[1], record[2], "region extracted from assembled contig", "conservative approach"))
                    if reads:
                        partial += len(reads)
                        for record in reads:
                            counter[lib_key] += 1
                            seq = Seq("".join(record[0].split("\n")[1:]))
                            desc = "%s | singleton %i bp" % (gene, len(seq))
                            unique_seq[lib_key].add(record[1])
                            SEQUENCES[lib_key].append(SeqRecord(seq, id = str(record[1]), description = desc))      
                            self.data[lib_key][gene].append((record[1], -1, -1, "singleton", "functional approach"))      
                    if assembled_genes:
                        func_goi += len(assembled_genes)
                        for record in assembled_genes:
                            counter[lib_key] += 1
                            seq = Seq("".join(record[0].split("\n")[1:]))
                            desc = "%s | assembled gene %i bp" % (gene, len(seq))
                            unique_seq[lib_key].add(record[1])
                            SEQUENCES[lib_key].append(SeqRecord(seq, id = str(record[1]), description = desc))    
                            self.data[lib_key][gene].append((record[1], -1, -1, "assembled gene", "functional approach"))  
        for lib, seqs in SEQUENCES.items():
            print lib
            print "unique sekvencija ima %i" % len(unique_seq[lib])
            print "sekvencija ima: %i" % counter[lib]
            print "zapisa ima: %i" % len(SEQUENCES[lib])
            FASTA = os.path.join(out, lib + "_roi.fasta")      
            SeqIO.write(seqs, FASTA, "fasta")
            
        
        
        conn.close() 
class BlastForNovelty(object): 
    pass
#TODO - u mrcini_new imas runBlast, prenesi u ovu klasu
class TableOverview(object):
    def __init__(self, files, data):
        self.files = files
        self.data = data
        self.__createTable()
    def __createTable(self):
        for file in self.files:
            f = open(os.path.join("/media/Transcend/TAGS/Jura_meta/MET5_6_7", file + ".csv"), "w")
            writer = UnicodeWriter(f)
            topics = ["Type", "Number annotated"]
            writer.writerow(topics)
            rows = [(key, str(len(val))) for key, val in self.data[file].items()]
            writer.writerows(rows)
            f.close()


if __name__ == "__main__":
    genes = GOIExtractor("amylomics_genes_of_interest.csv")  
    rowHeaders = ["seq_id", "start", "stop", "assigned annotation", "most similar to annotated", "percent identity to most similar database hit", "status", "remarks" ]
    csvFiles = []
    for res in glob.glob("/media/Transcend/TAGS/Jura_meta/MET5_6_7/*.xml"):
        hits = sorted(Blast(res).data, key = lambda x: x[2])
        stats = res.split(".")[0] + ".csv"
        lib = "_".join(os.path.basename(res).split("_")[0:-1])
        csvFiles.append(lib)
        f = open(os.path.join("/media/Transcend/TAGS/Jura_meta/MET5_6_7", stats), "w")
        writer = UnicodeWriter(f)
        writer.writerow(rowHeaders)
        for hit in hits:
            for gene, items in genes.data[lib].items():
                for item in items:
                    if hit[0] == item[0]:
                        row = [item[0], item[1], item[2], gene, hit[1], "%.2f" % hit[2], item[3], item[4]]
                        writer.writerow(row)
        f.close()
    TableOverview(csvFiles, genes.data)
    
        
        
        
                
            
    
        
    
            
          
            
        
    
