'''
Created on Sep 6, 2012

@author: astar
'''
import os, re
from pyparsing import Keyword
from CAZy.cazyData import CAZyMapper
from csvWriter import UnicodeWriter
from HMMER.parseHmmer_tab import HmmerTabParser as readParser
from HMMER.parseHmmerScaffold_tab import HmmerTabParser as contigParser
from collections import defaultdict
cazy = CAZyMapper()
from sqliteManager import sqliteManager
conn = sqliteManager("/media/Transcend/TAGS/Jura_meta/reads/metagenomes.sqlite")
path_in = "/media/Transcend/TAGS/Jura_meta/MET5_6_7"
pairs = []
skip = []
hmmer_files = [f for f in os.listdir(path_in) if f.endswith('.out')]
stats = ["KO", "CAZy", "number of contig hits", "number of genes assembled", "number of reads assigned to KO"]
for hmm in hmmer_files:
    first = "".join(hmm.split("_")[0:-1])
    for next_hmm in hmmer_files:
        second = "".join(next_hmm.split("_")[0:-1])
        if (first == second) and hmm != next_hmm and hmm not in skip:
            pairs.append((os.path.join(path_in, hmm), os.path.join(path_in, next_hmm)))
            skip.extend([hmm, next_hmm])
def __mapIdtoKegg(data, dna):
    koMap = defaultdict(list)
    if dna: annotation = contigParser(data, dna, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))
    else: annotation = readParser(data, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))
    for met_id, item in annotation.hmmerAnnotation.items():   
        [koMap[el[2]].append(el[0:2] + (met_id,) + el[3:]) for el in item]
    return koMap
         
for pair in pairs:
    if pair[0].find("contig") != -1: 
        contigs = pair[0]
        reads = pair[1]
        print "contigs: ", contigs
    else: 
        contigs = pair[1]
        reads = pair[0]
    contigDNA = os.path.splitext(contigs)[0] + ".fna"
    statFile = os.path.splitext(contigs)[0] + ".stat"
    koMap_ctg = __mapIdtoKegg(contigs, contigDNA)
    koMap_rds = __mapIdtoKegg(reads, None)
    f = open(statFile, "w")
    writer = UnicodeWriter(f)
    writer.writerow(stats)
    filtered = []
    for ko, item in koMap_rds.items():
        filtered.append(ko)
        cazyModules = ""
        if cazy.keggMap.get(ko, None):
            cazyModules = "|".join(cazy.keggMap[ko])
        ctg_num = 0
        if koMap_ctg.get(ko, None):
            ctg_num = len(koMap_ctg[ko])
        read_num = len(item)
        library = os.path.basename(reads)
        assembled_genes = conn.getContigCount(ko, library.split(".")[0])
        row = [ko, cazyModules, ctg_num, assembled_genes, read_num]
        writer.writerow(row)
    for ko, item in koMap_ctg.items():
        if ko not in filtered:
            cazyModules = ""
            if cazy.keggMap.get(ko, None):
                cazyModules = "|".join(cazy.keggMap[ko])
            ctg_num = len(item)
            row = [ko, cazyModules, ctg_num, 0, 0]
            writer.writerow(row)
    f.close()
            
conn.close()
#
#    
#    
#            
#    
#            