'''
Created on Mar 20, 2012

@author: astar
'''
import os, csv
from collections import defaultdict
from Bio import SeqIO
newblerDir = "/home/users/astar/TAGS/treangen-metAMOS-e602b91/newbler/bin"
data = '/home/users/astar/TAGS/glimmer-mg_02/scripts/HI7IROS01_02.hmmer_kegg58.out'
sff = '/home/users/astar/TAGS/assembly/HI7IROS01_02.sff'
from mongoStore import MongoConnect
from sffWrapper import AssemblyPipe
class Newbler:
    def __init__(self):
        AssemblyPipe(newblerDir, data, sff)
        for d in os.listdir(newblerDir):
            self.ko2contig2Reads = {}
            self._mapContigsToReads(d, os.path.join(os.path.abspath(d), "454ReadStatus.txt"))
            self.allcontigs = os.path.join(os.path.abspath(d), "454AllContigs.fna")
            
    def _mapContigsToReads(self, ko, allcontigmap):
        reader = csv.reader(open(allcontigmap, 'rb'), delimiter='\t')
        reader.next() 
        self.ko2contig2Reads[ko] = defaultdict(list)
        for row in reader:
            try:
                if len(row) > 2:
                    self.ko2contig2Reads[ko][row[2]].append(row[0])
                else:
                    self.ko2contig2Reads[ko][row[1]].append(row[0])
            except:
                print "####error in class Newbler Start:####"
                print row
                print "####error in class Newbler :End####" 
    def _store2Mongo(self):
        contig_dict = SeqIO.to_dict(SeqIO.parse(self.allcontigs, "fasta"))
        reads = SeqIO.index(sff, "sff")
        conn = MongoConnect("sample_A1_67_kegg")
        for ko, seqs in self.ko2contig2Reads.items():
            ko_id = ko
            for contig, reads in seqs.items():
                contig_id = contig
                if contig_dict.get(contig_id, None) == None:
                    fasta = None
                else:
                    fasta = contig_dict.get[contig_id].format("fasta")
                
                doc = {ko_id: {contig_id: {"fasta": fasta,
                                           "reads": [dict(acc, reads[acc].format("fasta")) for acc in reads]
                                           }                   
                              }     
                      }
                conn.dbh.insert
       
if __name__ == "__main__":
    Newbler()
        
         
    