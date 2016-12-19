'''
Created on Mar 20, 2012

@author: astar
'''
import os, csv
from collections import defaultdict
from Bio import SeqIO
from sffWrapper import AssemblyPipe
newblerDir = "/home/astar/BioApps/newbler/bin"
assemblyDir = '/home/astar/assembly'
dataDir = "/media/Transcend/TAGS/Jura_meta/MET5_6_7"
class Newbler:
    def __init__(self, hmmerData, sffData, assembly_tmp_dir):
        self.sff = sffData
        self.contig_dict = {}
        self.ko2contig2Reads = defaultdict(lambda: defaultdict(list))
        AssemblyPipe(newblerDir, hmmerData, sffData, assembly_tmp_dir)
        for d in os.listdir(assembly_tmp_dir):
            os.chdir(assembly_tmp_dir)
            if os.path.isdir(d):
                self.__mapContigsToReads(d, os.path.join(os.path.abspath(d), "454ReadStatus.txt"))
                allcontigs = os.path.join(os.path.abspath(d), "454AllContigs.fna")
                self.contig_dict[d] = SeqIO.to_dict(SeqIO.parse(allcontigs, "fasta"))
    def __mapContigsToReads(self, ko, allcontigmap):
        reader = csv.reader(open(allcontigmap, 'rb'), delimiter='\t')
        reader.next() 
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


       
if __name__ == "__main__":
    hmmer_files = [f for f in os.listdir(dataDir) if f.endswith('.out')]
    sff_files = [f for f in os.listdir(dataDir) if f.endswith('.sff')]
    for hmm in hmmer_files:
        hmmerFile = os.path.join(dataDir, hmm)
        for sff in sff_files:
            if hmm.split(".")[0] in sff.split("."): 
                sffFile = os.path.join(dataDir, sff)
                print hmm + " - " + sff
                assembly_tmp_dir = os.path.join(assemblyDir, hmm.split(".")[0])
                if not os.path.exists(assembly_tmp_dir):
                    os.makedirs(assembly_tmp_dir)
                Newbler(hmmerFile, sffFile, assembly_tmp_dir)
        
         
    