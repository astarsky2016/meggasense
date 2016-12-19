'''
Created on Mar 27, 2012

@author: astar
'''
import os, csv, optparse, re
from collections import defaultdict
from Bio import SeqIO
newblerDir = "/home/users/astar/TAGS/treangen-metAMOS-e602b91/newbler/bin"
data = '/home/users/astar/TAGS/glimmer-mg_02/scripts/HI7IROS01_02.hmmer_kegg58.out'
sff = '/home/users/astar/TAGS/assembly/HI7IROS01_02.sff'
assembly_dir = '/home/users/astar/TAGS/assembly/sample_A1_67_kegg_assembly'
from mongoStore import MongoConnect
def mapContigsToReads(ko, allcontigmap):
        reader = csv.reader(open(allcontigmap, 'rb'), delimiter='\t')
        reader.next() 
        ko2contig2Reads[ko] = defaultdict(list)
        for row in reader:
            try:
                if len(row) > 2:
                    ko2contig2Reads[ko][row[2]].append(row[0])
                else:
                    ko2contig2Reads[ko][row[1]].append(row[0])
            except:
                print "####error in class Newbler Start:####"
                print row
                print "####error in class Newbler :End####" 
def store2Mongo(allcontigs):
    contig_dict = SeqIO.to_dict(SeqIO.parse(allcontigs, "fasta"))
    reads = SeqIO.index(sff, "sff")
    conn = MongoConnect("annotations")
    for ko, seqs in ko2contig2Reads.items():
        ko_id = ko
        contigData = []
        for contig, readAcc in seqs.items():
            contigDoc = {}
            if contig_dict.get(contig, None) == None:
                fasta = None
            else:
                fasta = contig_dict[contig].format("fasta")
            contigDoc["assembly_fasta"] = fasta
            fastaId = ['read_fasta']*len(readAcc)
            contigDoc["reads"] = [dict([pair]) for pair in zip(fastaId, [reads[acc].format("fasta") for acc in readAcc])]
            contigData.append(contigDoc)           
        doc = {"ko_id": ko_id,
               "assembly": contigData
              }
        conn.dbh.fasta_base.insert(doc, safe=True)
        print "inserted"
parser = optparse.OptionParser()
parser.add_option('-i', dest='assembly_folders', help='folders with assembly data')
(options,args) = parser.parse_args()
dataDirs = options.assembly_folders
ko2contig2Reads = {}
f = open(dataDirs)
for d in f:
    if re.search('\w', d):
        dir = d.rstrip()
        ko = os.path.basename(dir)
        mapContigsToReads(ko, os.path.join(dir, "454ReadStatus.txt"))
        store2Mongo(os.path.join(dir, "454AllContigs.fna"))
f.close()