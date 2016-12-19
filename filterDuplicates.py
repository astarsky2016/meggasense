'''
Created on Oct 25, 2012

@author: astar
'''
import subprocess, StringIO, sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
sequences = "/media/Transcend/TAGS/Jura_meta/stats/A1-67_full_plate_roi.fa"
blast_cmd = "/home/astar/BioApps/ncbi-blast-2.2.26+/bin/blastn -query query.fa -subject subject.fa -outfmt 5"
f = open("/media/Transcend/TAGS/Jura_meta/stats/filtered.out", "w")
threshold = 0.95
count = 0
mem = dict([entry for entry in enumerate(SeqIO.parse(sequences, "fasta"))])
filtered = set(mem.keys())
print len(filtered)
print filtered
for query in mem.keys():
    rest = set(mem.keys())
    rest.remove(query)
    for subject in rest:
        query_tmp = open("query.fa", "w").write((mem[query].format("fasta")))
        subject_tmp = open("subject.fa", "w").write((mem[subject].format("fasta")))
        p = subprocess.Popen(blast_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
        tmp_xml = StringIO.StringIO(p)
        comparison = NCBIXML.read(tmp_xml)
        if len(comparison.alignments) > 0:
            for aln in comparison.alignments:
                for hsp in aln.hsps:
                    ratio = 0
                    duplicate = None
                    query_len = float(comparison.query_length)
                    subject_len = float(aln.length)
                    identities = float(hsp.identities)
                    alignment_length = float(hsp.align_length)
                    if query_len < subject_len: 
                        ratio = identities/query_len
                        duplicate = query
                    else: 
                        ratio = identities/subject_len
                        duplicate = subject
                    if ratio > threshold:
                        if identities/alignment_length >= ratio:
                            try:
                                filtered.remove(duplicate)
                            except KeyError:
                                print "already removed!"
                                continue
print len(filtered)
print filtered
f.write("\n".join(set([mem[el].id for el in filtered])))
FASTA = "/media/Transcend/TAGS/Jura_meta/stats/A1-67_full_plate_filtered_roi.fasta"
SeqIO.write([mem[el] for el in filtered], FASTA, "fasta")
f.close()
             



        
    