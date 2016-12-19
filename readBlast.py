'''
Created on Jan 11, 2012

@author: astar
'''
from Bio.Blast import NCBIXML
import sys
result_handle = open("/home/astar/Desktop/metaGeneMark.blast.xml")
blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        print alignment.title
        for hsp in alignment.hsps:
            print hsp.query
            print hsp.num_alignments
            result_handle.close()
            sys.exit(1)