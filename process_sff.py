'''
Created on Oct 5, 2012

@author: astar
'''
from Bio import SeqIO
sizes = [len(rec) for rec in SeqIO.parse("/media/Transcend/TAGS/Jura_meta/MET5_6_7/MET5_6_7.sff", "sff-trim")]
import pylab
pylab.hist(sizes, bins=100)
pylab.title("%i Library reads\nlengths: %i to %i" % (len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Read Lengths (bp, trimmed)")
pylab.ylabel("Number of Library Reads")
pylab.show()