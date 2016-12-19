'''
Created on Sep 15, 2012

@author: astar
'''
from Bio.Blast import NCBIXML
class Blast(object):
    def __init__(self, result):
        f = open(result)
        self.records = NCBIXML.parse(f)
        self.data = []
        self.__parse()
    def __parse(self):
        for record in self.records:
            seq_id = record.query.split()[0]
            try:
                definition = record.alignments[0].hit_def.split(">")[0].rstrip()
            except:
                continue       
            coverage = 0
            identity = 0
            for hsp in record.alignments[0].hsps:
                coverage += hsp.align_length
                identity += hsp.identities
            percentage = float(identity)/float(coverage) * 100
            self.data.append((int(seq_id), definition, percentage))
            #print "%.2f" % percentage
#            print record.alignments[0].hsps[0].identities      
#            print record.alignments[0].length