'''
Created on Dec 2, 2011

@author: astar
'''
import logging, sys, re, os
from myZODB import MyZODB, transaction
from genes.genesGrammar import GeneEntry
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('/home/astar/BioApps/TAGS/uploadGenes.log')
db = MyZODB('/home/astar/BioApps/TAGS/keggGenes.fs')
dbroot = db.dbroot
count = 0
monitor = 100
for line in open('/home/astar/BioApps/TAGS/genes/organisms.lst', 'r'):
    if re.search(r'\w', line):
        organism = os.path.join("/home/astar/BioApps/TAGS/genes", line.lstrip().rstrip())
        count += 1
        try:
            testString = open(organism, 'r').read().split('///')
        except IOError:
            print "##### no such file: " + organism + " #####"
            continue
        for string in testString:
            if re.search(r'\w', string):
                try:
                    x = GeneEntry(string.lstrip().rstrip() + "///")
                    dbroot[x.entry[0]] = x
                except Exception, err:
                    log.exception('Error from GeneEntry(): %s' % err)
                    print organism
                    print string
                    print "Unknown grammar at: " +  x.entry[0]
                    continue
        transaction.commit()
        if count == monitor: 
            print "finished loading %i genomes" % count
            monitor += count
db.close()
