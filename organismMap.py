'''
Created on Dec 7, 2011

@author: astar
'''
import logging, sys
from myZODB import MyZODB, transaction
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('/media/Transcend/TAGS/myZODB/mapGenomes.log')
db = MyZODB('/media/Transcend/TAGS/myZODB/keggGenomes.fs')
dbroot = db.dbroot
map = {}
for key in dbroot.keys():
    genome = dbroot[key]
    print genome.aliases[0] + "-" + genome.aliases[1]
    break
db.close()
