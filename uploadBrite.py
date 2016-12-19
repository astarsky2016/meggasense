'''
Created on Dec 14, 2011

@author: astar
'''
import logging, sys, re
from myZODB import MyZODB, transaction
from brite.koBriteGrammar import koBriteEntry
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('/media/Transcend/TAGS/myZODB/uploadBrite.log')
db = MyZODB('/media/Transcend/TAGS/myZODB/keggBrite.fs')
dbroot = db.dbroot
f = open('/media/Transcend/Working_KEGG/KEGG_core/ko00001.keg', 'r')
testString = f.read().split('#\n')
for string in testString:
    x = koBriteEntry(string)
    catA = x.A
    brite = dict()
    for i in range(len(x.rest)): 
        if i % 2 == 0:
            catB = x.rest[i]
            brite[catB] = dict()
        else:
            rest = x.rest[i]
            for i in range(len(rest)): 
                if i % 2 == 0:
                    listC = rest[i]
                    for el in listC:
                        catC = el
                        brite[catB][catC] = []
                else:
                    catD = rest[i]
                    brite[catB][catC] = catD
    dbroot[catA] = brite
    transaction.commit()
db.close()

