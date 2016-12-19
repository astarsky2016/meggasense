'''
Created on Dec 15, 2011

@author: astar
'''
import logging, sys, re
from storage.myZODB import MyZODB, transaction
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('/media/Transcend/TAGS/myZODB/testBrite.log')
db = MyZODB('/media/Transcend/TAGS/myZODB/keggBrite.fs')
dbroot = db.dbroot
for key in dbroot.keys():
    brite = dbroot[key]
    for el in brite.keys():
        if el == "Carbohydrate Metabolism":
            print el
            for cCat in brite[el]:
                print "\t" + cCat
                for d in brite[el][cCat]:
                    print "\t\t" + d
db.close()