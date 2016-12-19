'''
Created on Dec 2, 2011

@author: astar
'''
import logging, sys, re
from myZODB import MyZODB, transaction
from genes.koGrammar import KoEntry
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('/media/Transcend/TAGS/myZODB/uploadKOs.log')
db = MyZODB('/media/Transcend/TAGS/myZODB/keggKOs.fs')
dbroot = db.dbroot
f = open('/media/Transcend/Working_KEGG/genes/ko', 'r')
testString = f.read().split('///')
for string in testString:
    if re.search(r'\w', string):
        try:
            x = KoEntry(string.lstrip().rstrip() + "///")
            dbroot[x.ko_number] = x
            transaction.commit()
        except Exception, err:
            log.exception('Error from KoEntry(): %s' % err)
            db.close()
            f.close()
            sys.exit(1)
transaction.commit()
db.close()