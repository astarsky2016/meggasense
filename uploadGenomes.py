'''
Created on Dec 2, 2011

@author: astar
'''
import logging, sys
from myZODB import MyZODB, transaction
from genes.genomeGrammar import GenomeEntry
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('/media/Transcend/TAGS/myZODB/uploadGenomes.log')
db = MyZODB('/media/Transcend/TAGS/myZODB/keggGenomes.fs')
dbroot = db.dbroot
f = open('/media/Transcend/Working_KEGG/genes/genome', 'r')
testString = f.read().split('///')
for string in testString:
    try:
        x = GenomeEntry(string.lstrip().rstrip())
        dbroot[x.t_number] = x
    except Exception, err:
        log.exception('Error from GenomeEntry(): %s' % err)
        db.close()
        f.close()
        sys.exit(1)
transaction.commit()
db.close()
    

