'''
Created on Jul 3, 2012

@author: astar
'''
from csvWriter import UnicodeWriter
from pyparsing import Keyword
from HMMER.parseHmmerScaffold_tab import HmmerTabParser
import re
def __mapIdtoKegg():
    koMap = {}
    data = "/home/astar/jura_metagenomes/MET1_contigs.tab"
    x = HmmerTabParser(data, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))
    for cor_id, item in x.hmmerAnnotation.items():
        seqid = "|".join(cor_id.split("|")[0:2])
        print seqid
        koMap[seqid] = [el[2] for el in sorted(item, key=lambda x: x[3])]
    return koMap
def main():
    koMap = __mapIdtoKegg()
    print len(koMap)
    out = open('scatt_autoAnnootate.csv', 'w')
    writer = UnicodeWriter(out)
    writer.writerows([[key, ",".join(value)] for key, value in koMap.items()])
    out.close()
if __name__ == '__main__':
    main()      