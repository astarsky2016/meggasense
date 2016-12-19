'''
Created on Mar 26, 2012

@author: astar
'''
import csv, collections, re, sys
from nltk.tokenize import word_tokenize, sent_tokenize
from collections import defaultdict
class CAZyMapper:
    def __init__(self):
        cazyMapping = "/media/Transcend/TAGS/Jura_meta/dbCANhmm_mapping2_KEGGkoemit.csv"
        self.keggMap = defaultdict(list)
        reader = csv.reader(open(cazyMapping, 'rb'), delimiter=';')
        for row in reader:
            self.keggMap[row[1]].append(row[0])
        csvfile  = open('/media/Transcend/eclipseWorkspace/TAGS/src/CAZy/CAZY_families_descriptions_TABLE_cleaned_FINAL.csv', "rb")
        self.cazyDesc = collections.defaultdict(dict)
        reader = csv.reader(csvfile, delimiter=';')
        reader.next()
        for row in reader:
            try:
                self.cazyDesc[row[0]]['activities'] = self.__tokenize(row[1])
                self.cazyDesc[row[0]]['specificities'] = self.__tokenize(row[2])
            except:
                pass
        csvfile.close()
        
    def __tokenize(self, text):
        tokens = []
        for sent in sent_tokenize(text):
            tokens.extend([re.sub('^-', '', word) for word in word_tokenize(sent) if not re.match(r"[\x90-\xff]", word) and re.search("\w", word)])
        if tokens: 
            return " ".join(tokens)
        else:
            return None
if __name__ == "__main__":
    x = CAZyMapper()
    print x.cazyDesc['GH13']['specificities']
    for key, val in x.keggMap.items():
        print key, val
    print len(x.keggMap.items())

        
    