'''
Created on Jan 24, 2012

@author: astar
'''
from pyparsing import Keyword
import re, traceback, sys, collections
###################################################################################################            
###################################################################################################   
class HmmerTabParser:
    def __init__(self, filePath,  **kwargs):
        self.__data = open(filePath)
        for key in kwargs:
            if key == 'e_threshold': self.evalue = kwargs[key]
            if key == 'overlap_limit': self.overlap = kwargs[key]
            if key == 'infoLine': self.infoLine = kwargs[key]
            if key == 'startLine': self.startLine = kwargs[key]  
        self.__parseData()
    ############################################################################           
    ############################################################################  
    def __assignAnnotation(self, read, candidate):
        counter = 0
        overlap = 0
        breakExit = False
        for el in self.hmmerAnnotation.get(read, "None"):
            area = el[1] - el[0]
            if area > candidate[1] - candidate[0]: area = candidate[1] - candidate[0]
            if candidate[0] >= el[0]: 
                if candidate[1] >= el[1]:
                    if candidate[0] <= el[1]:
                        overlap = (el[1] - candidate[0])/area
                else:
                    overlap = 1
            else:
                if candidate[1] <= el[1]:
                    if candidate[1] >= el[0]:
                        overlap = (candidate[1] - el[0])/area
                else:
                    overlap = 1   
            if overlap > self.overlap:  
                breakExit = True
                break
            counter += 1
        if not breakExit: counter -= 1
        if self.hmmerAnnotation[read][counter][3] > candidate[3]:
            if overlap > self.overlap: 
                self.hmmerAnnotation[read][counter] = candidate
                return
            else:
                self.hmmerAnnotation[read].append(candidate)
                return
        else:
            if overlap < self.overlap:
                self.hmmerAnnotation[read].append(candidate)
                return
            else:
                return
    ############################################################################           
    ############################################################################
    def __parseData(self):
        start = False
        self.hmmerAnnotation = collections.defaultdict(list)
        for line in self.__data:
            if not start:
                try:
                    self.startLine.parseString(line)
                    start = True
                    continue
                except:
                    continue
            if self.infoLine.search(line):
                anno = re.split("\s+", line.rstrip())
                try:
                    read = anno[0]
                    alignment_start = float(anno[1])
                    alignment_end = float(anno[2])
                    score = float(anno[11])
                    e_val = float(anno[12])
                    ko = anno[6]
                    if self.hmmerAnnotation.get(read, "None") == "None":
                        self.hmmerAnnotation[read].append((alignment_start, alignment_end, ko, e_val, score))
                    else:
                        if (e_val < self.evalue): self.__assignAnnotation(read, (alignment_start, alignment_end, ko, e_val, score))
                except:
                    print "### print_exception:"
                    traceback.print_exc(file=sys.stdout)
                    print anno
                    print "###"
                    
                    
                    
if __name__ == "__main__":
#    data = "/media/Transcend/workDir/PaulPrograms/src/adigi.json"
#    x = HmmerTabParser(data, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))                    
#    for key, item in x.hmmerAnnotation.items():
#        print "#######" + key + "#######"
#        print [el[2] for el in item]
#        print item
    
    from collections import defaultdict    
    data = "/media/Transcend/TAGS/Jura_meta/reads/927-2E_test.out"
    x = HmmerTabParser(data, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))       
    inv = defaultdict(list)
    annotation = 'K00004'
    for key, val in x.hmmerAnnotation.items():
        for el in val: 
            inv[el[2]].append(key)
    f = open("acc.txt", "w")
    reads = "\n".join([read.split("_")[0] for read in inv[annotation]])
    f.write(reads)
    f.close()

        
    
            


    