'''
Created on Aug 8, 2012

@author: astar
'''
from pyparsing import Keyword
import re, traceback, sys, collections
from Bio import SeqIO
###################################################################################################            
###################################################################################################   
class HmmerTabParser:
    def __init__(self, filePath, contigDna, **kwargs):
        self.__data = open(filePath)
        self.__contigLengths = dict([(seq.id, len(seq)) for seq in SeqIO.parse(open(contigDna), "fasta")]) 
        for key in kwargs:
            if key == 'e_threshold': self.evalue = kwargs[key]
            if key == 'overlap_limit': self.overlap = kwargs[key]
            if key == 'infoLine': self.infoLine = kwargs[key]
            if key == 'startLine': self.startLine = kwargs[key]  
        self.__parseData()
    ############################################################################           
    ############################################################################  
    def __assignAnnotation(self, contig, candidate):
        counter = 0
        overlap = 0
        breakExit = False
        for el in self.hmmerAnnotation.get(contig, "None"):
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
        if self.hmmerAnnotation[contig][counter][3] > candidate[3]:
            if overlap > self.overlap: 
                self.hmmerAnnotation[contig][counter] = candidate
                return
            else:
                self.hmmerAnnotation[contig].append(candidate)
                return
        else:
            if overlap < self.overlap:
                self.hmmerAnnotation[contig].append(candidate)
                return
            else:
                return
    ############################################################################           
    ############################################################################
    def __dnaMetrics(self, frame, proteinStart, proteinEnd, contigName):
        orientation = ''
        forward, reverse = False, False
        if frame == '1': dnaStart, dnaEnd, orientation = proteinStart * 3 - 2, proteinEnd * 3, "forward"
        if frame == '2': dnaStart, dnaEnd, orientation = proteinStart * 3 - 1, proteinEnd * 3 + 1, "forward"
        if frame == '3': dnaStart, dnaEnd, orientation = proteinStart * 3, proteinEnd * 3 + 2, "forward"
        if frame == '4': dnaStart, dnaEnd, orientation = (self.__contigLengths[contigName] + 1)- proteinEnd * 3, self.__contigLengths[contigName] - (proteinStart - 1) * 3, "reverse" 
        if frame == '5': dnaStart, dnaEnd, orientation = self.__contigLengths[contigName]- proteinEnd * 3, self.__contigLengths[contigName] - (proteinStart - 1) * 3 - 1, "reverse" 
        if frame == '6': dnaStart, dnaEnd, orientation = self.__contigLengths[contigName]- proteinEnd * 3 - 1, self.__contigLengths[contigName] - (proteinStart - 1) * 3 - 2, "reverse"
        if dnaStart < 1: dnaStart = 1
        if dnaEnd > self.__contigLengths[contigName]: dnaEnd = self.__contigLengths[contigName]
        return (dnaStart, dnaEnd, orientation)          
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
                    contig = anno[0].split("_")[0]
                    frame = anno[0].split("_")[1]
                    alignment_start = float(anno[1])
                    alignment_end = float(anno[2])
                    score = float(anno[11])
                    e_val = float(anno[12])
                    ko = anno[6]
                    dnaCoords = self.__dnaMetrics(frame, alignment_start, alignment_end, contig)
                    if self.hmmerAnnotation.get(contig, "None") == "None":
                        self.hmmerAnnotation[contig].append((dnaCoords[0], dnaCoords[1], ko, e_val, score, dnaCoords[2]))
                    else:
                        if (e_val < self.evalue): self.__assignAnnotation(contig, (dnaCoords[0], dnaCoords[1], ko, e_val, score, dnaCoords[2]))
                except:
                    print "### print_exception:"
                    traceback.print_exc(hmmer_file=sys.stdout)
                    print anno
                    print "###"
                    
                    
                    
if __name__ == "__main__":
    data = "/media/Transcend/TAGS/Jura_meta/contigs/927-2E_test.out"
    contigDna = "/media/Transcend/TAGS/Jura_meta/contigs/927-2E_test.fna"
    x = HmmerTabParser(data, contigDna, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))                    
    for key, item in x.hmmerAnnotation.items():
        if key == 'contig00030':
            print "#######" + key + "#######"
            print [el[2] for el in item]
            for el in item:
                print el[0:3]
                print type(el[0])