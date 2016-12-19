'''
Created on Apr 12, 2013

@author: astar
'''
import csv, collections, sys
cr1 = csv.reader(open("/media/Transcend/workDir/PaulPrograms/src/zoophyteBase/gh13_cazy_motifs.txt","rb"))
cr2 = csv.reader(open("/media/Transcend/workDir/PaulPrograms/src/zoophyteBase/gh13_all_motifs.txt","rb"))
container2 = collections.defaultdict(list)
container1 = collections.defaultdict(list)
for el in cr2:
    container2[el[1]].append(el[0])
for el in cr1:
    container1[el[1]].append(el[0])
trehalose_hydrolase = set("1217 - 1012 - 1152 - 1013 - 1587 - 1586".split(" - "))
test = set(container1["1217 - 1012 - 1152 - 1013 - 1587 - 1586"])
complexity = sorted([(set(el[0].split(" - ")), el[1]) for el in container1.items()], key=lambda x: len(x[1]), reverse=True)
#counter = 0
#for pos1, el1 in enumerate(complexity):
#    for pos2, el2 in enumerate(complexity):
#        if pos1 != pos2:
#            smaller = len(el1[0])
#            if smaller > len(el2[0]): smaller = len(el2[0])
#            if float(len(el1[0].intersection(el2[0])))/float(smaller) > 0.75:
#                print el1
#                print el2
#                counter += 1
#                print "############################################"
#                if counter == 30: sys.exit(1)
#print "######################################################################"               
for el in sorted(container1.items(), key=lambda x: len(x[1]), reverse=True):
    candidate = set(el[0].split(" - "))
    if len(candidate.intersection(trehalose_hydrolase)) > 4:
        print el
