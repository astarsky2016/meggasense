'''
Created on Oct 17, 2012

@author: astar
'''
import re
from csvWriter import UnicodeWriter
inputFile = "/media/Transcend/TAGS/Jura_meta/MET5_6_7/454NewblerMetrics.txt"
f = open(inputFile)
reg1 = re.compile('\s+numberOfReads\s+=\s+(\d+)(, ){0,1}(\d+){0,1}.*')
reg2 = re.compile('\s+numAlignedReads\s+=\s+(\d+)(, ){0,1}(\d+\.\d+%){0,1}.*')
reg3 = re.compile('\s+numberSingleton\s+=\s+(\d+).*')
reg4 = re.compile('\s+numberOfContigs\s+=\s+(\d+).*')
reg5 = re.compile('\s+largestContigSize\s+=\s+(\d+).*')
reg6 = re.compile('\s+N50ContigSize\s+=\s+(\d+).*')
reg7 = re.compile('\s+numberOfBases\s+=\s+(\d+).*')
flag1, flag2 = False, False
columns = {}
for line in f:
    if reg1.match(line):
        columns[1] = ("Number of reads", reg1.match(line).group(3))
    if reg2.match(line):
        columns[2] = ("% aligned reads", reg2.match(line).group(3))
    if reg3.match(line):
        columns[3] = ("Number of singletons", reg3.match(line).group(1))
    if flag1:
        if reg4.match(line):
            columns[4] = ("Number of large contigs", reg4.match(line).group(1))
        if reg5.match(line):
            columns[5] = ("Largest contig (bp)", reg5.match(line).group(1)) 
            flag1 = False
        if reg6.match(line):
            columns[6] = ("N50 contig size (bp)", reg6.match(line).group(1))
    if flag2:
        if reg4.match(line):
            columns[7] = ("Number of All contigs", reg4.match(line).group(1))
        if reg7.match(line):
            columns[8] = ("Total number of bases", reg7.match(line).group(1))
    if re.search('largeContigMetrics', line):
        flag1 = True
    if re.search('allContigMetrics', line):
        flag2 = True 
f.close()
out = open("/media/Transcend/TAGS/Jura_meta/MET5_6_7/454NewblerMetrics.csv", "w")
writer = UnicodeWriter(out)
for el in range(1,9):
    writer.writerow(columns[el])
out.close()

    