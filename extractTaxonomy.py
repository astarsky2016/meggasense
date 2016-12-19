'''
Created on Apr 15, 2013

@author: astar
'''
from collections import defaultdict
names = open("/home/astar/BioApps/KronaTools-2.3/taxonomy/names.dmp")
nodes = open("/home/astar/BioApps/KronaTools-2.3/taxonomy/nodes.dmp")
taxDict = defaultdict(dict)
for line in names:
    elements = line.split("|")
    id = int(elements[0].rstrip())
    orgn = elements[1].lstrip().rstrip()
    taxDict["name"][id] = orgn
names.close()
nodesDict = {}
for line in nodes:
    elements = line.split("|")
    id = int(elements[0].rstrip())
    parent = elements[1].lstrip().rstrip()
    rank = elements[2].lstrip().rstrip()
    taxDict["parent"] = parent
    taxDict["rank"] = rank
nodes.close()
out = open("/home/astar/BioApps/KronaTools-2.3/taxonomy/nodes.dmp", "w")
for el in sorted(taxDict.keys()):
    print 
    
 
