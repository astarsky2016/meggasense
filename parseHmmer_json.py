'''
Created on Jan 24, 2012

@author: astar
'''
import simplejson, sys
f = open("/media/Transcend/workDir/PaulPrograms/src/adigi.json").readlines()
for line in f[0:100]:
    print line
f.close()
#f = simplejson.load(open("/media/Transcend/workDir/PaulPrograms/src/adigi.json"))
#kos = []
#for el in f:
#    print el['name'] + " : " + el['evalue']
#    kos.append(el['name'])
#    print el['seq']['name'] + " from: " + el['env']['from'] + " - to: " + el['env']['to']
#print len(f)
#print kos
    