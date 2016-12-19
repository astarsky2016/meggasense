'''
Created on Apr 18, 2012

@author: astar
'''
from keggDB.getData import KeggData
from csvWriter import UnicodeWriter
def main():
    out = open('kegg_summary.csv', 'w')
    keggKoAnnotation = KeggData('ko')
    writer = UnicodeWriter(out)
    for ko in keggKoAnnotation.dbroot.keys():
        row = []
        row.append(ko)
        keggDataObj = keggKoAnnotation.dbroot[ko]
#        row.append(keggDataObj.definition)
#        print keggDataObj.definition
#        print keggDataObj.class_element
#        pathway = ''
#        for el in keggDataObj.pathways:
#            pathway += " ".join(el) + " "
#        row.append(pathway)
        if "".join(keggDataObj.class_element) != "undetermined":
            row.append("|".join(keggDataObj.class_element))
        else:
            row.append("".join(keggDataObj.class_element))
        writer.writerow(row)
    out.close()
    keggKoAnnotation.close()
if __name__ == '__main__':
    main()
    