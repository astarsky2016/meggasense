'''
Created on Feb 2, 2012

@author: astar
'''
import os, sys, inspect
cmd_folder = os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0])
kegg_folder = "/".join(cmd_folder.split("/")[0:-3]) + "/KEGG/src"
if kegg_folder not in sys.path:
    sys.path.append(kegg_folder)
from storage.myZODB import MyZODB
class KeggData:
    data = {'ko': ('genes', 'koGrammar', 'KoEntry', '/media/Transcend/TAGS/myZODB/keggKOs.fs'), 
            'genome': ('genes', 'genomeGrammar', 'GenomeEntry', '/media/Transcend/TAGS/myZODB/keggGenomes.fs'),
            'brite': ('brite', 'koBriteGrammar', 'koBriteEntry', '/media/Transcend/TAGS/myZODB/keggBrite.fs')}
    def __init__(self, requestedData):
        exec('from %s.%s import %s' % KeggData.data[requestedData][0:-1])
        from genes.koGrammar import KoEntry
        self.__db = MyZODB(KeggData.data[requestedData][-1])
        self.dbroot = self.__db.dbroot
    def close(self):
        self.__db.close()
  
     