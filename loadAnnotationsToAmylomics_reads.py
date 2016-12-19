'''
Created on Mar 25, 2012

@author: astar
delete index http://localhost:8983/solr/update?stream.body=<delete><query>*:*</query></delete>&commit=true
'''
import re, sunburnt, os, sys
from lxml import etree
from StringIO import StringIO
from keggDB.getData import KeggData
from pyparsing import Keyword
from nltk.tokenize import word_tokenize, sent_tokenize
from parseHmmer_tab import HmmerTabParser
from collections import defaultdict
from CAZy.cazyData import CAZyMapper
from sqliteManager import sqliteManager
cazy = CAZyMapper()
keggKoAnnotation = KeggData('ko')  
conn = sqliteManager("/media/Transcend/TAGS/Jura_meta/reads/metagenomes.sqlite")
def __tokenize(text):
        tokens = []
        for sent in sent_tokenize(text):
            tokens.extend([re.sub('^-', '', word) for word in word_tokenize(sent) if not re.match(r"[\x90-\xff]", word) and re.search("\w", word)])
        if tokens: 
            return " ".join(tokens).encode('utf-8')
        else:
            return None
def __mapIdtoKegg(data):
    koMap = defaultdict(list)
    x = HmmerTabParser(data, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))
    for met_id, item in x.hmmerAnnotation.items():
        [koMap[el[2]].append(met_id) for el in item]
    return(dict([(key, item) for key, item in koMap.items()]))
def load(core, data, schema, hmm, library):
    koMap = __mapIdtoKegg(data)
    undetermined = 0
    solr_interface = sunburnt.SolrInterface(core, schema)
    #ids = ['kegg_orthologue', 'EC_number', 'CAZy_modules','read_hits', 'assembled_gene_hits', 'name','definition','class','pathway','pmid', 'global', 'specific', 'cazy_activities', 'cazy_specificities']   
    for ko, item in koMap.items():
        document = {}
        document["kegg_orthologue"] = ko
        document["read_hits"] = len(item)
        document["assembled_gene_hits"] = conn.getContigCount(ko, library)
        document["name"] = keggKoAnnotation.dbroot[ko].aliases
        document["definition"] = __tokenize(keggKoAnnotation.dbroot[ko].definition)
        ec = re.search('\[(.*)\]', keggKoAnnotation.dbroot[ko].definition)
        if ec: document["EC_number"] = ec.group(1)
        try:
            cazy_modules = cazy.keggMap[ko]
            document["CAZy_modules"] = tuple(cazy_modules)
            try:
                document["cazy_activities"] = __tokenize(" ".join([cazy.cazyDesc[module]['activities'] for module in cazy_modules]))
                document["cazy_specificities"] = __tokenize(" ".join([cazy.cazyDesc[module]['specificities'] for module in cazy_modules]))
            except:
                pass
        except:
            pass
        if keggKoAnnotation.dbroot[ko].definition == 'undetermined': undetermined += 1
        if keggKoAnnotation.dbroot[ko].class_element != 'undetermined':
            class_element = {}
            brite_global = {}
            brite_specific = {}
            for el in keggKoAnnotation.dbroot[ko].class_element:
                class_element[el.split('[')[0].rstrip()] = ''
                brite_global[el.split(';')[0].rstrip()] = ''
                brite_specific[el.split(';')[1].lstrip()] = ''
            document["class"] = tuple(class_element.keys())
            document["global"] = tuple(brite_global.keys())
            document["specific"] = tuple(brite_specific.keys())
        if keggKoAnnotation.dbroot[ko].pathways != 'undetermined':
            pathway = {}
            for path in keggKoAnnotation.dbroot[ko].pathways:
                pathway[path[1]] = ''
            document["pathway"] = tuple(pathway.keys())
        if keggKoAnnotation.dbroot[ko].references != 'undetermined':
            pubId = []
            for reference in keggKoAnnotation.dbroot[ko].references:
                pubId.append(reference[0])
            document["pmid"] = tuple(pubId)
        solr_interface.add(document)
    solr_interface.commit()
    print core + " - completed!"
    print undetermined
    print len(koMap)
def parseXML(xmlFile):
    coreMap = {}
    f = open(xmlFile)
    xml = f.read()
    f.close()
    tree = etree.parse(StringIO(xml))
    context = etree.iterparse(StringIO(xml))
    for action, elem in context:
        if elem.tag == 'core':
            coreMap[elem.values()[0]] = elem.values()[1]
    return coreMap
            
def main():
    coreURL = "http://localhost:8983/solr/"
    coreDir = "/media/Transcend/solr/apache-solr-3.6.0/example_amylomics/multicore"
    dataDir = "/media/Transcend/TAGS/Jura_meta/MET5_6_7" 
    configXML = "/media/Transcend/solr/apache-solr-3.6.0/example_amylomics/multicore/solr.xml"
    coreMap = parseXML(configXML)
    hmmer_files = [f for f in os.listdir(dataDir) if f.endswith('.out')]
    for hmm in hmmer_files:
        target = hmm.split(".")[0] + "_reads"
        core = coreURL + target
        schema = os.path.join(coreDir, coreMap[target], "conf/schema.xml")
        print schema
        print core
        data = os.path.join(dataDir, hmm)
        load(core, data, schema, hmm, target)
    keggKoAnnotation.close()  
    conn.close()
if __name__ == '__main__':
    main()
    
    
 
                    
        
