'''
Created on Mar 25, 2012

@author: astar
'''
import pickle, re, sunburnt, sys
from keggDB.getData import KeggData
from CAZy.cazyData import CAZyMapper
from nltk.tokenize import word_tokenize, sent_tokenize
from collections import defaultdict
def __tokenize(text):
        tokens = []
        for sent in sent_tokenize(text):
            tokens.extend([re.sub('^-', '', word) for word in word_tokenize(sent) if not re.match(r"[\x90-\xff]", word) and re.search("\w", word)])
        if tokens: 
            return " ".join(tokens)
        else:
            return None
def main():
    solr_interface = sunburnt.SolrInterface("http://localhost:1234/solr/core0/", "/media/Transcend/solr/lucidworks-solr-3.2.0_01/example/multicore/core0/conf/schema.xml")
    keggKoAnnotation = KeggData('ko')  
    cazys = CAZyMapper()
    pkl_file = open('read_data.pkl', 'rb')
    ko2contig2Reads = pickle.load(pkl_file)
    pkl_file.close()
    ids = ['ko_id','aliases','contig_num','read_num','ec_id','cazy_id','global_cat','process_cat','pathway_cat','ko_definition','cazy_activities','cazy_specificities']   
    categories = ['global_cat','process_cat','pathway_cat']
    counts = 0
    counts2 = 0
    for ko, val in ko2contig2Reads.items():
        counts += 1
        document = dict(zip(ids, [None]*11))
        keggDataObj = keggKoAnnotation.dbroot[ko]
        document["ko_id"] = ko
        ecMatch = re.search('\[(.*)\]', keggDataObj.definition)
        if ecMatch: document["ec_id"] = ecMatch.group(1)
        try:
            cid = cazys.keggMap[ko]
            counts2 += 1
            document["cazy_id"] = cid
            try:
                document["cazy_activities"] = cazys.cazyDesc[ko]['activities']
                document["cazy_specificities"] = cazys.cazyDesc[ko]['specificities']
            except:
                pass
        except:
            pass
        document["ko_definition"] = __tokenize(keggDataObj.definition)
        document['aliases'] = keggDataObj.aliases
        pathways =  defaultdict(list)
        for brite in keggDataObj.class_element:
            for el in enumerate(brite.split(";")):
                pathways[categories[el[0]]].append(el[1])
        for key, val in pathways.items():
            document[key] = tuple(val)
        document['contig_num'] = 0
        document['read_num'] = 0
        for contig, value in ko2contig2Reads[ko].items():
            document['read_num'] += len(value)
            if re.match('contig.*', contig): document['contig_num'] += 1
        solr_interface.add(document)
        solr_interface.commit()
    print counts   
    print counts2
if __name__ == '__main__':
    main()      
                    
        
