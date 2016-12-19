'''
Created on Mar 1, 2011

@author: astar
'''
from pyparsing import Word, Group, Suppress, OneOrMore, alphas, nums,\
alphanums, Keyword, ParserElement, \
LineEnd, SkipTo, ZeroOrMore, Optional, oneOf, Literal
from persistent import Persistent
import re
def trimFeature(feature):
    trimmedFeature = feature[0].rstrip()
    trimmedFeature = re.sub(r'\n', ' ', trimmedFeature)
    trimmedFeature = re.sub(r'\s+', ' ', trimmedFeature)
    return trimmedFeature   
#define features to extract
ParserElement.setDefaultWhitespaceChars(" \t")
ft_entry = Suppress(Keyword('ENTRY'))
ft_name = Suppress(Keyword('NAME'))
ft_definition = Suppress(Keyword('DEFINITION'))
ft_pathway = Suppress(Keyword('PATHWAY'))
ft_disease = Suppress(Keyword('DISEASE'))
ft_module = Suppress(Keyword('MODULE'))
ft_class = Suppress(Keyword('CLASS'))
ft_dblinks = Suppress(Keyword('DBLINKS'))
ft_genes = Suppress(Keyword('GENES'))
ft_reference = Suppress(Keyword('REFERENCE'))
ft_authors = Suppress(Keyword('AUTHORS'))
ft_title = Suppress(Keyword('TITLE'))
ft_journal = Suppress(Keyword('JOURNAL'))
number = Word(nums).setParseAction(lambda tokens: int(tokens[0]))
undetermined = "undetermined"
#define ENTRY element
ko_number = Word('K', alphanums)
ko_tag = Word(alphas)
entry = ft_entry + ko_number("ko_number") + ko_tag + LineEnd().suppress()
#define NAME element
aliases = SkipTo(oneOf("DEFINITION PATHWAY DISEASE MODULE CLASS DBLINKS GENES REFERENCE")).setParseAction(trimFeature)
name = Optional(ft_name + aliases("aliases"))
#define DEFINITION element
definition_text = SkipTo(oneOf("PATHWAY DISEASE MODULE CLASS DBLINKS GENES REFERENCE")).setParseAction(trimFeature)
definition = Optional(ft_definition + definition_text("definition"))
#define PATHWAY element
pathway_text = SkipTo(oneOf("DISEASE MODULE CLASS DBLINKS GENES REFERENCE")).setParseAction(lambda feature: [[el.lstrip().rstrip() for el in feat.lstrip().rstrip().split(" ", 1) if el.lstrip().rstrip()] for feat in feature[0].split("\n") if feat])
pathway = Optional(ft_pathway + pathway_text)("pathways")
#define DISEASE element
disease_text = SkipTo(oneOf("MODULE CLASS DBLINKS GENES REFERENCE")).setParseAction(lambda feature: [[el.lstrip().rstrip() for el in feat.split(" ", 1) if el.lstrip().rstrip()] for feat in feature[0].split("\n") if feat])
disease = Optional(ft_disease + disease_text)("diseases")
#define MODULE element
module_text = SkipTo(oneOf("CLASS DBLINKS GENES REFERENCE")).setParseAction(lambda feature: [[el.lstrip().rstrip() for el in feat.split(" ", 1) if el.lstrip().rstrip()] for feat in feature[0].split("\n") if feat])
module = Optional(ft_module + module_text)("diseases")
#define CLASS element
class_text = SkipTo(oneOf("DBLINKS GENES REFERENCE")).setParseAction(lambda feature: [feat.lstrip().rstrip() for feat in feature[0].split("\n") if feat])
class_element = Optional(ft_class + class_text)("class_element")
#define DBLINKS element
dblinks_text = SkipTo(oneOf("GENES REFERENCE")).setParseAction(lambda feature: [[el.lstrip().rstrip().split() for el in feat.lstrip().rstrip().split(":", 1) if el.lstrip().rstrip().split()] for feat in feature[0].split("\n") if feat])
dblinks = Optional(ft_dblinks + dblinks_text)("dblinks")
#define GENES element
genes_text = SkipTo(oneOf("REFERENCE ///")).setParseAction(lambda feature: [[el.lstrip().rstrip().split() for el in feat.lstrip().rstrip().split(":", 1) if el.lstrip().rstrip().split()] for feat in feature[0].split("\n") if feat])
genes = Optional(ft_genes + genes_text)("genes")
#define REFERENCE element
reference_id = ft_reference + Literal("PMID:").suppress() + number + Optional(OneOrMore(Word(alphanums+"()_.")).setParseAction(lambda tokens: " ".join(tokens))) +  LineEnd().suppress() 
authors = ft_authors + SkipTo(oneOf("TITLE JOURNAL")).setParseAction(trimFeature)
title = ft_title + SkipTo(oneOf("JOURNAL REFERENCE")).setParseAction(trimFeature)
journal = ft_journal + SkipTo("///").setParseAction(trimFeature)
reference = ZeroOrMore(Group(reference_id + authors + title + journal))("references")
#define KO entry element
annotations = ["ko_number", "aliases", "definition", "pathways", "diseases", "modules", "class_element", "dblinks", "genes", "references"]
annotation_exp = ["self." + annotation + " = " + undetermined for annotation in annotations]
ko_entry_elements = [entry, name, definition, pathway, disease, module, class_element, dblinks, genes, reference]
ko_entry_parse_expr = ko_entry_elements[0]
for feat in ko_entry_elements[1:]:
    ko_entry_parse_expr += feat
class KoEntry(Persistent):
    def __init__(self, ko):
        self.ko = ko
        self.__parseKo()
    def __parseKo(self):
        for exp in annotation_exp: exec exp
        parsedResultList = ko_entry_parse_expr.parseString(self.ko)
        for key in parsedResultList.keys(): 
            exp = "self." + key + " = parsedResultList[key]"
            exec exp    
if __name__ == "__main__":
    import logging, sys, re
    logging.basicConfig(level=logging.WARNING)
    log = logging.getLogger('/media/Transcend/Working_KEGG/KEGG_core/log')
    f = open('/media/Transcend/Working_KEGG/genes/ko', 'r')
    testString = f.read().split('///')
    count = 0
    for string in testString:
        if re.search(r'\w', string):
            try:
                #count += 1
                x = KoEntry(string.lstrip().rstrip() + "///")
                elements = (x.ko_number, x.aliases, x.definition, x.pathways, x.diseases, x.modules, x.class_element, x.dblinks, x.genes, x.references)
                for el in elements:
                    if el != "undetermined":
                        pass
                        #print el
                        #print "\n".join(x.references)
                if x.ko_number == "K02040": 
                    print x.definition
                    print x.pathways
                    print x.class_element
                    sys.exit(1)
            except Exception, err:
                log.exception('Error from KoEntry():')
                sys.exit(1)