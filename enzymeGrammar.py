'''
Created on Dec 12, 2011

@author: astar
'''
from pyparsing import Word, Group, Suppress, OneOrMore, alphas, nums,\
alphanums, Keyword, ParserElement, \
LineEnd, SkipTo, ZeroOrMore, Optional, oneOf, Literal
from persistent import Persistent
import re
#define features to extract
ParserElement.setDefaultWhitespaceChars(" \t")
ft_entry = Suppress(Keyword('ENTRY'))
ft_name = Suppress(Keyword('NAME'))
ft_class = Suppress(Keyword('CLASS'))
ft_sysname = Suppress(Keyword('SYSNAME'))
ft_reaction = Suppress(Keyword('REACTION'))
ft_allreac = Suppress(Keyword('ALL_REAC'))
ft_substrate = Suppress(Keyword('SUBSTRATE'))
ft_product = Suppress(Keyword('PRODUCT'))
ft_cofactor = Suppress(Keyword('COFACTOR'))
ft_inhibitor = Suppress(Keyword('INHIBITOR'))
ft_effector = Suppress(Keyword('EFFECTOR'))
ft_comment = Suppress(Keyword('COMMENT'))
ft_reference = Suppress(Keyword('REFERENCE'))
ft_authors = Suppress(Keyword('AUTHORS'))
ft_title = Suppress(Keyword('TITLE'))
ft_journal = Suppress(Keyword('JOURNAL'))
ft_pathway = Suppress(Keyword('PATHWAY'))
ft_orthology = Suppress(Keyword('ORTHOLOGY'))
ft_genes = Suppress(Keyword('GENES'))
ft_dblinks = Suppress(Keyword('DBLINKS'))
number = Word(nums).setParseAction(lambda tokens: int(tokens[0]))
undetermined = "undetermined"
#define ENTRY element
ec_number = Word(nums+".-")
ec_tag = Keyword("EC").suppress()
desc = OneOrMore(Word(alphas))
entry = ft_entry + ec_tag + ec_number("ec") + desc + LineEnd().suppress()
#define NAME element
aliases = SkipTo(oneOf("CLASS SYSNAME REACTION ALL_REAC SUBSTRATE PRODUCT COFACTOR INHIBITOR EFFECTOR COMMENT REFERENCE PATHWAY ORTHOLOGY GENES DBLINKS")).setParseAction(lambda feature: [[feat.lstrip().rstrip().replace(";", "") for feat in feature[0].split("\n") if feat]])
name = Optional(ft_name + aliases("aliases"))
#define CLASS element
class_text = SkipTo(oneOf("SYSNAME REACTION ALL_REAC SUBSTRATE PRODUCT COFACTOR INHIBITOR EFFECTOR COMMENT REFERENCE PATHWAY ORTHOLOGY GENES DBLINKS")).setParseAction(lambda feature: [[feat.lstrip().rstrip().replace(";", "") for feat in feature[0].split("\n") if feat]])
class_element = Optional(ft_class + class_text("class"))
#define SYSNAME element
sysname_text = SkipTo(oneOf("REACTION ALL_REAC SUBSTRATE PRODUCT COFACTOR INHIBITOR EFFECTOR COMMENT REFERENCE PATHWAY ORTHOLOGY GENES DBLINKS")).setParseAction(lambda feature: [[feat.lstrip().rstrip().split(":") for feat in feature[0].split("\n") if feat]])
sysname = Optional(ft_sysname + sysname_text("sysname"))
#define REACTION element
reaction_text = SkipTo(oneOf("ALL_REAC SUBSTRATE PRODUCT COFACTOR INHIBITOR EFFECTOR COMMENT REFERENCE PATHWAY ORTHOLOGY GENES DBLINKS")).setParseAction(lambda feature: [[feat.lstrip().rstrip().replace(";", "") for feat in feature[0].split("\n") if feat]])
reaction = Optional(ft_reaction + reaction_text("reaction"))
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
    f = open('/media/Transcend/Working_KEGG/KEGG_core/ko_genes', 'r')
    testString = f.read().split('///')
    count = 0
    for string in testString:
        if re.search(r'\w', string):
            try:
                count += 1
                x = KoEntry(string.lstrip().rstrip() + "///")
                elements = (x.ko_number, x.aliases, x.definition, x.pathways, x.diseases, x.modules, x.class_element, x.dblinks, x.genes, x.references)
                for el in elements:
                    print el
                if count == 21: sys.exit(1)
            except Exception, err:
                log.exception('Error from KoEntry():')
                sys.exit(1)