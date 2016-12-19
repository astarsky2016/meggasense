'''
Created on Jan 26, 2011

@author: astar
'''
from pyparsing import Word, Group, Suppress, OneOrMore, alphas, nums,\
alphanums, Keyword, commaSeparatedList, ParserElement, \
LineEnd, SkipTo, ZeroOrMore, Optional, Literal, oneOf
from persistent import Persistent
def trimFeature(feature):
        import re
        trimmedFeature = feature[0].rstrip().lstrip()
        trimmedFeature = re.sub(r'\n', ' ', trimmedFeature)
        trimmedFeature = re.sub(r'\s+', ' ', trimmedFeature)
        return trimmedFeature   
#define features to extract
ParserElement.setDefaultWhitespaceChars(" \t")
ft_entry = Suppress(Keyword('ENTRY'))
ft_name = Suppress(Keyword('NAME'))
ft_definition = Suppress(Keyword('DEFINITION'))
ft_annotation = Suppress(Keyword('ANNOTATION'))
ft_taxonomy = Suppress(Keyword('TAXONOMY'))
ft_lineage = Suppress(Keyword('LINEAGE'))
ft_data_source = Suppress(Keyword('DATA_SOURCE'))
ft_original_db = Suppress(Keyword('ORIGINAL_DB'))
ft_keywords = Suppress(Keyword('KEYWORDS'))
ft_disease = Suppress(Keyword('DISEASE'))
ft_comment = Suppress(Keyword('COMMENT'))
ft_chromosome = Suppress(Keyword('CHROMOSOME'))
ft_plasmid = Suppress(Keyword('PLASMID'))
ft_sequence = Suppress(Keyword('SEQUENCE'))
ft_length = Suppress(Keyword('LENGTH'))
ft_statistics = Suppress(Keyword('STATISTICS'))
ft_number_of_dna = Suppress(Keyword('Number of nucleotides:'))
ft_number_of_prot_genes = Suppress(Keyword('Number of protein genes:'))
ft_number_of_rna_genes = Suppress(Keyword('Number of RNA genes:'))
ft_number_of_contigs = Suppress(Keyword('Number of contigs:'))
ft_number_of_singletons = Suppress(Keyword('Number of singletons:'))
ft_reference = Suppress(Keyword('REFERENCE'))
ft_authors = Suppress(Keyword('AUTHORS'))
ft_title = Suppress(Keyword('TITLE'))
ft_journal = Suppress(Keyword('JOURNAL'))
number = Word(nums).setParseAction(lambda tokens: int(tokens[0]))
undetermined = "undetermined"
#define ENTRY element
t_number = Word('T', nums)
genome_integrity = ZeroOrMore(Word(alphas))
genome_integrity.setParseAction(lambda tokens: " ".join(tokens))
entry = Optional(ft_entry + t_number("t_number") + genome_integrity("genome_integrity") + LineEnd().suppress())
#define NAME element
aliases = commaSeparatedList(OneOrMore(Word(alphanums+"._' ")))
name = Optional(ft_name + aliases("aliases") + LineEnd().suppress())
#define DEFINITION element
definition_text = SkipTo(oneOf("ANNOTATION TAXONOMY LINEAGE DATA_SOURCE ORIGINAL_DB KEYWORDS DISEASE COMMENT CHROMOSOME PLASMID STATISTICS REFERENCE"))
definition_text.setParseAction(trimFeature)
definition = Optional(ft_definition + definition_text("definition"))
#define ANNOTATION element
annotation_content = OneOrMore(Word(alphas))
annotation = Optional(ft_annotation + annotation_content("annotation").setParseAction(lambda features: " ".join(features)) + LineEnd().suppress())
#define TAXONOMY element
tax_id = Literal("TAX:") + number("taxonomy")
taxonomy = Optional(ft_taxonomy + tax_id + LineEnd().suppress())
#define LINEAGE element
lineage = Optional(ft_lineage + SkipTo(oneOf("DATA_SOURCE ORIGINAL_DB KEYWORDS DISEASE COMMENT CHROMOSOME PLASMID STATISTICS REFERENCE"))("lineage").setParseAction(trimFeature))
#define DATA_SOURCE element
data_source_db = OneOrMore(Word(alphanums+"():-_."))
data_source = Optional(ft_data_source + data_source_db("data_source").setParseAction(lambda tokens: " ".join(tokens)) + LineEnd().suppress())
#define ORIGINAL_DB element
original_db = SkipTo(oneOf("KEYWORDS DISEASE COMMENT CHROMOSOME PLASMID STATISTICS REFERENCE")).setParseAction(lambda feature: [[feat.lstrip().rstrip() for feat in feature[0].lstrip().rstrip().split("\n")]])
original_dbs = Optional(ft_original_db + original_db("original_dbs")) 
#define KEYWORDS element
keyword = OneOrMore(Word(alphas+",;")).setParseAction(lambda tokens: " ".join(tokens))
keywords = Optional(ft_keywords + keyword("keywords") + LineEnd().suppress())
#define DISEASE element
h_number = Word('H', nums)
disease_desc = OneOrMore(Word(alphas+"()")).setParseAction(lambda tokens: " ".join(tokens))
disease = Optional(ft_disease + h_number("h_number") + disease_desc("disease") + LineEnd().suppress()) 
#define COMMENT element
comment = Optional(ft_comment + SkipTo(oneOf("CHROMOSOME PLASMID STATISTICS REFERENCE"))("comment").setParseAction(trimFeature)) 
#define CHROMOSOME element
chromosome = ZeroOrMore(Group(ft_chromosome + OneOrMore(Word(alphanums+";:()_")).setParseAction(lambda tokens: " ".join(tokens)) + LineEnd().suppress() \
+ ft_sequence + Literal('RS:').suppress() + Word(alphanums+"_") + LineEnd().suppress() + ft_length + number + LineEnd().suppress()))("chromosomes")
#define PLASMID element
plasmid = ZeroOrMore(Group(ft_plasmid + OneOrMore(Word(alphanums+";:()_")).setParseAction(lambda tokens: " ".join(tokens)) + LineEnd().suppress() \
+ ft_sequence + Literal('RS:').suppress() + Word(alphanums+"_") + LineEnd().suppress() + ft_length + number + LineEnd().suppress()))("plasmids")
#define STATISTICS element
statistics = Group(Optional(ft_statistics + Optional(ft_number_of_dna + number("dna_length") + LineEnd().suppress()) + \
Optional(ft_number_of_prot_genes + number("number_of_proteins") + LineEnd().suppress()) + \
Optional(ft_number_of_rna_genes + number("number_of_rna_genes") + LineEnd().suppress()) + \
Optional(ft_number_of_contigs + number("number_of_contigs") + LineEnd().suppress()) + \
Optional(ft_number_of_singletons + number("number_of_singletons") + LineEnd().suppress())))("statistics")
#define REFERENCE element
reference_id = ft_reference + Literal("PMID:").suppress() + number + Optional(OneOrMore(Word(alphanums+"()_.")).setParseAction(lambda tokens: " ".join(tokens))) +  LineEnd().suppress() 
authors = ft_authors + SkipTo(oneOf("TITLE JOURNAL")).setParseAction(trimFeature)
title = ft_title + SkipTo(oneOf("JOURNAL REFERENCE")).setParseAction(trimFeature)
journal = ft_journal + OneOrMore(Word(alphanums+":-(),;_|")).setParseAction(lambda tokens: " ".join(tokens)) + LineEnd().suppress()
reference = ZeroOrMore(Group(reference_id + authors + title + journal))("references")
#define GENOME entry element
genome_entry_elements = [entry, name, definition, annotation, taxonomy, lineage, data_source, original_dbs, 
                keywords, disease, comment, chromosome, plasmid, statistics, reference]
genome_entry_parse_expr = genome_entry_elements[0]
for feat in genome_entry_elements[1:]:
    genome_entry_parse_expr += feat
class GenomeEntry(Persistent):
    def __init__(self, genome):
        self.genome = genome
        self.__parseGenome()
    def __parseGenome(self):
        annotations = ["t_number", "genome_integrity", "aliases", "definition", "annotation", "taxonomy", "lineage", "data_source", 
                       "original_dbs", "keywords", "h_number", "disease", "comment", "chromosomes", "plasmids", "statistics", "references"]
        for exp in ["self." + annotation + " = " + undetermined for annotation in annotations]: exec exp
        parsedResultList = genome_entry_parse_expr.parseString(self.genome)
        for key in parsedResultList.keys(): 
            exp = "self." + key + " = parsedResultList[key]"
            exec exp    
if __name__ == "__main__":
    import logging, sys
    logging.basicConfig(level=logging.WARNING)
    log = logging.getLogger('/media/Transcend/Working_KEGG/KEGG_core/log')
    f = open('/media/Transcend/Working_KEGG/KEGG_core/genes_genome', 'r')
    testString = f.read().split('///')
    for string in testString:
        try:
            x = GenomeEntry(string.lstrip().rstrip())
            elements = (x.t_number, x.genome_integrity, x.aliases, x.definition, x.annotation, x.taxonomy, x.lineage, x.data_source, x.original_dbs, 
                        x.keywords, x.h_number, x.disease, x.comment, x.chromosomes, x.plasmids, x.statistics, x.references)
            for el in elements:
                print el
        except Exception, err:
            log.exception('Error from GenomeEntry():')
            sys.exit(1)
      
        
       