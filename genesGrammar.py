'''
Created on Nov 29, 2011

@author: astar
'''
from pyparsing import Word, Group, Suppress, OneOrMore, alphas, nums,\
alphanums, Keyword, ParserElement, LineEnd, SkipTo, ZeroOrMore, Optional, oneOf
from persistent import Persistent
def trimFeature(feature):
    import re
    trimmedFeature = feature[0].rstrip()
    trimmedFeature = re.sub(r'\n', ' ', trimmedFeature)
    trimmedFeature = re.sub(r'\s+', ' ', trimmedFeature)
    return trimmedFeature   
#define features to extract
ParserElement.setDefaultWhitespaceChars(" \t")
ft_entry = Suppress(Keyword('ENTRY'))
ft_name = Suppress(Keyword('NAME'))
ft_definition = Suppress(Keyword('DEFINITION'))
ft_orthology = Suppress(Keyword('ORTHOLOGY'))
ft_pathway = Suppress(Keyword('PATHWAY'))
ft_class = Suppress(Keyword('CLASS'))
ft_position = Suppress(Keyword('POSITION'))
ft_motif= Suppress(Keyword('MOTIF'))
ft_dblinks = Suppress(Keyword('DBLINKS'))
ft_structure = Suppress(Keyword('STRUCTURE'))
ft_aaseq = Suppress(Keyword('AASEQ'))
ft_ntseq = Suppress(Keyword('NTSEQ'))
number = Word(nums).setParseAction(lambda tokens: int(tokens[0]))
undetermined = "undetermined"
#define ENTRY element
entry = ft_entry + OneOrMore(Word(alphanums+"_.-'"))("entry") + LineEnd().suppress()
#define NAME element
name_text = SkipTo(oneOf("DEFINITION ORTHOLOGY PATHWAY CLASS POSITION MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [[feat.lstrip().rstrip() for feat in feature[0].lstrip().rstrip().split(",")]])
name = Optional(ft_name + name_text("name"))
#define DEFINITION element
definition_text = SkipTo(ft_orthology | ft_pathway | ft_class | ft_position | ft_motif | ft_dblinks | ft_structure | ft_aaseq | ft_ntseq).setParseAction(trimFeature)
definition = Optional(ft_definition + definition_text("definition"))
#define ORTHOLOGY element
orthology = Optional(ft_orthology + SkipTo(oneOf("PATHWAY CLASS POSITION MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [[el.lstrip().rstrip() for el in feat.lstrip().rstrip().split(" ", 1) if el.lstrip().rstrip()] for feat in feature[0].split("\n") if feat]))("ko_orthology")
#define PATHWAY element
pathway_text = SkipTo(oneOf("CLASS POSITION MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [[el.lstrip().rstrip() for el in feat.lstrip().rstrip().split(" ", 1) if el.lstrip().rstrip()] for feat in feature[0].split("\n") if feat])
pathway = Optional(ft_pathway + pathway_text)("pathways")
#define CLASS element
class_text = SkipTo(oneOf("POSITION MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [feat.lstrip().rstrip() for feat in feature[0].split("\n") if feat])
class_element = Optional(ft_class + class_text)("class_element")
#define POSITION element
position = Optional(ft_position + SkipTo(oneOf("MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")))("positions")
#define MOTIF element
motifs = Optional(ft_motif + SkipTo(oneOf("DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [feat.lstrip().rstrip().split() for feat in feature[0].split("\n") if feat]))("motifs")
#define DBLINKS element
dblinks_text = SkipTo(oneOf("STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [[el.lstrip().rstrip().split() for el in feat.lstrip().rstrip().split(":", 1) if el.lstrip().rstrip().split()] for feat in feature[0].split("\n") if feat])
dblinks = Optional(ft_dblinks + dblinks_text)("dblinks")
#define STRUCTURE element
structure = Optional(ft_structure + SkipTo(oneOf("AASEQ NTSEQ")).setParseAction(lambda feature: [feat.lstrip().rstrip().split() for feat in feature[0].split("\n")]))("structures")
#define AASEQ element
aaseq = Optional(ft_aaseq + number + LineEnd().suppress() + SkipTo(Keyword("NTSEQ")).setParseAction(lambda features: "".join([feat.lstrip().rstrip() for feat in features[0].split("\n")])))("aaseq")
#define NTSEQ element
ntseq = Group(ft_ntseq + number + LineEnd().suppress() + SkipTo("///").setParseAction(lambda features: "".join([feat.lstrip().rstrip() for feat in features[0].split("\n")])))("ntseq")
#define GENE entry element
annotations = ["entry", "name", "definition", "orthology", "pathways", "class_element", "motifs", "dblinks", "structures"]
annotation_exp = ["self." + annotation + " = " + undetermined for annotation in annotations]
annotation_exp.extend(["self.positions = undetermined", "self.complement = undetermined"])
gene_entry_elements = [entry, name, definition, orthology, pathway, class_element, position, motifs, dblinks, structure]
gene_entry_parse_expr = gene_entry_elements[0]
for feat in gene_entry_elements[1:]:
    gene_entry_parse_expr += feat
class GeneEntry(Persistent):
    def __init__(self, gene):
        self.__parseGene(gene)     
    def __extractPositions(self, feature):
        import re
        if feature:
            complement = False
            matchPos = re.findall(r'\d+', feature)
            matchComp = re.search(r'complement', feature, re.IGNORECASE)
            if matchComp: complement = True
            positions = []
            while matchPos:
                positions.append(matchPos[0:2])
                matchPos = matchPos[2:]
            self.positions = positions
            self.complement = complement
        else:
            complement = undetermined
            self.positions = undetermined
            self.complement = undetermined
    def __parseGene(self, gene):
        for exp in annotation_exp: exec exp
        parsedResultList = gene_entry_parse_expr.parseString(gene)
        for key in parsedResultList.keys(): 
            if key == 'positions':
                if parsedResultList[key]: self.__extractPositions(parsedResultList[key][0])
                else: self.__extractPositions('')
            else:
                exp = "self." + key + " = parsedResultList[key]"
                exec exp  
if __name__ == "__main__":
    import logging, re, os, sys
    logging.basicConfig(level=logging.WARNING)
    log = logging.getLogger('/media/Transcend/Working_KEGG/KEGG_core/log')
    for line in open('/media/Transcend/Working_KEGG/KEGG_core/genes/organisms_test.lst', 'r'):
        count = 0
        if re.search(r'\w', line):
            organism = os.path.join("/media/Transcend/Working_KEGG/genes", line.lstrip().rstrip())
            try:
                testString = open(organism, 'r').read().split('///')
            except IOError:
                print "#####no such file: " + organism + "#####"
                continue
            for string in testString:
                if re.search(r'\w', string):
                    try:
                        count += 1
                        geneDoc = {}
                        x = GeneEntry(string.lstrip().rstrip() + "///")
                        print dir(x)
                        elements = (x.entry, x.name, x.definition, x.orthology, x.pathways, x.class_element, x.complement, x.positions, 
                                    x.motifs, x.dblinks, x.structures)
                        for el in elements:
                            print el
                        if count == 2: sys.exit(1)
                    except Exception, err:
                        log.exception('Error from GeneEntry():')
                        print organism
                        print string
                        continue