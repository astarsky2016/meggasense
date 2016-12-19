'''
Created on Nov 28, 2011

@author: astar
'''
a = """ENTRY       AT4G01540         CDS       A.thaliana
NAME        NTM1
DEFINITION  NTM1 (NAC WITH TRANSMEMBRANE MOTIF1); transcription activator/ transcription factor
POSITION    4
MOTIF       Pfam: NAM
            PROSITE: NLS_BP NAC
DBLINKS     NCBI-GI: 145332955
            NCBI-GeneID: 828149
            MIPS: AT4G01540.1
            TAIR: AT4G01540
            UniProt: A8MQY1 C0SVG8
AASEQ       423
            MMKGLIGYRFSPTGEEVINHYLKNKLLGKYWLVDEAISEINILSHKPSKDLPKLARIQSE
            DLEWYFFSPIEYTNPNKMKMKRTTGSGFWKPTGVDREIRDKRGNGVVIGIKKTLVYHEGK
            SPHGVRTPWVMHEYHITCLPHHKRKYVVCQVKYKGEAAEISYEPSPSLVSDSHTVIAITG
            EPEPELQVEQPGKENLLGMSVDDLIEPMNQQEEPQGPHLAPNDDEFIRGLRHVDRGTVEY
            LFANEENMDGLSMNDLRIPMIVQQEDLSEWEGFNADTFFSDNNNNYNLNVHHQLTPYGDG
            YLNAFSGYNEGNPPDHELVMQENRNDHMPRKPVTGTIDYSSDSGSDAGSISTTVKQEIPR
            AVDAPMNNESSLVKTEKKGLFIVEDAMERNRKKPRFIYLMKMIIGNIISVLLPVKRLIPV
            KKL
NTSEQ       1272
            atgatgaaaggtctgattgggtatagatttagtccgacgggagaggaagtgatcaaccat
            tacctaaagaacaaacttctgggtaagtattggctcgttgatgaagctattagcgagatc
            aacatcttgagtcacaaacccagcaaggatttgcctaagttagctaggatccaatcggaa
            gatcttgaatggtatttcttctctccgattgagtacacgaacccgaataagatgaaaatg
            aagaggacgacaggttctgggttttggaaacctactggtgttgatcgggaaattagggat
            aaaagaggaaatggtgttgtgatagggattaagaagacgcttgtgtaccatgaaggtaag
            agtcctcatggagttagaactccttgggttatgcacgagtatcacatcacttgcttgcct
            catcataagaggaaatatgttgtctgccaagtaaagtataagggtgaagctgcagaaatt
            tcatatgagccaagtccctctttggtatccgattcgcataccgtcatagcgattaccgga
            gaaccggaacctgagcttcaggttgagcagccaggtaaagaaaatctcttgggtatgtct
            gtagatgatttgatagaaccaatgaaccaacaagaggagccacaaggtcctcacttagct
            ccgaatgatgatgagtttatacgtggattgaggcatgttgatcgagggacggttgaatat
            ttgtttgccaatgaagaaaacatggatggtttgtctatgaatgacttgagaatcccaatg
            atcgtccaacaagaggatctctctgagtgggagggatttaacgcagacacctttttcagc
            gacaacaacaataactataaccttaacgtgcatcatcaactaacgccttacggcgatggc
            tatttgaatgcattttcgggttataacgaagggaatcctcccgatcacgaattagtgatg
            caagagaaccgcaacgatcacatgccaaggaaacctgtgacagggaccattgattatagc
            agcgatagtggcagtgatgctggatccatatctacaacggtgaaacaagaaatcccaaga
            gctgttgatgcacccatgaacaatgagtcatctttggtgaaaacagagaagaaaggcttg
            tttattgtagaggacgcaatggagagaaaccgcaagaaaccacgatttatctatctcatg
            aagatgatcataggcaacatcatatcggttttactacccgtcaaaagattgatcccggtg
            aagaagttatga
///"""
def extractPositions(feature):
    import re
    complement = False
    matchPos = re.findall(r'\d+', feature)
    matchComp = re.search(r'complement', feature, re.IGNORECASE)
    if matchComp: complement = True
    positions = []
    while matchPos:
        positions.append(matchPos[0:2])
        matchPos = matchPos[2:]
    print positions
    print complement
from pyparsing import *
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
#define ENTRY element
entry = ft_entry + OneOrMore(Word(alphanums+"_."))("entry") + LineEnd().suppress()
#define NAME element
name_text = SkipTo(oneOf("DEFINITION ORTHOLOGY PATHWAY CLASS POSITION MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [[feat.lstrip().rstrip() for feat in feature[0].lstrip().rstrip().split(",")]])
name = Optional(ft_name + name_text("name"))
#define DEFINITION element
definition_text = SkipTo(ft_orthology | ft_pathway | ft_class | ft_position | ft_motif | ft_dblinks | ft_structure | ft_aaseq | ft_ntseq)
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
#position = Optional(ft_position + Word(nums)("pos_from") + Literal("..") + Word(nums)("pos_to") + LineEnd().suppress())
position = Optional(ft_position + SkipTo(oneOf("MOTIF DBLINKS STRUCTURE AASEQ NTSEQ")))("positions")
motifs = Optional(ft_motif + SkipTo(oneOf("DBLINKS STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [feat.lstrip().rstrip().split() for feat in feature[0].split("\n") if feat]))("motifs")
#define DBLINKS element
dblinks_text = SkipTo(oneOf("STRUCTURE AASEQ NTSEQ")).setParseAction(lambda feature: [[el.lstrip().rstrip().split() for el in feat.lstrip().rstrip().split(":", 1) if el.lstrip().rstrip().split()] for feat in feature[0].split("\n") if feat])
dblinks = Optional(ft_dblinks + dblinks_text)("dblinks")
structure = Optional(ft_structure + SkipTo(oneOf("AASEQ NTSEQ")).setParseAction(lambda feature: [feat.lstrip().rstrip().split() for feat in feature[0].split("\n")]))("structures")
#define AASEQ element
aaseq = Optional(ft_aaseq + number + LineEnd().suppress() + SkipTo(Keyword("NTSEQ")).setParseAction(lambda features: "".join([feat.lstrip().rstrip() for feat in features[0].split("\n")])))("aaseq")
#define NTSEQ element
ntseq = Group(ft_ntseq + number + LineEnd().suppress() + SkipTo("///").setParseAction(lambda features: "".join([feat.lstrip().rstrip() for feat in features[0].split("\n")])))("ntseq")
genes = entry + name + definition + orthology + pathway + class_element + position + motifs + dblinks + structure + aaseq + ntseq
parsedResultList = genes.parseString(a)
for key in parsedResultList.keys():
    print parsedResultList[key]
print "##############"
print parsedResultList.entry
print parsedResultList.name
print parsedResultList.definition.rstrip()
print parsedResultList.ko_orthology
print parsedResultList.pathways
print parsedResultList.class_element
print parsedResultList.complement
extractPositions(parsedResultList.positions[0])
print parsedResultList.motifs
print parsedResultList.dblinks
print parsedResultList.structures
print parsedResultList.aaseq
print parsedResultList.ntseq
