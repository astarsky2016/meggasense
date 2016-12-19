'''
Created on Dec 12, 2011

@author: astar
'''
from persistent import Persistent
from pyparsing import Word, Group, Suppress, OneOrMore, alphas, nums, alphanums, Keyword, ParserElement, \
LineEnd, SkipTo, ZeroOrMore, Optional, oneOf, Literal
import pyparsing
from collections import defaultdict
ParserElement.setDefaultWhitespaceChars(" \t")
suppressElement1 = Literal("<b>").suppress()
suppressElement2 = Literal("</b>").suppress()
ko_number = Word('K', alphanums)("ko_number")
tagA = Keyword("A").suppress()
tagB = Keyword("B").suppress()
tagC = Keyword("C").suppress()
tagD = Keyword("D").suppress()
lineA = tagA + suppressElement1 + OneOrMore(Word(alphas)).setParseAction(lambda features: " ".join(features)) + suppressElement2 + LineEnd().suppress()
lineB = tagB + suppressElement1 + OneOrMore(Word(alphas)).setParseAction(lambda features: " ".join(features)) + suppressElement2 + LineEnd().suppress()
lineC = tagC + SkipTo(LineEnd()) + LineEnd().suppress()
lineD = tagD + ko_number + SkipTo(LineEnd()).suppress() + LineEnd().suppress()
ko_brite_parse_expr = Optional(SkipTo("A<b>")).suppress() + lineA("A") + Group(OneOrMore(tagB + LineEnd().suppress() + 
                                     lineB + Group(OneOrMore(Group(OneOrMore(lineC)) + Group(ZeroOrMore(lineD))))))("hierarchy")
class koBriteEntry(Persistent):
    def __init__(self, brite):
        self.brite = brite
        self.__parseBrite()
        self.hierarchy = defaultdict(lambda: defaultdict(defaultdict))
    def __parseBrite(self):
        parsedResultList = ko_brite_parse_expr.parseString(self.brite)
        self.A = parsedResultList.A[0]
        self.rest = parsedResultList.hierarchy
                  
if __name__ == "__main__":
    import logging, sys, re, traceback
    logging.basicConfig(level=logging.WARNING)
    log = logging.getLogger('/media/Transcend/Working_KEGG/KEGG_core/log')
    f = open('/media/Transcend/Working_KEGG/KEGG_core/ko00001.keg', 'r')
    testString = f.read().split('#\n')
    for string in testString:
        try:
            x = koBriteEntry(string)
        except pyparsing.ParseException:
            traceback.print_exc(file=sys.stdout)
            sys.exit(1)