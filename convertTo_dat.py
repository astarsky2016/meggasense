'''
Created on Jan 26, 2012

@author: astar
'''
import optparse, os
def convert(hmmLib):
    header = "# STOCKHOLM 1.0"
    sep = "//"
    f = open(hmmLib, 'r')
    (base,ext) = os.path.splitext(options.hmm_file)
    dat = open('%s.dat' % base, 'w')
    try:
        id = ''
        ml = ''
        records = f.read().split("//")
        for hmm in records:
            for line in hmm.split("\n"):
                if line.startswith("NAME"):
                    id = "#=GF ID   " + line.split(" ")[2]
                if line.startswith("LENG"):
                    ml = "#=GF ML   " + line.split(" ")[2]
            dat.write("\n".join([header, id, ml, sep]))
            dat.write("\n")
    finally:
        f.close()
if __name__ == "__main__":
    usage = 'usage: %prog [options] -i <hmm file>'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', dest='hmm_file', help='hmm file')
    (options,args) = parser.parse_args()
    if not options.hmm_file:
        parser.error('Must provide hmm file with -i')    
    (base,ext) = os.path.splitext(options.hmm_file)
    out = open('%s.dat' % base, 'w')
    convert(options.hmm_file)