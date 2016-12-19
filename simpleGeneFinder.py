from Bio import SeqIO
import sys, getopt, re
from collections import defaultdict
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
def findGenes(inputfile, outputfile):
    bpNumber = 54
    SD_buffer_bp = 15
    starters = ["ATG", "GTG", "TTG", "CAT", "CAC", "CAA"]
    stoppers = ["TGA", "TAG","TAA", "TCA", "CTA", "TTA"]
    SD_motif = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
    SD_motif.add_instance(Seq("GGAG",SD_motif.alphabet))
    SD_motif.add_instance(Seq("GAGG",SD_motif.alphabet))
    SD_motif.add_instance(Seq("AGGA",SD_motif.alphabet))
    startCodonNumber = defaultdict(int)
    stopCodonNumber = defaultdict(int)
    SDNumber = defaultdict(int)
    for read in SeqIO.parse(inputfile, "fasta"):
        counter = 0
        for startCodon in starters:
            hits = [match.start() for match in re.finditer(re.escape(startCodon), read.seq.tostring())]
            for hit in hits:
                extend = 0
                start = 0
                if counter < 3:
                    extend = hit + bpNumber
                    start = hit
                else:
                    if hit > bpNumber: start = hit - bpNumber 
                    extend = hit
                startCodonNumber[read.seq[start:extend].tostring()] += 1
                SD_start = hit - SD_buffer_bp
                if SD_start < 0: SD_start = 0
                SD_stop = hit + SD_buffer_bp + 3
                SD_buffer = read.seq[SD_start:SD_stop]
                if [(pos, seq) for pos,seq in SD_motif.search_instances(SD_buffer)]:
                    print "found SD gene"
                    SDNumber[read.seq[start:extend].tostring()] += 1
            counter += 1
        counter = 0
        for stopCodon in stoppers:
            hits = [match.start() for match in re.finditer(re.escape(stopCodon), read.seq.tostring())]
            for hit in hits:
                extend = 0
                start = 0
                if counter < 3:
                    if hit > bpNumber:
                        start = hit - bpNumber
                        extend = hit 
                else:
                    start = hit + 3
                    extend = start + bpNumber
                stopCodonNumber[read.seq[start:extend].tostring()] += 1
            counter += 1
        
    f = open(outputfile, "w")
    f.write("number of start codons is: %i\n" % sum([value for value in startCodonNumber.values()]))
    f.write("number of stop codons is: %i\n" % sum([value for value in stopCodonNumber.values()]))
    f.write("number of SD start codons is: %i\n" % sum([value for value in SDNumber.values()]))
    f.close()
                    
                    

    
def main():
#    inputfile = ''
#    outputfile = ''
#    try:
#        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
#    except getopt.GetoptError:
#        print 'test.py -i <inputfile> -o <outputfile>'
#        sys.exit(2)
#    for opt, arg in opts:
#        if opt == '-h':
#            print 'test.py -i <inputfile> -o <outputfile>'
#            sys.exit()
#        elif opt in ("-i", "--ifile"):
#            inputfile = arg
#        elif opt in ("-o", "--ofile"):
#            outputfile = arg
    inputfile = "HI7IROS01_02.fasta"
    outputfile = "summary.out"
    print 'Input file is: "%s"' % inputfile
    print 'Output file is: "%s"' % outputfile
    print "processing ..."
    
    findGenes(inputfile, outputfile)
if __name__ == "__main__":
    main()
    #for seq in SeqIO.parse("/media/Transcend/Epicoccum/DBLargeContigs.fas", "fasta"):
