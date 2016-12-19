'''
Created on Jan 23, 2012

@author: astar
'''
from optparse import OptionParser
import subprocess, os, sys
toolBox_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
sff_extract = '/home/astar/BioApps/sff_extract_0_2_11.py'
####################################################################
def main():
    parser = OptionParser(usage="usage: %prog [options] <sff_file>",
                          version="%prog 1.0")
    parser.add_option("-p", "--paired",
                      action="store_true",
                      dest="paired_end",
                      default=False,
                      help="to extract the paired-end sequences, needs SSAHA2 in path nad linker fasta file")
    parser.add_option("-c", "--clip",
                      action="store_true", # optional because action defaults to "store"
                      dest="clip",
                      default=False,
                      help="clip the fragments recommended by the 454 software")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("must provide name of the sff file to convert!")
#    else:
#        params = dict(c='', l='')
#        sff = ''
#        linker = ''
#        if len(args) >= 2:
#            if options.clip: params['c'] = '-c'
#            sff = ' '.join(args)
#        else:
#            sff = args[0]
#            if options.paired_end: 
#                params['l'] = '-l'
#                linker = args[0]
#                sff = args[1]
#        comm = 'python %s %s %s %s %s' % (sff_extract, params['c'], params['l'], linker, sff)
#        print comm
    print options
    print args

if __name__ == '__main__':
    main()
