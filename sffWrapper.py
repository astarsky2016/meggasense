'''
Created on Mar 22, 2012

@author: astar
'''
import csv, os, tempfile
from Bio import SeqIO
from collections import defaultdict
from HMMER.parseHmmer_tab import HmmerTabParser
from pyparsing import Keyword
import re, subprocess, sys
class Chdir:         
      def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)

      def __del__( self ):
        os.chdir( self.savedPath )
   
class AssemblyPipe:  
    def __init__(self, binDir, annotationFile, sffFile):         
        self.newblerDir = binDir
        self.data = annotationFile
        self.sff = sffFile
        self.kos = HmmerTabParser(self.data, e_threshold = 10, overlap_limit = 0.5, infoLine = re.compile("\S"), startLine = Keyword("# <seq id>"))   
        self._dictinvert()
        self._extractAndAssemble()
    def _dictinvert(self):
        inv = defaultdict(list)
        for k, v in self.kos.hmmerAnnotation.items():
            for ko in v: inv[ko[2]].append(k)
        self.kos = inv
    def _extractAndAssemble(self):
        #assembly_tmp_dir = tempfile.mkdtemp()
        fnull = open(os.devnull, 'w')  
        assembly_tmp_dir = "/home/users/astar/TAGS/assembly/sample_A1_67_kegg_assembly"
        changeDir = Chdir(assembly_tmp_dir)
        os.symlink(self.sff, os.path.join(assembly_tmp_dir, os.path.basename(self.sff)))
        readsff = SeqIO.index(self.sff, "sff")
        cmd = os.path.join(self.newblerDir, 'runAssembly')
        for func, readList in self.kos.items():
            f = open("selection.sff", "wb")
            writer = SeqIO.SffIO.SffWriter(f)
            writer.write_file((readsff[read] for read in readList))
            f.close()
            cmd_expr = "%s -o %s selection.sff" % (cmd, func)
            p = subprocess.Popen(cmd_expr, shell=True, stdout = fnull)
            os.waitpid(p.pid, 0)
        del changeDir
        fnull.close
        print "assembly completed!"
    #os.removedirs(assembly_tmp_dir)
  
    
    