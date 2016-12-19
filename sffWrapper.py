'''
Created on Mar 22, 2012

@author: astar
'''
from collections import defaultdict
from pyparsing import Keyword
import re, subprocess, os
from Bio import SeqIO
class Chdir:         
    def __init__(self, newPath):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)
    def __del__(self):
        os.chdir(self.savedPath)
class AssemblyPipe:  
    def __init__(self, binDir, hmmerData, sffFile, assembly_tmp_dir):         
        self.newblerDir = binDir
        self.assembly_dir = assembly_tmp_dir
        self.sff = sffFile
        self.kos = hmmerData   
        self.__dictinvert()
        self.__extractAndAssemble()
    def __dictinvert(self):
        inv = defaultdict(list)
        for key, val in self.kos.items():
            for el in val: 
                inv[el[2]].append(key)
        self.kos = inv
    def __extractAndAssemble(self):
        fnull = open(os.devnull, 'w')  
        changeDir = Chdir(self.assembly_dir)
        cmd = os.path.join(self.newblerDir, 'runAssembly')
        readsff = SeqIO.index(self.sff, "sff")
        for func, readList in self.kos.items():
            f = open("selection.sff", "wb")
            writer = SeqIO.SffIO.SffWriter(f)
            writer.write_file((readsff[read.split("_")[0]] for read in readList))
            f.close()
            cmd_expr = "%s -o %s selection.sff" % (cmd, func)
            p = subprocess.Popen(cmd_expr, shell=True, stdout = fnull)
            os.waitpid(p.pid, 0)
        del changeDir
        fnull.close
        print "assembly completed!"
  
    
    