'''
Created on Mar 20, 2012

@author: astar
'''
import os, sys, time
from subprocess import Popen
assembly_dir = '/home/users/astar/TAGS/assembly/sample_A1_67_kegg_assembly'
def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
class MongoPar:
    def __init__(self):
        dirList = [os.path.join(assembly_dir,d) for d in os.listdir(assembly_dir) if os.path.isdir(os.path.join(assembly_dir,d))]
        self.max_proc = 112
        self.cmds = []
        counter = 1
        for chunk in chunks(dirList, 112):
            name = os.path.join('/home/users/astar/tmp', "file" + str(counter))
            file = open(name, "w")
            file.write("\n".join(chunk))
            file.close()
            self.cmds.append("python /home/users/astar/mrcina/mongoChunkLoader.py -i %s" % name)
            counter += 1
        self.__exec_par()
   
    def __exec_par(self):
        total = len(self.cmds)
        finished = 0
        running = 0
        p = []
    
        while finished + running < total:
            # launch jobs up to max
            while running < self.max_proc and finished+running < total:
                p.append(Popen(self.cmds[finished+running], shell=True))
                #print 'Running %d' % p[running].pid
                running += 1
    
            # are any jobs finished
            new_p = []
            for i in range(len(p)):
                if p[i].poll() != None:
                    running -= 1
                    finished += 1
                else:
                    new_p.append(p[i])
    
            # if none finished, sleep
            if len(new_p) == len(p):
                time.sleep(1)
            p = new_p
    
        # wait for all to finish
        for i in range(len(p)):
            p[i].wait()
       
if __name__ == "__main__":
    MongoPar()
        
         
    
