'''
Created on Aug 2, 2012

@author: astar
'''
import subprocess, os, glob
outDir = "/home/users/astar/TAGS/glimmer-mg_02/scripts/matis_za_analizu"
def runGlimmer(fastaFile, project):
    fasta = os.path.join(outDir, fastaFile)
    annotate_cmd = "/home/users/astar/TAGS/glimmer_mg/scripts/glimmer_mg.py --indel -p 60 %s %s" % (fasta, project)
    p = subprocess.Popen(annotate_cmd, shell=True)
    finished = p.wait()
def main():
    for fasta in glob.glob("*.fasta"):
        print "working on %s ..." % fasta
        project = os.path.splitext(fasta)[0]
        runGlimmer(fasta, project)
    print "finished!"
if __name__ == "__main__":
    main()