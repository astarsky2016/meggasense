'''
Created on Aug 1, 2012

@author: astar
'''
import subprocess, os, glob
outDir = "/home/users/astar/TAGS/glimmer-mg_02/scripts/matis_za_analizu"
def annotate(query, out):
    query = os.path.join(outDir, query)
    out = os.path.join(outDir, out)
    annotate_cmd = """/home/users/astar/hmmer-3.0-linux-intel-x86_64/PfamScan/pfam_scan.pl -fasta %s -dir 
    /home/users/astar/TAGS/hmmerDB/kegg58 -custom kegg_58 -cpu 100 -outfile %s""" % (query, out)
    p = subprocess.Popen(annotate_cmd, shell=True)
    finished = p.wait()
def main():
    for fasta in glob.glob("*.fasta"):
        outFile = os.path.splitext(fasta)[0] + "_reads.hmmer_kegg58.out"
        annotate(fasta, outFile)