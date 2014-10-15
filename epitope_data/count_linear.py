import sys, os
import numpy as np
from Bio import SeqIO

def count_nonlinear(infile, outfile):
    lines = list(SeqIO.parse(open(infile, 'r'), 'fasta'))
    
    counts = np.zeros(550)
    
    for line in lines:
        for site in enumerate(line.seq):
            if site[0] + 1 > 16 and site[1] != '-':
                counts[site[0] - 16] += 1
            
    out_file = open(outfile, 'w')
    out_file.write('site\tcount\n')
    for element in enumerate(counts):
        out_file.write(str(element[0] + 1) + '\t' + str(element[1]) + '\n')
    out_file.close
        
def make_edges(infile, outfile):
    lines = list(SeqIO.parse(open(infile, 'r'), 'fasta'))
    
    edges = []
    
    for line in lines:
        sites = []
        for site in enumerate(line.seq):
            if site[0] + 1 > 16 and site[1] != '-':
                sites.append(site[0] - 16)
                
        for site1 in range(len(sites) - 1):
            for site2 in range(site1 + 1, len(sites)):
                edges.append(str(sites[site1]) + '\t' + str(sites[site2]) + '\n')
            
    
    out_file = open(outfile, 'w')
    out_file.write('n1\tn2\n')
    for element in edges:
        out_file.write(element)
    out_file.close

def main():
    infile = sys.argv[1]
    outfile_counts = sys.argv[2]
    outfile_edges = sys.argv[3]
    count_nonlinear(infile, outfile_counts)
    make_edges(infile, outfile_edges)

main()