import sys, os
import numpy as np

def count_nonlinear(infile, outfile):
    lines = list(open(infile, 'r').readlines())
    
    counts = np.zeros(550)
    
    for line in lines:
        split_line = ((line.strip()).split('\t')[1]).split(',')
        for element in split_line:
            counts[int(element[1:len(element)]) - 1] += 1
            
    out_file = open(outfile, 'w')
    out_file.write('site\tcount\n')
    for element in enumerate(counts):
        out_file.write(str(element[0] + 1) + '\t' + str(element[1]) + '\n')
    out_file.close
        
def make_edges(infile, outfile):
    lines = list(open(infile, 'r').readlines())
    
    edges = []
    
    for line in lines:
        split_line = ((line.strip()).split('\t')[1]).split(',')
        for element1 in range(len(split_line) - 1):
            for element2 in range(element1 + 1, len(split_line)):
                edges.append(split_line[element1][1:len(split_line[element1])] + '\t' + split_line[element2][1:len(split_line[element2])] + '\n')
    
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
