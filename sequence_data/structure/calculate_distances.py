#!/usr/bin/python

from Bio.PDB import *

def getResNumber(residue):
	id = residue.get_id()
	assert(id[0]==' ') # just make sure we're working with a properly numbered PDB file
	assert(id[2]==' ')
	return id[1]


pdbfile='4fnk_monomer_renumbered.pdb'

parser=PDBParser()

structure=parser.get_structure('', pdbfile)
model=structure[0]

distances = []

for chain in model:
    for r1 in chain:
        for r2 in chain:
            if getResNumber(r2) != getResNumber(r1):
                distances.append([getResNumber(r1), getResNumber(r2), float(r1['CA']-r2['CA'])])

this_map = {}
for i in range(1, 318):
    this_map[i+8+16] = i
for i in range(319, 490):
    this_map[i+8+16+3] = i

epitopes = list(open('../../../manuscript/numbering_table.csv', 'r').readlines())
sites = []
for i in enumerate(epitopes):
    split_line = i[1].split(',')
    if i[0] == 0:
        continue
    if split_line[15] != 'NA':
        if int(split_line[15]) > 0 and int(split_line[3]) in this_map.keys():
            sites.append(this_map[int(split_line[3])])
print(sites)

this_map = {}
for i in range(1, 318):
    this_map[i] = i+8+16
for i in range(319, 490):
    this_map[i] = i+8+16+3

out_file1 = open('4fnk_monomer_distances.dat', 'w')
out_file1.write('res1,res2,distance\n')
out_file2 = open('4fnk_monomer.edges', 'w')
out_file2.write('n1,n2\n')
for distance in distances:
    out_file1.write(str(distance[0]) + ',' + str(distance[1]) + ',' + str(distance[2]) + '\n')
    if distance[0] in sites and distance[1] in sites and distance[2] < 10:# and (distance[0] - distance[1])**2 > 1:
        out_file2.write(str(this_map[distance[0]]) + ',' + str(this_map[distance[1]]) + '\n')
out_file1.close()
out_file2.close()
