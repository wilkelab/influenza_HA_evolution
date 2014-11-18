import numpy as np
from Bio.PDB import *

def calc_distances_aas(model):
    count = 1
    sites = []
    for a in model.get_residues():
        if a.get_resname() == 'SIA':
            continue
        
        distances = []
        for b in model.get_residues():
            if b.get_resname() != 'SIA':
                distance = b['CA'] - a['CA']
                distances.append(distance)
        sites.append(distances)
        print(count)
        count += 1

    out_file = open('distances.dat', 'w')
    for a in sites:
        for b in a:
          if b == a[-1]:
                out_file.write(str(b))
          else:
                out_file.write(str(b) + ',')
        out_file.write('\n')
    out_file.close()
                
def calc_distances_sia():

    location = np.zeros(3)

    for a in model.get_residues():
        if a.get_resname() == "SIA":
            for a1 in a:
                location += a1.get_coord()
              
    for a in model.get_residues():
        if a.get_resname() == "SIA":
            for b in model.get_residues():
                if b.get_resname() != 'SIA':
                    distances.write(str(np.sqrt(sum((location - b['CA'].get_coord())**2))) + '\n')
                
def main():
    parser = PDBParser()
    structure = parser.get_structure('temp', '4fnk_SA_bound.pdb')
    calc_distances_aas(structure[0])
    
main()
