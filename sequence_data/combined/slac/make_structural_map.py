###################################################################################
## This script was written by Austin Meyer in the lab of Dr. Claus Wilke at the  ##
## University of Texas at Austin.  For any questions or concerns you can email   ##
## Austin at austin.g.meyer@gmail.com.                                           ##
###################################################################################

###################################################################################
## This script takes command line arguments in the following form                ##
## python translate_align_mapping.py fasta_file pdb_structure_file output_type   ##
##                                                                               ##
## Here is an example run to output the aligned and reverse translated sequences ##
## python translate_align_mapping.py fasta.fasta structure.pdb seqs              ##
##                                                                               ##
## Here is an example run to return the sequence-structure map                   ##
## python translate_align_mapping.py fasta.fasta structure.pdb map               ##
###################################################################################

######################################################################
## Some of these includes are probably unecessary.                  ##
######################################################################
import sys, subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

######################################################################
## This function accepts the full unfettered sequence file full pdb ##
## file, and calls the above translate_align function to do its     ##
## job.  Then, it calls the rsa determining function above and the  ##
## header function.  It combines the information to output a single ##
## file with the codon, residue, and rsa.  There are a bunch of     ##
## other columns just to satisfy the format of the script Mike has  ##
## already written.                                                 ##
######################################################################

def main():

  args =  sys.argv
  alignment =  args[1]
  pdb      =  args[2]

  aas = get_aa_fromPDB( pdb )

  make_combined_alignment(alignment, aas)
    
  return 0
  
######################################################################
## This function returns the structure attribute of a PDB file      ##
## from the PDB parser module.                                      ##
######################################################################

def parsePDBStructure( pdb_id ):
    parser = PDBParser()
    structure = parser.get_structure('test_rsa', pdb_id)
    return structure

######################################################################
## This function allows me to renumber the residues in a chain      ##
## to fix a particular odd numbering problem in Neuraminidase       ##
######################################################################

def get_aa_fromPDB( pdb_id ):
  structure = parsePDBStructure( pdb_id )
  polypeptides = ""
  ppb=PPBuilder()
  for pp in ppb.build_peptides(structure):
    polypeptides += pp.get_sequence()
  return polypeptides

def make_combined_alignment(alignment, aas):
  aa_file = open('aas.tmp', 'w')
  aa_file.write('>ref_seq\n' + str(aas) + '\n')
  aa_file.close()

  subprocess.call('mafft --add aas.tmp ' + str(alignment) + ' > new_alignment.tmp', shell=True) 

  lines = list(SeqIO.parse(open('new_alignment.tmp', 'r'), 'fasta'))

  for line in lines:
    if line.id == 'ref_seq':
      outfile = open('structural_map.txt', 'w')
      for site in line.seq:
        outfile.write(str(site) + '\n')
      outfile.close()
      break

## Execute the main function

if __name__ == "__main__":
    main()

