pdb = '4fnk_monomer'
#cmd.load(pdb + '.pdb')

cmd.hide('all')
cmd.show('cartoon')
cmd.cartoon('tube')
cmd.set('cartoon_tube_radius', 0.9)

# open the file of new values (just 1 column of numbers, one for each alpha carbon)
inFile = open("alltree/fulltree.correlations", 'r')
 
# create the global, stored array
stored = []
 
# read the new B factors from file
for line in inFile.readlines(): stored.append( float(line) )
 
max_b = max(stored)
min_b = min(stored)

# close the input file
inFile.close()
 
# clear out the old B Factors
cmd.alter("%s and n. CA"%pdb, "b=0.0")
 
# update the B Factors with new properties
cmd.alter("%s and n. CA"%pdb, "b=stored.pop(0)")
 
# color the protein based on the new B Factors of the alpha carbons
cmd.spectrum("b", "rainbow", "%s and n. CA"%pdb, minimum=min_b, maximum=max_b)

cmd.ray("775", "2400")
#cmd.png("images/%s.png"%pdb[0:4])