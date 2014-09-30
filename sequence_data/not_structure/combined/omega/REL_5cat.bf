/*This batch file uses Sergei`s function to build a likelihood fxn etc.*/


global k;
global t;
global w;

Cat_Num=5;
Cat_Name="w";

BuildCategory(Cat_Num, Cat_Name);

/* Include relevant functions */
#include "gy94_header.ibf";
#include "REL.mdl";

//OPTIMIZATION_PRECISION = 0.000001;
LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
//MAXIMUM_ITERATIONS_PER_VARIABLE = 100000;

/* Read in the data */

DataSet 		raw_data = ReadDataFile ("ha.nuc");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter (raw_data,3,"", "","TAA,TAG,TGA");

/* Collect observed frequencies into vectors */
HarvestFrequencies(dataFreq,filt_data,3,1,1);


/* Get the codon frequencies from the nucleotide frequency vectors */
codonFreq = BuildCodonFrequencies(dataFreq);

/* Define the model and tree */
Model GY94 = (GY94RateMatrix, codonFreq, 1);
UseModel(GY94);
Tree    Tree01 = DATAFILE_TREE;


LikelihoodFunction  LikFnPart = (filt_data, Tree01);

Optimize (paramValues, LikFnPart);

fprintf (stdout, LikFnPart);

ConstructCategoryMatrix(catmat,LikFnPart,COMPLETE);
ConstructCategoryMatrix (priors, LikFnPart, WEIGHTS);

fprintf(stdout, catmat);
fprintf(stdout, priors);


/*Compute posterior probabilities*/

result = {Rows(catmat),Columns(catmat)}; //rows are categories and columns are sites
for( i=0; i<Rows(catmat); i=i+1)
{
  for(j=0; j<Columns(catmat); j=j+1)
  {
    //fprintf(stdout,"\n",catmat[i][j], " ", priors[i][0],"\n");
    result[i][j] = catmat[i][j];
  }
}



/*Normalize the probabilities*/

for( j=0; j<Columns(catmat); j=j+1)
{
  sum = 0;
  for( i=0; i<Rows(catmat); i=i+1)
  {
    sum = sum+result[i][j];
  }
  for( i=0; i<Rows(catmat); i=i+1)
  {
    result[i][j] = result[i][j]/sum;
  }
}

fprintf(stdout,result); 


