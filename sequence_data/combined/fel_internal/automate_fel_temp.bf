/* Automates FEL with the following input. Written 10/1/2014 AGM. */

BASEDIR = "change_path/";
datafile="nucleotide.fasta";
output="model.dat";
treefile="nucleotide.tree";
sites="sites.dat";

inputRedirect = {};
inputRedirect["01"]="Universal";             //Genetic code
inputRedirect["02"]="New Analysis";          //New analysis
inputRedirect["03"]=BASEDIR+datafile;        //Fasta file, full path
inputRedirect["04"]="Default";               //Use HKY85 and MG94xHKY85
inputRedirect["05"]=BASEDIR+treefile;        //Tree
inputRedirect["06"]=BASEDIR+output;          //Output
inputRedirect["07"]="Estimate dN/dS only";   //Only estimate w
inputRedirect["08"]="Two rate FEL";          //Estimate one rate
inputRedirect["09"]="0.01";                   //Significance level
inputRedirect["10"]="Internal Only";                   //What parts of the tree
inputRedirect["11"]=BASEDIR+sites;           //sitewise output

ExecuteAFile ("QuickSelectionDetection.bf", inputRedirect);
