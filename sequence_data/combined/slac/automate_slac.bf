/* Automates FEL with the following input. Written 10/1/2014 AGM. */

BASEDIR = "/Users/austin/Desktop/slac/";
datafile="nucleotide.fasta";
output="run.log";
treefile="nucleotide.tree";
sites="sites.dat";

inputRedirect = {};
inputRedirect["01"]="Universal";             //Genetic code
inputRedirect["02"]="New Analysis";          //New analysis
inputRedirect["03"]=BASEDIR+datafile;        //Fasta file, full path              
inputRedirect["04"]="Default";               //Use Default
inputRedirect["05"]=BASEDIR+treefile;        //Tree
inputRedirect["06"]=BASEDIR+output;          //Output
inputRedirect["07"]="Estimate + CI";              
inputRedirect["08"]="Single Ancestor Counting";
inputRedirect["09"]="Tips vs Internals";
inputRedirect["10"]="Count";  
inputRedirect["11"]=BASEDIR+sites; 

ExecuteAFile ("QuickSelectionDetection.bf", inputRedirect);
