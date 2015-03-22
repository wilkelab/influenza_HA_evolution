/* Automates FEL with the following input. Written 10/1/2014 AGM. */

BASEDIR = "/home/austin/Work/influenza_HA_evolution/sequence_data/not_structure/change/slac/";
datafile="nucleotide.fasta";
output="run.log";
treefile="nucleotide.tree";
sites="sites.dat";

inputRedirect = {};
inputRedirect["01"]="Universal";             //Genetic code
inputRedirect["02"]="New Analysis";          //New analysis
inputRedirect["03"]=BASEDIR+datafile;        //Fasta file, full path
inputRedirect["04"]="Custom";               
inputRedirect["05"]="Default";               //Use HKY85 and MG94xHKY85
inputRedirect["05"]=BASEDIR+treefile;        //Tree
inputRedirect["06"]=BASEDIR+output;          //Output
inputRedirect["07"]="Single Ancestor Counting";   //Only estimate w
inputRedirect["08"]="Tips vs Internals";          //Estimate one rate
inputRedirect["09"]="Count";                   //Significance level
inputRedirect["10"]=BASEDIR+sites;           //sitewise output

ExecuteAFile ("QuickSelectionDetection.bf", inputRedirect);
