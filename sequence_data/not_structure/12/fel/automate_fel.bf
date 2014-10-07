/* Automates FEL with the following input. Written 10/1/2014 AGM. */

BASEDIR = "/Users/austin/Google Drive/Data/influenza_HA_evolution/sequence_data/not_structure/12/";
datafile="nucleotide.fasta";
output="run.log";
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
inputRedirect["08"]="One rate FEL";          //Estimate one rate
inputRedirect["09"]="0.1";                   //Significance level
inputRedirect["10"]=BASEDIR+sites;           //sitewise output

ExecuteAFile ("QuickSelectionDetection.bf", inputRedirect);
